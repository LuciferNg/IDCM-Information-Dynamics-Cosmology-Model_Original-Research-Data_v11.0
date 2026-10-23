#!/usr/bin/env python3
"""
M4 R3 — lattice 高精度值重判 (2026-10-24)
==========================================
動機（SESSION_START §四 R3）：/tmp scripts 用 PDG 引用值 m_u=2.16±0.07、
m_d=4.67±0.48、m_s=93.4±2.5；其中 m_d 嘅 ±0.48（10.3%）令
  s/d=20（dev 0.000%）、span≈12（dev 0.042%）兩個「命中」全部落入誤差帶 ⟹ 判後驗。
本輪換 FLAG 2024 lattice 值（誤差細 1 個數量級）重判同一批命中。

數據來源（實檔核查，非二手引用）：
  FLAG Review 2024, arXiv:2411.04268（PDF 已下載 cache/scratch/m4_flag/flag2024.pdf,
  md5 3eedc5a7a51d3cbda1c460e00f1e2a9c）:
    Table 1 / eq.(34): Nf=2+1+1  mud=3.427(51), ms=93.46(58), ms/mud=27.227(81)
    Table 1 / eq.(32): Nf=2+1    mud=3.387(39), ms=92.4(1.0), ms/mud=27.42(12)
    Table 12:          Nf=2+1+1  mu=2.14(8), md=4.70(5), mu/md=0.465(24)
                       Nf=2+1    mu=2.27(9), md=4.67(9), mu/md=0.485(19)
  重夸克維持 PDG 2024（MS̄）：mc(mc)=1.27, mb(mb)=4.18, mt(pole)=172.69 GeV
  ※ FLAG 只給 mc(3GeV)=0.989(10)，尺度唔同 ⟹ 唔混用（跨尺度比較會污染）

測試（全部事前定義）：
  T1  c/s vs γ20/γ5      （/tmp 命中 dev 0.0424%）
  T2  s/d vs hop_max 20   （/tmp 命中 dev 0.000%）
  T3  span_Q/span_ν vs 12 （/tmp 命中 dev 0.042%）
  T4  m_d/m_u vs 結構常數池（新增：lattice 誤差令此比首次可用）
每個測試報：值、σ（誤差傳播）、dev、dev/σ、density PN（同 09 protocol）
判準：dev 在 <1σ 內且 density PN>5% ⟹ 只算「一致」，唔算閉合
"""
import json
import numpy as np

L = 88
n = np.arange(1, 45, dtype=float)
gamma = 4.0*np.sin(n*np.pi/L)**2
pool = (gamma[:, None]/gamma[None, :]).ravel()
lp = np.sort(np.log(pool))

def dev_of(t, k=4):
    i = np.searchsorted(lp, np.log(t))
    lo = max(0, i-k); hi = min(lp.size, i+k)
    cand = np.exp(lp[lo:hi])
    dv = np.abs(cand/t - 1.0)*100.0
    j = np.argmin(dv)
    return dv[j], cand[j]

def mc_family(thr, lo, hi, Nrat, NMC=200_000, seed=20261026):
    rng = np.random.default_rng(seed)
    hit = np.zeros(NMC, dtype=bool)
    for _ in range(Nrat):
        u = rng.uniform(np.log(lo), np.log(hi), NMC)
        t = np.exp(u)
        idx = np.searchsorted(lp, u)
        best = np.full(t.size, np.inf)
        for k in range(-3, 4):
            ii = np.clip(idx+k, 0, lp.size-1)
            best = np.minimum(best, np.abs(np.exp(lp[ii])/t - 1.0)*100.0)
        hit |= (best < thr)
    return float(np.mean(hit)*100.0)

# ── 質量集 ──
m_t, m_b, m_c = 172.69, 4.18, 1.27          # GeV
# 中微子（O257 R7 定案值，eV）
m1, m2, m3 = 0.020908, 0.022645, 0.053586

datasets = {
    "PDG2024(/tmp 用值)": dict(mu=0.00216, du=0.00007, md=0.00467, dd=0.00048, ms=0.0934, ds=0.0025),
    "FLAG24 Nf=2+1+1":   dict(mu=0.00214, du=0.00008, md=0.00470, dd=0.00005, ms=0.09346, ds=0.00058),
    "FLAG24 Nf=2+1":     dict(mu=0.00227, du=0.00009, md=0.00467, dd=0.00009, ms=0.0924, ds=0.0010),
}

print("="*80)
print("M4 R3: lattice 高精度重判（FLAG 2024 Nf=2+1+1 / 2+1 vs PDG）")
print("="*80)
print(f"{'dataset':<20} {'mu(MeV)':>12} {'md(MeV)':>12} {'ms(MeV)':>12}  相對誤差 u/d/s")
for k, d in datasets.items():
    print(f"{k:<20} {d['mu']*1e3:8.3f}±{d['du']*1e3:5.3f} {d['md']*1e3:8.3f}±{d['dd']*1e3:5.3f} "
          f"{d['ms']*1e3:8.3f}±{d['ds']*1e3:5.3f}  "
          f"{d['du']/d['mu']*100:4.1f}%/{d['dd']/d['md']*100:4.1f}%/{d['ds']/d['ms']*100:4.1f}%")

g20g5 = gamma[19]/gamma[4]
print(f"\n事前結構數: γ20/γ5 = {g20g5:.9f}  (γ20={gamma[19]:.9f}, γ5={gamma[4]:.9f});  hop_max = 12;  s/d 目標 = 20")

# ── 誤差傳播（獨立高斯，MC 2e5）──
def prop(fn, d, NMC=200_000, seed=777):
    rng = np.random.default_rng(seed)
    mu = rng.normal(d['mu'], d['du'], NMC)
    md = rng.normal(d['md'], d['dd'], NMC)
    ms = rng.normal(d['ms'], d['ds'], NMC)
    v = fn(mu, md, ms)
    return float(np.mean(v)), float(np.std(v))

tests = {}
print("\n" + "-"*80)
print("T1  c/s vs γ20/γ5")
print("-"*80)
for k, d in datasets.items():
    val = m_c/d['ms']
    sig = val*(d['ds']/d['ms'])
    dev = abs(g20g5/val - 1)*100
    print(f"  {k:<20} c/s = {val:.6f} ± {sig:.6f}  γ20/γ5 = {g20g5:.6f}  dev {dev:.4f}% = {dev/(sig/val*100):.2f}σ")
    tests[f"T1 c/s {k}"] = dict(value=val, sigma=sig, target=g20g5, dev_pct=dev, nsigma=dev/(sig/val*100))

print("\n" + "-"*80)
print("T2  s/d vs hop_max = 20（/tmp 記錄 dev 0.000%，因 m_d 誤差 10.3%）")
print("-"*80)
for k, d in datasets.items():
    val, sig = prop(lambda mu, md, ms: ms/md, d)
    dev = abs(val/20.0 - 1)*100
    print(f"  {k:<20} s/d = {val:.5f} ± {sig:.5f}  dev vs 20 = {dev:.4f}% = {dev/(sig/val*100):.2f}σ  "
          f"1σ 帶 [{val-sig:.3f}, {val+sig:.3f}]")
    tests[f"T2 s/d {k}"] = dict(value=val, sigma=sig, target=20.0, dev_pct=dev, nsigma=dev/(sig/val*100))

print("\n" + "-"*80)
print("T3  span_Q/span_ν vs 12 = hop_max（span_Q = log10(m_t/m_u), span_ν = log10(m3/m1)）")
print("-"*80)
span_nu = np.log10(m3/m1)
print(f"  span_ν = log10({m3}/{m1}) = {span_nu:.6f}（O257 R7 值；m_ν 誤差 <0.1%）")
for k, d in datasets.items():
    spanQ, sigQ = prop(lambda mu, md, ms: np.log10(m_t/mu), d)  # 同單位 GeV，比值無量綱
    R = spanQ/span_nu
    sigR = sigQ/span_nu
    dev = abs(R/12.0 - 1)*100
    print(f"  {k:<20} span_Q = {spanQ:.6f} ± {sigQ:.6f}  R = {R:.6f} ± {sigR:.6f}  "
          f"dev vs 12 = {dev:.4f}% = {dev/(sigR/R*100):.2f}σ  1σ [{R-sigR:.3f}, {R+sigR:.3f}]")
    tests[f"T3 span {k}"] = dict(value=R, sigma=sigR, target=12.0, dev_pct=dev, nsigma=dev/(sigR/R*100))

print("\n" + "-"*80)
print("T4  m_d/m_u vs 結構常數池（新增能力：lattice 令 m_d 誤差 1.1%）")
print("-"*80)
struct = {"L5/F5=11/5": 11/5, "L7/F7=29/13": 29/13, "B=4?": 4.0, "sqrt(F5)=√5": np.sqrt(5),
          "L5/B=11/4": 11/4, "F7/F5=13/5": 13/5, "N_gen²=9? ": 9.0, "L5/F7*? ": 11/13}
for k, d in datasets.items():
    val, sig = prop(lambda mu, md, ms: md/mu, d)
    print(f"  {k:<20} m_d/m_u = {val:.5f} ± {sig:.5f}")
    row = []
    for sname, sv in struct.items():
        dev = abs(val/sv - 1)*100
        if dev < 1.0:
            row.append(f"{sname} dev {dev:.3f}% ({dev/(sig/val*100):.2f}σ)")
    print("      候選 <1%: " + ("; ".join(row) if row else "冇"))
    tests[f"T4 md/mu {k}"] = dict(value=val, sigma=sig, candidates=row)

# T4 density 閘：事前族 = 結構整數比池 K={1,2,3,4,5,7,8,11,13,21,29,44,76,88}（p/q，1≤比≤3）
K = np.array([1, 2, 3, 4, 5, 7, 8, 11, 13, 21, 29, 44, 76, 88], dtype=float)
fam = np.array([p/q for p in K for q in K if 1.0 <= p/q <= 3.0])
print(f"\n  T4 density 族: 整數比池 {len(fam)} 值（域 1-3）")
rng_f = np.random.default_rng(20261027)
NMC_F = 1_000_000
tt = np.exp(rng_f.uniform(np.log(1.0), np.log(3.0), NMC_F))
for k in datasets:
    dev = tests[f"T4 md/mu {k}"].get("candidates") or []
    val = tests[f"T4 md/mu {k}"]["value"]
    near = fam[np.abs(fam/val - 1)*100 < 0.1]
    thr = 0.032 if any("L5" in s for s in dev) else None
    if thr is None:
        print(f"  {k:<20} 無 <1% 候選 → 免 density")
        continue
    dn = np.min(np.abs(fam[None, :]/tt[:, None] - 1)*100, axis=1)
    p1 = float(np.mean(dn < thr)*100)
    print(f"  {k:<20} m_d/m_u={val:.5f}: 族內 <0.1% 接近者 {len(near)} 個 "
          f"{np.round(near,4).tolist()} | 單比 P1(thr {thr}%) = {p1:.2f}%")
    tests[f"T4 md/mu {k}"]["family_near_01pct"] = near.tolist()
    tests[f"T4 md/mu {k}"]["family_P1_pct"] = p1

print("\n" + "-"*80)
print("Density 閘（09 protocol）對 T1/T2/T3 命中")
print("-"*80)
qd = datasets["FLAG24 Nf=2+1+1"]
qrat = {}
vals = {"u": qd['mu'], "d": qd['md'], "s": qd['ms'], "c": m_c, "b": m_b, "t": m_t}
qs = list(vals)
for i, a in enumerate(qs):
    for b in qs[i+1:]:
        qrat[f"{a}/{b}"] = vals[a]/vals[b]
lo, hi = min(qrat.values()), max(qrat.values())
for nm, key in [("T1 c/s (FLAG 2+1+1)", "T1 c/s FLAG24 Nf=2+1+1"),
                ("T2 s/d (FLAG 2+1+1)", "T2 s/d FLAG24 Nf=2+1+1")]:
    dev = tests[key]["dev_pct"]
    if dev < 1.0:
        pn = mc_family(dev, lo, hi, len(qrat))
        tests[key]["density_PN_pct"] = pn
        print(f"  {nm}: dev {dev:.4f}% → 同域同 N={len(qrat)} 家族級 PN = {pn:.2f}%  "
              + ("🔴 隨機級" if pn > 5 else "🟡"))
    else:
        print(f"  {nm}: dev {dev:.4f}% ≥1% → 免 MC")

print("\n" + "="*80)
print("判讀：dev<1σ 只代表「同觀測一致」；PN>5% ⟹ 池密度足以解釋，唔可升格。")
print("="*80)

with open("data/11_lattice_flag.json", "w", encoding="utf-8") as f:
    json.dump({"datasets": datasets, "tests": tests}, f, ensure_ascii=False, indent=1, default=str)
print("→ data/11_lattice_flag.json")
