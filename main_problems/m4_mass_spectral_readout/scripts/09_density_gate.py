#!/usr/bin/env python3
"""
M4 R1-corr — 統一 density 閘 (2026-10-24)
==========================================
問題（R1 重建時發現）：01-08 用咗三種唔同 MC null protocol：
  02: 家族級（15 比至少一命中），thr 0.1%，質量域 log10(2e-3..173) → 97.9%
  04: 單比級，thr 0.05%，質量域 log-均勻 e..τ → 16.8%
  05: 單比級，thr 0.1%，差比域 0.1-10 → 43.3%（記錄寫 8.7%，唔可復現）
⟹ 三個 sector 嘅「density 殺」數字互不可比，亦冇一條用「命中自身 dev」做 thr。

本腳本 = 統一 protocol（事前固定）：
  1. 池 P = {γ_a/γ_b | a,b ∈ 1..44}（1936 值，K_chain L=88 譜）
  2. 每次命中用**該命中自身 dev** 做 threshold（唔可以放寬）
  3. 域 = 該 sector 受測比值嘅 [min,max]（log-均勻抽樣）——即測試實際覆蓋嘅比值範圍
  4. 兩個判準都報：
       單比級 P1 = 一個隨機目標落入 thr 嘅機率
       家族級 PN = N 個獨立隨機目標至少一個命中（N = 該 sector 受測比數）
  5. NMC = 10^6（每一 thr 級）；向量化，log 空間最近鄰 + 鄰居精確 dev 校正
輸出：data/09_density_gate.json（機器讀）＋ stdout 表
"""
import json
import numpy as np

# ── K_chain 譜 ──
L = 88
n = np.arange(1, 45, dtype=float)
gamma = 4.0*np.sin(n*np.pi/L)**2
pool = (gamma[:, None]/gamma[None, :]).ravel()
lp = np.sort(np.log(pool))
print("="*78)
print("M4 R1-corr: 統一 density 閘（單一 protocol 重判全部 sector 命中）")
print("="*78)
print(f"池: γ_a/γ_b, a,b ∈ 1..44 → {pool.size} 值;  range [{pool.min():.6g}, {pool.max():.6g}]")

# ── 最近鄰 dev（log 空間定位 + 鄰居精確校正） ──
def dev_of(t, k=4):
    """回傳 target t 喺池中嘅最小相對 dev（%）同該池值。"""
    i = np.searchsorted(lp, np.log(t))
    lo = max(0, i-k); hi = min(lp.size, i+k)
    cand = np.exp(lp[lo:hi])
    dv = np.abs(cand/t - 1.0)*100.0
    j = np.argmin(dv)
    return dv[j], cand[j]

# 驗證最近鄰 vs 全池暴力（抽 200 個 target）
rng_v = np.random.default_rng(4242)
for _ in range(200):
    t = 10**rng_v.uniform(np.log10(0.01), np.log10(300))
    d1, _ = dev_of(t)
    d2 = np.min(np.abs(pool/t - 1.0)*100.0)
    assert abs(d1-d2) < 1e-9, (t, d1, d2)
print("最近鄰實作 vs 全池暴力: 200/200 一致 ✅")

def mc_single(thr, lo, hi, NMC=1_000_000, seed=20261024):
    rng = np.random.default_rng(seed)
    u = rng.uniform(np.log(lo), np.log(hi), NMC)
    t = np.exp(u)
    idx = np.searchsorted(lp, u)
    best = np.full(t.size, np.inf)
    for k in range(-3, 4):
        ii = np.clip(idx+k, 0, lp.size-1)
        cand = np.exp(lp[ii])
        dv = np.abs(cand/t - 1.0)*100.0
        best = np.minimum(best, dv)
    return float(np.mean(best < thr)*100.0)

def mc_family(thr, lo, hi, Nrat, NMC=200_000, seed=20261025):
    """Nrat 個獨立隨機比（同一域），至少一個命中。"""
    rng = np.random.default_rng(seed)
    hit = np.zeros(NMC, dtype=bool)
    for _ in range(Nrat):
        u = rng.uniform(np.log(lo), np.log(hi), NMC)
        t = np.exp(u)
        idx = np.searchsorted(lp, u)
        best = np.full(t.size, np.inf)
        for k in range(-3, 4):
            ii = np.clip(idx+k, 0, lp.size-1)
            cand = np.exp(lp[ii])
            best = np.minimum(best, np.abs(cand/t - 1.0)*100.0)
        hit |= (best < thr)
    return float(np.mean(hit)*100.0)

# ── sector 數據 ──
# PDG 2024 MS̄ @2GeV（同 01-08 遷入腳本用嘅值，保持一致）
quark_pdg = {"u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18, "t": 172.69}
m_e, m_mu, m_tau = 0.510998950e-3, 0.1056583755, 1.77686
hadron_ratios = {
    "π−π0 / n−p": 4.5936/1.29333205,
    "π−π0 / d−u": 4.5936/2.3,
    "n−p / d−u": 1.29333205/2.3,
    "K / π−π0": 3.937/4.5936,
    "K / n−p": 3.937/1.29333205,
    "K / d−u": 3.937/2.3,
}
hadron_abs = {  # Δm/m_e 形式（05 [C]）
    "(π−π0)/m_e": 4.5936/0.510998950,
    "(n−p)/m_e": 1.29333205/0.510998950,
    "(K)/m_e": 3.937/0.510998950,
    "(d−u)/m_e": 2.3/0.510998950,
}

sectors = {}
qnames = list(quark_pdg)
qratios = {}
for i, qi in enumerate(qnames):
    for j, qj in enumerate(qnames):
        if i < j:
            qratios[f"{qi}/{qj}"] = quark_pdg[qi]/quark_pdg[qj]
sectors["quark(PDG,15比)"] = qratios
sectors["lepton(3比)"] = {"μ/e": m_mu/m_e, "τ/μ": m_tau/m_mu, "τ/e": m_tau/m_e}
sectors["hadron差比(6比)"] = hadron_ratios
sectors["hadron絕對(4比)"] = hadron_abs

out = {"pool_size": int(pool.size), "sectors": {}}
print("\n" + "-"*78)
print("每命中: dev → 池身份 | 單比 MC P1 | 家族 MC PN（同域同 N）")
print("-"*78)
for sname, ratios in sectors.items():
    lo, hi = min(ratios.values()), max(ratios.values())
    rows = []
    print(f"\n[{sname}]  域 [{lo:.4f}, {hi:.4f}]  N={len(ratios)}")
    for rname, t in sorted(ratios.items(), key=lambda kv: kv[1]):
        d, p = dev_of(t)
        rows.append({"name": rname, "target": t, "dev_pct": d, "pool_value": p})
        if d < 1.0:  # 只對 <1% 命中做 MC（其餘冇判別意義）
            p1 = mc_single(d, lo, hi)
            pn = mc_family(d, lo, hi, len(ratios))
            rows[-1]["P1_pct"] = p1
            rows[-1]["PN_pct"] = pn
            flag = "🔴 隨機級" if pn > 5 else ("🟡" if pn > 1 else "🟢低密度")
            print(f"  {rname:<14} = {t:>11.6f}  dev {d:8.4f}%  池 {p:.6f} | P1 {p1:6.2f}%  PN {pn:6.2f}%  {flag}")
        else:
            print(f"  {rname:<14} = {t:>11.6f}  dev {d:8.4f}%  池 {p:.6f} | (dev≥1% 免 MC)")
    # 家族級 baseline：任何 dev 都用固定 0.1% / 1% 作參考
    base = {}
    for thr in (0.1, 0.05, 0.02):
        base[str(thr)] = {"P1": mc_single(thr, lo, hi), "PN": mc_family(thr, lo, hi, len(ratios))}
    out["sectors"][sname] = {"domain": [lo, hi], "N": len(ratios), "hits": rows, "baseline_by_thr": base}
    print("  參考密度（固定 thr）: " + " | ".join(
        f"{k}%: P1 {v['P1']:.2f}% PN {v['PN']:.2f}%" for k, v in base.items()))

print("\n" + "="*78)
print("判讀規則（事前）: PN>5% ⟹ 該命中隨機級（後驗，唔升格）; PN<1% ⟹ 池密度不足解釋")
print("              任何 sector 命中若要升格，必須 PN<1% 且樣本外盲測過")
print("="*78)

with open("data/09_density_gate.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("→ data/09_density_gate.json")
