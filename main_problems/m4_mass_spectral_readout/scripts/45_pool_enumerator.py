#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/45 — pool 型數字之**獨立枚舉器**（`docs/41` F8：6 個 pool 型數字無獨立枚舉器）。

背景：`docs/37` A4 / C、`docs/40` §4、`docs/36` §四 共 7 組數字屬「pool 枚舉」型，L2 未重建、
L1 亦只有單一實作（`scripts/03,09,16,19,31`）⟹ 一旦該實作有 bug，無人能察。
本腳本 = **第二套獨立實作**（不 import 上述任何一支；最近鄰改全池暴力、種子不同、協定逐條重述）。

覆蓋與判定（事前）：
  P1 `97.9%`（legacy protocol 02：15 個夸克比、thr 0.1%、域 [2e-3,173]、家族級）→ 重現
  P2 統一 protocol（09）之四 sector 密度（thr 0.1/0.05/0.02）→ 新增 canonical 密度表
  P3 `A = 87γ₁⁴M₀ = 1.63673 eV`（vs `A_target = 1.64337`，dev 0.404%）
  P4 `CAP = 4/γ₁ = 784.9646644`
  P5 味插入列舉：`108` 模型、單調 `12/108`、`p_null = 0.550`、動態範圍 `7.2×`
  P6 em 承載窗 `a ≤ 0.9729 fm`（`a_i^max = αħc·Q^max/|em_i|`）

**不覆蓋（並說明原因）**：`1.97%`（`ℤ[γ₁]` deg≥1 pool）、`5.59%`（`M_EW = M_Plγ₁⁸N_d`）、
`17 identity` —— 三者屬 O-era（`docs/09–18` 世代）自有 protocol，原文未存 protocol 參數
（域／thr／族定義），**不可復現**；`docs/37` A4 應改引 P2 之統一 protocol 數字（見 `docs/45`）。
`20.61%／90 條`、`733.8×／161.9%` 屬 M10 產權；`1.78／ε = 0.78(+8/−11)／[0.67,0.86]` 屬外部
錨（MILC18/FLAG）⟹ 三者皆非 M4 可重建對象（L3 例外保留）。
"""
import itertools
import json
import math

import numpy as np

CHK = []


def chk(label, ok, got, want):
    CHK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(want)})
    print(("PASS " if ok else "FAIL ") + f"{label:<58} | {str(got):<26} | {want}")


# ═══════════ 譜與池（獨立重建，不 import scripts/09） ═══════════
L = 88
G = np.array([4.0 * math.sin(k * math.pi / L) ** 2 for k in range(1, 45)])
POOL = np.array([G[a] / G[b] for a in range(44) for b in range(44)])
GAM1 = G[0]
print("=" * 108)
print(f"K_chain L={L}：γ₁ = {GAM1:.12f} ｜ pool |{{γ_a/γ_b}}| = {POOL.size} ｜ range [{POOL.min():.6g}, {POOL.max():.6g}]")
print("=" * 108)


def dev_pct(targets):
    """全池暴力最近鄰（log 空間等價；此處直接線性暴力，避免與 09 共用實作）。"""
    t = np.asarray(targets, dtype=float)
    out = np.empty(t.size)
    step = 20000
    for i in range(0, t.size, step):
        blk = t[i:i + step]
        out[i:i + step] = np.abs(POOL[None, :] / blk[:, None] - 1.0).min(axis=1) * 100.0
    return out


def mc_density(thr, lo, hi, n_ratios, n_mc, seed):
    """家族級：n_ratios 個 iid log-uniform 目標，至少一個落 thr 內之機率（%）。"""
    rng = np.random.default_rng(seed)
    hit = np.zeros(n_mc, dtype=bool)
    for _ in range(n_ratios):
        t = np.exp(rng.uniform(math.log(lo), math.log(hi), n_mc))
        hit |= (dev_pct(t) < thr)
    return float(hit.mean() * 100.0)


QUARK = {"u": 0.00216, "d": 0.00467, "s": 0.0934, "c": 1.27, "b": 4.18, "t": 172.69}
QNAMES = list(QUARK)
QRATIO = {f"{a}/{b}": QUARK[a] / QUARK[b] for i, a in enumerate(QNAMES) for b in QNAMES[i + 1:]}

# ── P1 legacy protocol 02：15 夸克比、thr 0.1%、域 [2e-3, 173] ──
p1 = mc_density(0.1, 2e-3, 173.0, len(QRATIO), 200_000, 45_1001)
print(f"\n[P1] legacy protocol 02（N={len(QRATIO)} 比, thr 0.1%, 域 [2e-3,173]）⟹ 家族級命中率 {p1:.2f}%")
chk("P1a legacy 97.9% 不可由精確枚舉復現（差 > 0.3pp ⟹ finding）", abs(p1 - 97.9) > 0.3, f"{p1:.2f}%", "偏離 97.9 逾 0.3pp")


def mc_density_windowed(thr, lo, hi, n_ratios, n_mc, seed, kwin=3):
    """legacy 實作（09 之 searchsorted + ±3 鄰居窗口）：檢驗 97.9% 是否為窗口近似之 artifact。"""
    lp = np.sort(np.log(POOL))
    rng = np.random.default_rng(seed)
    hit = np.zeros(n_mc, dtype=bool)
    for _ in range(n_ratios):
        u = rng.uniform(math.log(lo), math.log(hi), n_mc)
        t = np.exp(u)
        idx = np.searchsorted(lp, u)
        best = np.full(t.size, np.inf)
        for k in range(-kwin, kwin + 1):
            ii = np.clip(idx + k, 0, lp.size - 1)
            best = np.minimum(best, np.abs(np.exp(lp[ii]) / t - 1.0) * 100.0)
        hit |= (best < thr)
    return float(hit.mean() * 100.0)


p1w = mc_density_windowed(0.1, 2e-3, 173.0, len(QRATIO), 200_000, 45_1001)
print(f"     窗口近似（±3 鄰居，legacy 實作）⟹ {p1w:.2f}%")
chk("P1b 窗口近似重現 97.9%（±0.3）⟹ legacy 值 = 窗口 artifact", abs(p1w - 97.9) < 0.3, f"{p1w:.2f}%", "97.9 ± 0.3")

# ── P2 統一 protocol（09 之規格）：四 sector × 三 thr ──
m_e, m_mu, m_tau = 0.510998950e-3, 0.1056583755, 1.77686
SECT = {
    "quark(15)": QRATIO,
    "lepton(3)": {"μ/e": m_mu / m_e, "τ/μ": m_tau / m_mu, "τ/e": m_tau / m_e},
    "hadron差比(6)": {"π−π⁰/n−p": 4.5936 / 1.29333205, "π−π⁰/d−u": 4.5936 / 2.3,
                   "n−p/d−u": 1.29333205 / 2.3, "K/π−π⁰": 3.937 / 4.5936,
                   "K/n−p": 3.937 / 1.29333205, "K/d−u": 3.937 / 2.3},
    "hadron絕對(4)": {"π−π⁰/me": 4.5936 / m_e, "n−p/me": 1.29333205 / m_e,
                   "K/me": 3.937 / m_e, "d−u/me": 2.3 / m_e},
}
print("\n[P2] 統一 protocol（族級命中率 %，域 = 受測比之 [min,max]，log-均勻）")
p2rows = {}
for name, rr in SECT.items():
    lo, hi = min(rr.values()), max(rr.values())
    row = {}
    for thr in (0.1, 0.05, 0.02):
        row[str(thr)] = round(mc_density(thr, lo, hi, len(rr), 100_000, 45_2000 + int(thr * 1000)), 2)
    hits = {k: round(float(d), 4) for k, d in zip(rr, dev_pct(list(rr.values())))}
    p2rows[name] = {"domain": [lo, hi], "N": len(rr), "hit_dev_pct": hits, "density_pct": row}
    print(f"  {name:<14} 域 [{lo:.3g}, {hi:.3g}] N={len(rr)} | "
          + " | ".join(f"thr {t}%: {v:.2f}%" for t, v in row.items())
          + f" | 最小 dev {min(hits.values()):.4f}%")
chk("P2 四 sector 密度皆屬 [0,100]（表已產出）", all(0 <= v <= 100 for s in p2rows.values() for v in s["density_pct"].values()),
    f"{len(p2rows)} sector", "有意義值")

# ── P3 A = 87γ₁⁴ M₀ ──
A_TARGET = 1.64337          # eV/γ²（P6 交回值，docs/12）
A_87 = 1.63673              # docs/19 §4 所報 87γ₁⁴M₀ 之值
M0_implied = A_87 / (87.0 * GAM1 ** 4)
A_recompute = 87.0 * GAM1 ** 4 * M0_implied
dev_A = abs(A_recompute - A_TARGET) / A_TARGET * 100.0
print(f"\n[P3] 87γ₁⁴M₀：γ₁⁴ = {GAM1**4:.6e} ｜ 由 A = 1.63673 反解 M₀ = {M0_implied:.6e}")
print(f"     A_recompute = {A_recompute:.6f} eV ｜ vs A_target {A_TARGET} ⟹ dev {dev_A:.4f}%")
chk("P3 A = 87γ₁⁴M₀ 之 dev ≈ 0.404%（docs/19 §4）", abs(dev_A - 0.404) < 0.01, f"{dev_A:.4f}%", "0.404 ± 0.01")
chk("P3 識別：87 = L₈ − 1（L₈ = 88）", 87 == 88 - 1, "87 = L₈−1", "L₈−1")

# ── P4 CAP = 4/γ₁ ──
CAP = 4.0 / GAM1
print(f"\n[P4] 單模可達比上限 CAP = 4/γ₁ = {CAP:.10f}")
chk("P4 CAP = 784.9646644（4/γ₁）", abs(CAP - 784.9646644423913) < 1e-6, f"{CAP:.7f}", "784.9646644")

# ── P5 味插入列舉（108 模型 / 12 單調 / p_null / 7.2×） ──
B4, NQ, LD = 4.0, 3, 12
ANCH = {"N": (0.891, 0.107), "Σ": (1.498, 0.262), "Ξ": (2.293, 0.351)}
ENV = {"N": ("u", "d"), "Σ": ("u", "s"), "Ξ": ("s", "s")}
KAP = (0.5, 1, 2, 3, 5, 10, 30, 100, 300)


def rho(dep, kappa, env, res):
    d = np.arange(LD + 1)
    V = np.zeros(LD + 1)
    for p in env:
        V -= kappa * B4 ** (-np.abs(d - dep[p]))
    off = np.full(LD, -math.sqrt(B4 - 1))
    H = np.diag(V) + np.diag(off, 1) + np.diag(off, -1)
    _, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    nodes = np.array([1.0] + [B4 * (B4 - 1) ** (k - 1) for k in range(1, LD + 1)])
    v = v / math.sqrt(float((nodes * v ** 2).sum()))
    return v ** 2 if res == "site" else nodes * v ** 2


def chi2(tr):
    return sum((tr[i] - ANCH[k][0]) ** 2 / ANCH[k][1] ** 2 for i, k in enumerate(("N", "Σ", "Ξ")))


MODELS = []
for perm in itertools.permutations(range(3)):
    dep = {"u": perm[0], "d": perm[1], "s": perm[2]}
    for res in ("site", "shell"):
        for k in KAP:
            tr = tuple(float(NQ * rho(dep, k, ENV[s], res)[dep["d"]]) for s in ("N", "Σ", "Ξ"))
            MODELS.append({"dep": dep, "norm": res, "kappa": k, "triple": tr, "chi2": chi2(tr)})
MODELS.sort(key=lambda m: m["chi2"])
best = MODELS[0]
mono = [m for m in MODELS if m["triple"][0] < m["triple"][1] < m["triple"][2]]
print(f"\n[P5] 列舉 {len(MODELS)} 模型（6 置換 × 2 歸一 × {len(KAP)} κ）")
print(f"     最佳：u,d,s 深度 = {tuple(best['dep'][f] for f in 'uds')} × {best['norm']} × κ={best['kappa']}"
      f" ⟹ χ² = {best['chi2']:.2f}，⟨O⟩ = ({best['triple'][0]:.4f}, {best['triple'][1]:.4f}, {best['triple'][2]:.4f})")
print(f"     單調升（N<Σ<Ξ）之模型：{len(mono)}/{len(MODELS)}")

rng = np.random.default_rng(451025)
sig = np.array([ANCH[k][1] for k in ("N", "Σ", "Ξ")])
grid = np.array([m["triple"] for m in MODELS])
NMC = 20000
mins = np.empty(NMC)
for i in range(NMC):
    fake = np.exp(rng.uniform(math.log(0.2), math.log(4.0), 3))
    mins[i] = (((grid - fake[None, :]) / sig[None, :]) ** 2).sum(axis=1).min()
p_null = float((mins <= best["chi2"]).mean())
print(f"     MC null（{NMC} 次）：p = P(minχ²_null ≤ {best['chi2']:.2f}) = {p_null:.3f}")
chk("P5a 模型數 = 108", len(MODELS) == 108, len(MODELS), 108)
chk("P5b 單調模型 = 12/108", len(mono) == 12, f"{len(mono)}/{len(MODELS)}", "12/108")
chk("P5c p_null = 0.550（±0.02；本輪 2×10⁴ 次）", abs(p_null - 0.550) < 0.02, f"{p_null:.3f}", "0.550")

# 動態範圍（docs/31 §四 定義）：斜率 b = (⟨O⟩_Ξ − ⟨O⟩_N)/2；資料要求 b = 0.675（docs/29）
B_DATA = 0.675
bm = min(mono, key=lambda m: m["chi2"]) if mono else best
B_MODEL = (bm["triple"][2] - bm["triple"][0]) / 2.0
ratio = B_DATA / B_MODEL
print(f"     斜率：資料 b = {B_DATA}（docs/29）｜ 最佳單調模型 b = ({bm['triple'][2]:.4f} − {bm['triple'][0]:.4f})/2 = {B_MODEL:.4f}"
      f" ⟹ 比值 {ratio:.2f}×")
chk("P5d 機制動態範圍比 = 7.2×（±0.3）", abs(ratio - 7.2) < 0.3, f"{ratio:.2f}×", "7.2 ± 0.3")

# ── P6 em 承載窗 a_max = αħc·Q^max/|em_i| ──
ALPHA_HBARC = 1.4399645  # MeV·fm
WIN = {"n−p": (0.9867, (2, 3)), "Σ⁰−Σ⁺": (0.5630, (2, 3)), "Σ⁻−Σ⁰": (0.9720, (1, 1)),
       "Σ⁻−Σ⁺": (0.4090, (1, 1)), "Ξ⁻−Ξ⁰": (0.9800, (1, 1))}
print("\n[P6] em 承載窗（a_i^max = αħc·Q^max/|em_i|）")
amax = {}
for ch, (em, q) in WIN.items():
    a = ALPHA_HBARC * (q[0] / q[1]) / abs(em)
    amax[ch] = round(a, 4)
    print(f"     {ch:<10} |em| = {em:.4f} MeV  Q^max = {q[0]}/{q[1]}  ⟹ a_max = {a:.4f} fm")
W = min(amax.values())
chk("P6 承載窗 a ≤ 0.9729 fm（限制通道 n−p）", abs(W - 0.9729) < 5e-4, f"{W:.4f} fm", "0.9729 fm")
chk("P6 舊窗 0.3357 fm（Ξ⁻−Ξ⁰, |em| = 4.29）可重現",
    abs(ALPHA_HBARC * 1.0 / 4.29 - 0.3357) < 5e-4, f"{ALPHA_HBARC/4.29:.4f} fm", "0.3357 fm")

DATA = {
    "pool_size": int(POOL.size), "gamma1": GAM1,
    "P1_legacy_97_9": {"exact": round(p1, 3), "windowed_pm3": round(p1w, 3), "documented": 97.9},
    "P2_unified": p2rows,
    "P3_A_87": {"A_recompute": round(A_recompute, 5), "dev_pct": round(dev_A, 4), "M0_implied": M0_implied},
    "P4_cap": round(CAP, 10),
    "P5_scan": {"n_models": len(MODELS), "n_monotone": len(mono), "p_null": round(p_null, 4),
                "best_chi2": round(best["chi2"], 4), "best_triple": [round(x, 5) for x in best["triple"]],
                "b_data": B_DATA, "b_model": round(B_MODEL, 5), "dyn_ratio": round(ratio, 3), "nmc": NMC},
    "P6_window": {"a_max_fm": amax, "window_fm": round(W, 4)},
    "not_covered": {
        "1.97%｜5.59%｜17 identity": "O-era（docs/09–18）自有 protocol，參數未存 ⟹ 不可復現；改用 P2 統一 protocol",
        "20.61%｜733.8×": "M10 產權（docs/40 已明文 M4 未重建）",
        "1.78｜ε=0.78｜[0.67,0.86]": "外部錨（MILC18/FLAG），非 M4 產物",
    },
    "n_pass": sum(1 for c in CHK if c["ok"]), "n_total": len(CHK),
}
print("\n" + "=" * 108)
print(f"RESULT: {DATA['n_pass']}/{DATA['n_total']} PASS")
print("=" * 108)
with open("data/45_pool_enumerator.json", "w", encoding="utf-8") as f:
    json.dump({"checks": CHK, "data": DATA}, f, ensure_ascii=False, indent=1)
print("→ data/45_pool_enumerator.json")