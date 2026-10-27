#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/04 — R4 譜讀數唯一化（P4）

(A) L0 全譜 γ_a/γ_b（a≠b, 1..44）→ 最近命中、E 閘、MC null 10^7
(B) L0_sel = O209 選擇集 Sel(𝒮) = {1,2,3,4,5,7,9,11,13}（9 模式，引用不重推）
    → 72 值；並列 O209 唯一分配下嘅參考比值（τ↔3, μ↔5, e↔13）
(C) L1（× 結構整數 ±1）、L2（雙 γ 比）→ dof 階梯 ＋ look-elsewhere 位元預算
(D) holdout 目標（p/e、n/e）→ 樣本外檢定
(E) 判決

零 import 上游；譜 γ_n = 4sin²(πn/88) 只引用定義。
"""
import json
import os

import numpy as np
from mpmath import mp, mpf, log, sqrt

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

checks = []


def chk(name, cond, detail=""):
    checks.append({"check": name, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + name + ("  " + detail if detail else ""))


def gamma(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


print("=" * 74)
print("M9 scripts/04 — R4 譜讀數唯一化")
print("=" * 74)

# 目標：三個輕子比（CODATA 2022 direct entries）＋ 兩個 holdout
TARGETS = {
    "mu/e": (mpf("206.7682827"), mpf("0.0000046"), "lepton"),
    "tau/mu": (mpf("16.8170"), mpf("0.0011"), "lepton"),
    "tau/e": (mpf("3477.23"), mpf("0.23"), "lepton"),
    "p/e": (mpf("1836.152673426"), mpf("0.000000032"), "holdout"),
    "n/e": (mpf("1838.68366200"), mpf("0.00000074"), "holdout"),
}

G = {n: gamma(n) for n in range(1, 45)}
GFL = {n: float(G[n]) for n in G}


def ref_values(pairs):
    """pairs: list of (a,b) → γ_a/γ_b float"""
    return np.array([GFL[a] / GFL[b] for a, b in pairs], dtype=np.float64)


# ════════════════════════════════════════════════════════════════════
# (A) L0 全譜
# ════════════════════════════════════════════════════════════════════
print("\n[A] L0 = γ_a/γ_b（a≠b, 1..44）")
L0_pairs = [(a, b) for a in range(1, 45) for b in range(1, 45) if a != b]
L0 = ref_values(L0_pairs)
L0_sorted = np.sort(L0)
print(f"  族大小 = {len(L0)} 值；值域 [{L0.min():.4g}, {L0.max():.4g}]")

l0_rows = []
for k, (v, u, kind) in TARGETS.items():
    vf = float(v)
    i = int(np.argmin(np.abs(np.log(L0_sorted / vf))))
    cand = L0_sorted[i]
    dev_rel = cand / vf - 1
    sigma = abs(dev_rel) * vf / float(u)
    a, b = L0_pairs[int(np.argmin([abs(GFL[x] / GFL[y] - vf) for x, y in L0_pairs]))]
    l0_rows.append({"target": k, "kind": kind, "obs": mp.nstr(v, 11),
                    "sigma_obs": mp.nstr(u, 3), "best_pair": [a, b],
                    "best_value": f"{cand:.10g}", "dev_rel": f"{dev_rel:.6e}",
                    "dev_pct": f"{dev_rel * 100:.6f}", "sigma": f"{sigma:.4g}"})
    print(f"  {k:7s} obs={float(v):>12.6f}  L0 最近 = γ{a}/γ{b} = {cand:.8g}"
          f"  dev={dev_rel * 100:+.6f}%  ({sigma:.3g}σ, σ_obs={float(u):.2g})")

# E 閘 ＋ MC
ln_range = float(np.log(L0.max() / L0.min()))
rho = len(L0) / ln_range
print(f"  log 覆蓋 = {ln_range:.4f}；密度 ρ = {rho:.1f} / ln-unit")
NMC = 10_000_000
rng = np.random.default_rng(20261027)
t_mc = np.exp(rng.uniform(np.log(1.0), np.log(3477.0), NMC))
lnL0 = np.log(L0_sorted)
mc_l0 = []
for tau in (2.2e-8, 6.5e-5, 1e-5, 1e-4, 1e-3):
    idx = np.searchsorted(lnL0, np.log(t_mc))
    idx = np.clip(idx, 1, len(lnL0) - 1)
    d = np.minimum(np.abs(lnL0[idx] - np.log(t_mc)),
                   np.abs(lnL0[idx - 1] - np.log(t_mc)))
    p = float(np.mean(d <= np.log1p(tau)))
    E = min(1.0, 2 * tau * rho)
    mc_l0.append({"tau": tau, "p_mc": p, "E_local": E})
    print(f"  　τ={tau:<9g} P_MC = {p:.6f}   E_local = {E:.6f}")

# 5σ 所需容差（look-elsewhere 校正）
def tau_for_p(p_goal, n_family, span):
    """p_goal = 2·τ·N/span ⟹ τ"""
    return p_goal * span / (2 * n_family)


tau5_L0 = tau_for_p(2.87e-7, len(L0), ln_range)
print(f"  ★ L0 要 5σ 需 τ ≤ {tau5_L0:.3e}（相對）—— 而最接近嘅真值 dev = "
      f"{min(abs(float(r['dev_rel'])) for r in l0_rows[:3]):.3e}")

# ════════════════════════════════════════════════════════════════════
# (B) L0_sel：O209 選擇集（引用，不重推）
# ════════════════════════════════════════════════════════════════════
print("\n[B] L0_sel = Sel(𝒮) = {1,2,3,4,5,7,9,11,13}（O209 選擇律；本 M 只引用）")
SEL = [1, 2, 3, 4, 5, 7, 9, 11, 13]
L0s_pairs = [(a, b) for a in SEL for b in SEL if a != b]
L0s = ref_values(L0s_pairs)
print(f"  族大小 = {len(L0s)} 值；值域 [{L0s.min():.4g}, {L0s.max():.4g}]")

sel_rows = []
for k, (v, u, kind) in TARGETS.items():
    vf = float(v)
    j = int(np.argmin(np.abs(np.log(L0s / vf))))
    cand = L0s[j]
    dev_rel = cand / vf - 1
    sigma = abs(dev_rel) * vf / float(u)
    sel_rows.append({"target": k, "obs": mp.nstr(v, 11), "best_value": f"{cand:.8g}",
                     "dev_rel": f"{dev_rel:.6e}", "dev_pct": f"{dev_rel * 100:.6f}",
                     "sigma": f"{sigma:.4g}"})
    print(f"  {k:7s} L0_sel 最近 = {cand:.8g}  dev = {dev_rel * 100:+.4f}%  ({sigma:.3g}σ)")

# O209 唯一分配參考值（τ↔3, μ↔5, e↔13；只作 γ-ratio 讀法之參考，唔重推 |Δm|^α）
ref_pairs = {"mu/e ← γ5/γ13": (5, 13), "tau/mu ← γ3/γ5": (3, 5), "tau/e ← γ3/γ13": (3, 13)}
assign_rows = []
print("  O209 唯一分配下之 γ-ratio 讀值（純譜比，非 |Δm|^α 律）：")
for name, (a, b) in ref_pairs.items():
    r = float(G[a] / G[b])
    tgt = {"mu/e": float(TARGETS["mu/e"][0]), "tau/mu": float(TARGETS["tau/mu"][0]),
           "tau/e": float(TARGETS["tau/e"][0])}[name.split()[0]]
    dev = r / tgt - 1
    assign_rows.append({"reading": name, "value": f"{r:.8g}", "obs": tgt,
                        "dev_pct": f"{dev * 100:.4f}", "factor_off": f"{abs(tgt / r):.4g}"})
    print(f"    {name:18s} = {r:10.6g}   vs obs {tgt:10.6f}   dev {dev * 100:+9.2f}%"
          f"  (差 {abs(tgt / r):8.4g}×)")
chk("B1 O209 唯一分配下 γ-ratio 讀值全部錯 >50%",
    all(abs(float(r["dev_pct"])) > 50 for r in assign_rows),
    "; ".join(r["dev_pct"] for r in assign_rows))

# ════════════════════════════════════════════════════════════════════
# (C) dof 階梯 ＋ 位元預算
# ════════════════════════════════════════════════════════════════════
print("\n[C] dof 階梯（加自由度會唔會「解釋」到？）")
S_STRUCT = [3, 4, 5, 7, 11, 12, 13, 16, 22, 29, 44, 45, 47, 55, 88, 90]
L1 = np.concatenate([L0 * s for s in S_STRUCT] + [L0 / s for s in S_STRUCT])
L1 = np.sort(np.unique(np.round(np.log(L1), 14)))
L1 = np.exp(L1)
print(f"  L1 = L0 × s^±1，s ∈ {S_STRUCT}（{len(S_STRUCT)} 個）→ {len(L1)} 值")

L2 = np.unique(np.round(np.log(np.outer(L0, L0).ravel()), 14))
L2 = np.exp(L2)
print(f"  L2 = L0 × L0 → {len(L2)} 值")

ladder = []
for name, fam in [("L0", L0_sorted), ("L1", L1), ("L2", L2)]:
    lnf = np.log(fam)
    rng_ln = float(lnf[-1] - lnf[0])
    rho_f = len(fam) / rng_ln
    row = {"level": name, "n_values": int(len(fam)), "log_range": rng_ln,
           "density": rho_f}
    for k, (v, u, kind) in TARGETS.items():
        vf = float(v)
        i = int(np.argmin(np.abs(lnf - np.log(vf))))
        dev = fam[i] / vf - 1
        row[f"dev_pct_{k}"] = f"{dev * 100:.6f}"
        row[f"sigma_{k}"] = f"{abs(dev) * vf / float(u):.4g}"
    # MC：同一 10^7 log-uniform 目標
    for tau in (1e-5, 1e-4):
        idx = np.clip(np.searchsorted(lnf, np.log(t_mc)), 1, len(lnf) - 1)
        d = np.minimum(np.abs(lnf[idx] - np.log(t_mc)),
                       np.abs(lnf[idx - 1] - np.log(t_mc)))
        row[f"p_mc_{tau:g}"] = float(np.mean(d <= np.log1p(tau)))
    tau5 = tau_for_p(2.87e-7, len(fam), rng_ln)
    row["tau_needed_5sigma"] = f"{tau5:.3e}"
    ladder.append(row)
    print(f"  {name}: N={len(fam):>8d}  ρ={rho_f:>10.2f}/ln  5σ 需 τ≤{tau5:.3e}  "
          f"| μ/e 最好 dev={row['dev_pct_mu/e']:>11s}% ({row['sigma_mu/e']}σ)"
          f"  P_MC(1e-5)={row['p_mc_1e-05']:.4f}  P_MC(1e-4)={row['p_mc_0.0001']:.4f}")

# look-elsewhere 位元預算
print("\n  位元預算（5σ ≈ 21.5 bit）")
bits = []
for row in ladder:
    span_bits = float(np.log2(3477.0))
    fam_bits = float(np.log2(row["n_values"]))
    for tau in (1e-5, 1e-4):
        got = span_bits - fam_bits - float(np.log2(1 / (2 * tau)))
        bits.append({"level": row["level"], "tau": tau, "bits": got})
        print(f"    {row['level']} τ={tau:<7g} 淨資訊 = {got:+.1f} bit"
              f"（需 21.5）⟹ {'合格' if got >= 21.5 else '不足'}")
chk("C1 任何 dof 級別在 τ≥1e-5 之淨資訊 < 21.5 bit",
    all(b["bits"] < 21.5 for b in bits))

# ════════════════════════════════════════════════════════════════════
# (D) holdout（樣本外）
# ════════════════════════════════════════════════════════════════════
print("\n[D] holdout 目標（p/e, n/e）— 同一凍結語法")
hold = []
for k in ("p/e", "n/e"):
    v, u, kind = TARGETS[k]
    vf = float(v)
    out = {"target": k, "obs": mp.nstr(v, 12), "sigma_obs": mp.nstr(u, 3)}
    for name, fam in [("L0", L0_sorted), ("L0_sel", np.sort(L0s))]:
        lnf = np.log(fam)
        i = int(np.argmin(np.abs(lnf - np.log(vf))))
        dev = fam[i] / vf - 1
        out[name] = {"best": f"{fam[i]:.8g}", "dev_pct": f"{dev * 100:.6f}",
                     "sigma": f"{abs(dev) * vf / float(u):.4g}"}
        print(f"  {k:5s} {name:7s} 最近 = {fam[i]:10.6g}  dev = {dev * 100:+9.4f}%"
              f"  ({abs(dev) * vf / float(u):.3g}σ)")
    hold.append(out)

# ════════════════════════════════════════════════════════════════════
# (E) 判決
# ════════════════════════════════════════════════════════════════════
verdict = {
    "L0_zero_dof": "三個輕子比值全部落空（um 級 >1e6σ）⟹ σ 級譜比身份不存在",
    "L0_sel_O209": "唯一分配下 γ-ratio 讀值錯 1–2.7 個數量級 ⟹ 輕子質量唔可能係 γ 比值",
    "dof_ladder": "L1/L2 可命中但淨資訊 < 21.5 bit（look-elsewhere 未過）⟹ 零權重",
    "holdout": "p/e、n/e 同樣落空 ⟹ 唔係輕子特有，係譜比語言整體無法達 σ 級",
    "remaining_dof": "即使某比命中，O209 分配仍屬上游（本 M 不碰）⟹ 本 M 不主張任何勝出公式",
    "verdict": "🔴 輕子扇區無 σ 級結構比值（誠實負）",
}
print("\n[E] 判決")
for k, v in verdict.items():
    print(f"  {k:16s}: {v}")

out = os.path.join(DATA, "04_spectral.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump({"L0": l0_rows, "L0_MC": mc_l0, "L0_density": rho,
               "tau_needed_5sigma_L0": tau5_L0,
               "L0_sel": sel_rows, "assignment_reading": assign_rows,
               "dof_ladder": ladder, "bit_budget": bits, "holdout": hold,
               "verdict": verdict, "checks": checks,
               "pass": sum(1 for c in checks if c["pass"]), "total": len(checks)},
              f, ensure_ascii=False, indent=2)
print("\n" + "=" * 74)
print(f"PASS {sum(1 for c in checks if c['pass'])}/{len(checks)}   -> {out}")
print("=" * 74)