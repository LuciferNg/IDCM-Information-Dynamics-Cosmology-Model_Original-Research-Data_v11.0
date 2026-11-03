#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
20_r2_composite_mods_push.py — R2 續推：{18,30,36} ↔ (u,d,s) 之最終裁決

R2 上一輪判 🟡（集合係 /45 加倍影像；γ 全序 canonical；接 u,d,s 屬命名步）。
本輪追問三件事：
  §1 {18,30,36} 之**內在**刻畫（不靠 /45 影像）：/90 三波共振（γ_a+γ_b=γ_c）之三個收束模
  §2 Galois 結構：γ36 = σ(γ18)、γ27 = σ(γ9)；γ30 = γ9+γ27 = 3；恆等式 γ36 = γ18/γ9
  §3 質量比檢定：六個 S₃ 排列全掃（禁以單一排列事後挑選）
  §4 裁決
輸出 data/20_r2_composite_mods_push.json
"""
from fractions import Fraction as Fr
import itertools
import json
import pathlib
import mpmath as mp

mp.mp.dps = 45
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OK, BAD = [], []


def chk(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f"  |  {detail}" if detail else ""))
    (OK if cond else BAD).append(label)
    return cond


g90 = lambda n: 4 * mp.sin(mp.pi * n / 90) ** 2

print("=" * 82)
print("§1  內在刻畫：/90 三波共振之「收束模」集合")
print("=" * 82)
res = []
for a in range(1, 44):
    for b in range(a + 1, 44):
        for c in range(b + 1, 45):          # c ≤ 44：排除 Nyquist γ₄₅ = 4（γ_a+γ_{45−a}=4 係恆等，非共振）
            if abs(g90(a) + g90(b) - g90(c)) < mp.mpf("1e-40"):
                res.append((a, b, c))
print(f"  γ_a+γ_b = γ_c（a<b<c ≤ 45）之解：{len(res)} 組")
for t in res:
    print(f"    {t}   γ{t[0]} + γ{t[1]} = γ{t[2]}   （{float(g90(t[0])+g90(t[1])):.10f}）")
sums = sorted({c for _, _, c in res})
print(f"  收束模（c）集合 = {sums}")
chk("★ {18,30,36} = 三波共振之收束模集合（內在刻畫，不需 /45 影像）", sums == [18, 30, 36])
chk("共振組數 = 3（每組對應一個收束模）", len(res) == 3, f"{res}")

print()
print("=" * 82)
print("§2  Galois 結構與恆等式")
print("=" * 82)
V = {n: g90(n) for n in (9, 18, 27, 30, 36)}
print(f"  γ₉  = {mp.nstr(V[9], 12)}   = (3−√5)/2")
print(f"  γ₁₈ = {mp.nstr(V[18], 12)}   = (5−√5)/2")
print(f"  γ₂₇ = {mp.nstr(V[27], 12)}   = (3+√5)/2")
print(f"  γ₃₀ = {mp.nstr(V[30], 12)}   = 3")
print(f"  γ₃₆ = {mp.nstr(V[36], 12)}   = (5+√5)/2")
chk("γ₃₆ = σ(γ₁₈)（√5 → −√5）", abs(V[36] - (5 + mp.sqrt(5)) / 2) < mp.mpf("1e-40")
    and abs(V[18] - (5 - mp.sqrt(5)) / 2) < mp.mpf("1e-40"))
chk("γ₂₇ = σ(γ₉)", abs(V[27] - (3 + mp.sqrt(5)) / 2) < mp.mpf("1e-40")
    and abs(V[9] - (3 - mp.sqrt(5)) / 2) < mp.mpf("1e-40"))
chk("γ₁₈ + γ₃₆ = 5（Galois 對之和為有理）", abs(V[18] + V[36] - 5) < mp.mpf("1e-40"))
chk("γ₉ + γ₂₇ = 3 = γ₃₀（共振）", abs(V[9] + V[27] - 3) < mp.mpf("1e-40"))
chk("★ 恆等式 γ₃₆ = γ₁₈ / γ₉（精確）", abs(V[36] - V[18] / V[9]) < mp.mpf("1e-40"),
    f"{mp.nstr(V[18]/V[9], 12)} vs {mp.nstr(V[36], 12)}")

print()
print("=" * 82)
print("§3  質量比檢定（六個 S₃ 排列全掃；PDG 2024 MS-bar 2 GeV）")
print("=" * 82)
# PDG：m_u = 2.16 (+0.49/−0.26)、m_d = 4.67 (+0.48/−0.17)、m_s = 93.4 (+8.6/−3.4) MeV
MU, MU_s = 2.16, 0.49
MD, MD_s = 4.67, 0.48
MS, MS_s = 93.4, 8.6
print(f"  m_u = {MU} ± {MU_s}｜m_d = {MD} ± {MD_s}｜m_s = {MS} ± {MS_s} MeV")
obs = {"d/u": MD / MU, "s/d": MS / MD, "s/u": MS / MU}
obs_s = {"d/u": obs["d/u"] * ((MU_s / MU) ** 2 + (MD_s / MD) ** 2) ** 0.5,
         "s/d": obs["s/d"] * ((MD_s / MD) ** 2 + (MS_s / MS) ** 2) ** 0.5,
         "s/u": obs["s/u"] * ((MU_s / MU) ** 2 + (MS_s / MS) ** 2) ** 0.5}
for k in obs:
    print(f"  {k} = {obs[k]:.4f} ± {obs_s[k]:.4f}")

cands = {"18": V[18], "30": V[30], "36": V[36]}
best = None
rows = []
for perm in itertools.permutations(cands):
    u, d, s = (cands[k] for k in perm)
    pred = {"d/u": d / u, "s/d": s / d, "s/u": s / u}
    z = {k: (float(pred[k]) - obs[k]) / obs_s[k] for k in obs}
    chi2 = sum(v * v for v in z.values())
    rows.append((chi2, perm, z))
    if best is None or chi2 < best[0]:
        best = (chi2, perm, z)
rows.sort()
for chi2, perm, z in rows[:3]:
    print(f"  (u,d,s) = {perm} ⟹ d/u {z['d/u']:+.2f}σ、s/d {z['s/d']:+.2f}σ、s/u {z['s/u']:+.2f}σ  χ² = {chi2:.1f}")
chk("★ 最佳排列仍被排除（χ² > 9 ⟹ >3σ 合併）", best[0] > 9, f"最佳 χ² = {best[0]:.1f}")
chk("★ s/d 預言 1.206 vs 實測 20.0（差 16.6×）—— 質量讀法致命", abs(rows[0][2]["s/d"]) > 3,
    f"{rows[0][2]['s/d']:+.2f}σ")
chk("d/u 雖合（|z| < 1）但 m_u 之不確定度達 23% ⟹ 無鑑別力", abs(rows[0][2]["d/u"]) < 1
    and MU_s / MU > 0.2, f"σ(m_u)/m_u = {MU_s/MU:.2f}")

print()
print("=" * 82)
print("§4  裁決")
print("=" * 82)
print("""  ① {{18,30,36}} 有**內在**刻畫：恰為 /90 三波共振之三個收束模（非僅 /45 影像）—— 集合層 ✅
  ② 排序 canonical（γ 嚴格遞增）；且 γ₃₆ = γ₁₈/γ₉（精確）
  ③ 但**質量讀法判死**：最佳 S₃ 排列 χ² = {:.1f}；s/d 預言 1.206 vs 實測 20.0（16.6×）
  ④ 唯一合者 d/u（0.02σ）之觀測不確定度 23% ⟹ 三個比之中一個「合」屬巧合可期
  ⟹ R2 最終判：🔴（質量／味標籤讀法）；僅餘「集合 + 全序」之純結構陳述（無物理內容）""".format(best[0]))
chk("R2 由 🟡（命名步）改判 🔴（質量讀法否證）", best[0] > 9)

print()
print("=" * 82)
print(f"§5  自審：PASS {len(OK)}、FAIL {len(BAD)}")
print("=" * 82)
for b in BAD:
    print(f"  ✗ {b}")

(DATA / "20_r2_composite_mods_push.json").write_text(json.dumps({
    "resonance_triples": [list(t) for t in res],
    "closing_modes": sums,
    "identity_g36_eq_g18_over_g9": str(mp.nstr(V[18] / V[9], 12)),
    "obs_ratios": {k: round(v, 4) for k, v in obs.items()},
    "obs_sigma": {k: round(v, 4) for k, v in obs_s.items()},
    "best_perm": list(best[1]), "best_chi2": round(best[0], 2),
    "best_z": {k: round(v, 2) for k, v in best[2].items()},
    "verdict": "R2 red: mass/flavour reading refuted (best S3 permutation chi2 = %.1f)" % best[0],
    "n_pass": len(OK), "n_fail": len(BAD),
}, indent=1, ensure_ascii=False), encoding="utf-8")
print("\n[data] 20_r2_composite_mods_push.json 已寫")
print("== done rc=%d" % (0 if not BAD else 1))
raise SystemExit(0 if not BAD else 1)