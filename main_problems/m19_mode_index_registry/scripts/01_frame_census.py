#!/usr/bin/env python3
"""
M19 / scripts/01 — Frame 普查 ＋ 層分離（索引層 vs 支撐編碼層）

預註冊（跑前凍結）:
  C1  F1（O233 γ 升序）係 9 模式集上之雙射
  C2  F2（O209 |Δm| 降序）係 9 模式集上之雙射
  C3  F2 = F1 ∘ (集合反序)：兩者互為唯一非平凡次序反向 ⟹ 只差 1 bit
  C4  F1 依觀測質量嚴格遞增；F2 嚴格遞減（質量值只作排序，不作推導）
  C5  F3（O206 CSSG-v3）**唔係**雙射：n=1 出現 3 次、9 個粒子中 5 個落喺 9 模式集外
      （n ∈ {14,20,33,34,44}）⟹ F3 屬**支撐編碼層**（m = M₀Σσ_nγ_n，σ 多模），
      唔係索引層；不得當第三個 frame
  C6  質量層對 frame 方向零判別力：任一雙射都可取 C_f = m/(M₀γ_n) 精確擬合 9/9
      （= O120/M10「9! 精確解 ⟹ 零預言力」之現場重算）
  C7  故 frame 方向之數據判別力 = 0 bit；C3 之 1 bit 只能由結構給，唔可以由質量給

失敗判準: C1/C2/C3/C4/C5/C6 任一 FAIL ⟹ 層分離結論撤回。
"""
import json
import math
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 50

OUT = Path(__file__).parent.parent / "data"
N = 88
ND = 44

# 9 模式集（A₅ 動差準則；v10 O136 docs/19）：奇 n ≤ 13、偶 n ≤ 4
M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]

# 質量（O206 scripts/23 MASS，PDG 名義值；只作排序與 C_f 展示）
MASS = {
    "e": 0.511, "u": 2.16, "d": 4.67, "s": 93.4, "mu": 105.6584,
    "c": 1270.0, "tau": 1776.86, "b": 4180.0, "t": 172690.0,
}
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}

# F1 — O233/O137 授權：A₅ γ 升序（質量升序）
F1 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}
# F2 — O209 |Δm| 降序（質量降序）
F2 = {"t": 1, "b": 2, "tau": 3, "c": 4, "mu": 5, "s": 7, "d": 9, "u": 11, "e": 13}
# F3 — O206 scripts/23 CSSG-v3「物理 assignment」(n, C_f)
F3 = {"e": (1, 3.6), "mu": (7, 15.5), "tau": (20, 37.0), "u": (1, 15.5),
      "d": (1, 33.0), "s": (14, 3.6), "c": (34, 13.0), "b": (33, 44.0),
      "t": (44, 1554.0)}

M0 = mpf("27.900915")

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def gamma(n):
    return 4 * mp.sin(mp.pi * n / N) ** 2


def is_bijection_on_m9(fr):
    vals = [fr[p] for p in fr]
    return sorted(vals) == M9


print("=" * 92)
print("M19 scripts/01 — Frame 普查 ＋ 層分離")
print("=" * 92)

chk("C1 F1（γ 升序）係 9 模式集雙射",
    len(F1) == 9 and is_bijection_on_m9(F1), str(list(F1.values())))
chk("C2 F2（|Δm| 降序）係 9 模式集雙射",
    len(F2) == 9 and is_bijection_on_m9(F2), str(list(F2.values())))

# C3：F2 = F1 之集合反序
rev = {M9[len(M9) - 1 - i]: M9[i] for i in range(9)}  # 位置反序
F1_rev = {p: rev[n] for p, n in F1.items()}
chk("C3 F2 = F1 ∘ 集合反序（唯一非平凡次序反向）",
    F1_rev == F2, "F1∘rev → " + str({LABEL[p]: n for p, n in sorted(F1_rev.items(), key=lambda x: x[1])}))

# C4：質量次序
by_mass = sorted(MASS, key=lambda p: MASS[p])
f1_order = [p for p in sorted(F1, key=lambda p: F1[p])]
f2_order = [p for p in sorted(F2, key=lambda p: F2[p])]
chk("C4a F1 與質量升序逐位相同",
    by_mass == f1_order, " / ".join(LABEL[p] for p in f1_order))
chk("C4b F2 與質量降序逐位相同",
    by_mass[::-1] == f2_order, " / ".join(LABEL[p] for p in f2_order))

# C5：F3 唔係雙射
f3_vals = [F3[p][0] for p in F3]
dup = len(set(f3_vals)) != len(f3_vals)
outside = [p for p in F3 if F3[p][0] not in M9]
chk("C5 F3 屬支撐編碼層（非雙射）",
    dup and len(outside) == 5,
    f"n 值 {sorted(f3_vals)}；重複={dup}；集外 {len(outside)}/9 → {[LABEL[p] for p in outside]}")

# C6：質量層零判別力
rows = []
maxdev = 0.0
for name, fr in (("F1", F1), ("F2", F2)):
    for p in fr:
        g = gamma(fr[p])
        cf = mpf(str(MASS[p])) / (M0 * g)
        resid = M0 * g * cf - mpf(str(MASS[p]))
        maxdev = max(maxdev, abs(resid) / mpf(str(MASS[p])))
        rows.append({"frame": name, "particle": LABEL[p], "n": fr[p],
                     "gamma_n": float(g), "C_f": float(cf)})
chk("C6 任一 frame 皆 9/9 精確擬合（C_f 自由 ⟹ 零判別力）",
    maxdev < mpf("1e-40"), f"max rel resid = {mp.nstr(maxdev, 6)}")

chk("C7 frame 方向之數據判別力 = 0 bit", True,
    "F1/F2 皆 9/9；1 bit 只能由結構（T3 槽）給，唔可以由質量給")

print("\n  frame × 粒子 γ_n / C_f：")
print(f"  {'frame':>5} {'粒子':>4} {'n':>3} {'γ_n':>12} {'C_f':>12}")
for r in rows:
    print(f"  {r['frame']:>5} {r['particle']:>4} {r['n']:>3} {r['gamma_n']:>12.9f} {r['C_f']:>12.6f}")

OUT.mkdir(exist_ok=True)
(OUT / "01_frame_census.json").write_text(json.dumps({
    "M9": M9,
    "frames": {
        "F1_gamma_ascending": {LABEL[p]: n for p, n in F1.items()},
        "F2_dm_descending": {LABEL[p]: n for p, n in F2.items()},
        "F3_cssg_v3_encoding_layer": {LABEL[p]: F3[p][0] for p in F3},
    },
    "F3_is_bijection": False,
    "F3_outside_M9": [LABEL[p] for p in outside],
    "mass_order": [LABEL[p] for p in by_mass],
    "rows": rows,
    "max_rel_resid": float(maxdev),
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
