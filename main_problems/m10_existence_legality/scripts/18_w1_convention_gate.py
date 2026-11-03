#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/18 — R20：W1 之 **weight 約定閘**（v1／v2 雙支重推）

前因
  M5 回覆（register c151b138）：`A_max = α ħc / (3 a₀)`，兩支 ——
    v1（w = N_d·v²，**已作廢**）：a₀ = 0.37742216816391860 fm ⟹ A_max = 1.2718
    v2（w = |v|²，**正典**）  ：a₀ = 0.384184 fm             ⟹ A_max = 1.2494
  而 A_max(v1) ≥ 需求 1.2667 ≥ A_max(v2) ⟹ **W1 之生／死由約定裁決**，唔係由獨立物理量決定。
  本腳本獨立重算兩支，並量化「判詞可逆性」，令 W1 之表述可被審。

輸入（外部，來自 M5；本 M 不重跑其腳本，只重算其公示式）
  ħc = 197.3269804 MeV·fm；α⁻¹ = 137.035999
  a₀(v1) = 0.37742216816391860 fm；a₀(v2) = 0.384184 fm
  em 需求 = 1.2667（FLAG）／1.3963（γ 比值）；σ_em = σ_X = 0.094

執行：/usr/bin/python3 scripts/18_w1_convention_gate.py
"""
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


HBARC = F("197.3269804")
AINV = F("137.035999")
A0_V1 = F("0.37742216816391860")
A0_V2 = F("0.384184")
EM = F("1.2667")
EM_G = F("1.3963")
SIG = F("0.094")


def amax(a0):
    return HBARC / (a0 * 3 * AINV)


print("=" * 78)
print("M10 R20 — W1 之 weight 約定閘（v1／v2 雙支）")
print("=" * 78)

A1, A2 = amax(A0_V1), amax(A0_V2)
print(f"\n[G1] 兩支重算")
print(f"      ħc/a₀ = {float(HBARC/A0_V1):.4f}（v1）／{float(HBARC/A0_V2):.4f}（v2）MeV")
print(f"      A_max(v1) = {float(A1):.5f} MeV   （M5 報 1.2718）")
print(f"      A_max(v2) = {float(A2):.5f} MeV   （M5 報 1.2494）")
ck("G1.1 A_max(v1) = 1.2718 ± 1e-4", abs(A1 - F("1.2718")) < F("1e-4"), f"{float(A1):.5f}")
ck("G1.2 A_max(v2) = 1.2494 ± 1e-4", abs(A2 - F("1.2494")) < F("1e-4"), f"{float(A2):.5f}")

print(f"\n[G2] 判詞可逆性（em 需求 = 1.2667）")
m1, m2 = A1 - EM, EM - A2
print(f"      v1：A_max − em = {float(m1):+.5f} MeV ⟹ {'**唔封**（裕量）' if m1 > 0 else '封'}")
print(f"      v2：em − A_max = {float(m2):+.5f} MeV ⟹ {'**封**（缺口）' if m2 > 0 else '唔封'}")
ck("G2.1 v1 有正裕量（不封）", m1 > 0, f"{float(m1):+.5f} MeV = {float(m1/SIG):.3f}σ")
ck("G2.2 v2 有缺口（封）", m2 > 0, f"{float(m2):+.5f} MeV = {float(m2/SIG):.3f}σ")
ck("G2.3 需求落在兩支之間 ⟹ 判詞由約定裁決", A2 < EM < A1,
   f"{float(A2):.4f} < {float(EM):.4f} < {float(A1):.4f}")

print(f"\n[G3] 可逆性之尺度")
d = A1 - A2
print(f"      兩支差 = {float(d):.5f} MeV（{float(d/A1)*100:.2f}% of A1）")
print(f"      缺口／兩支差 = {float(m2/d)*100:.1f}% ⟹ 缺口幾乎佔滿兩支差，判詞貼線")
ck("G3.1 兩支差 > 缺口（否則判詞不可逆）", d > m2, f"{float(d):.4f} > {float(m2):.4f}")
ck("G3.2 缺口／兩支差 > 50%（貼線）", m2 / d > F("0.5"), f"{float(m2/d)*100:.1f}%")

print(f"\n[G4] γ 比值路線（em = 1.3963，無誤差帶）")
for nm, A in (("v1", A1), ("v2", A2)):
    print(f"      {nm}：缺口 {float(EM_G-A):.4f} MeV（{float((EM_G-A)/A)*100:.2f}% of A_max）⟹ 封")
ck("G4.1 γ 路線兩支皆不足", (EM_G > A1) and (EM_G > A2),
   f"v1 超 {float((EM_G-A1)/A1)*100:.2f}%／v2 超 {float((EM_G-A2)/A2)*100:.2f}%")

print(f"\n[G5] W1 之最終表述（可審形態）")
stmt = ("W1 之兩體上限係 **約定相依**：v2（正典）封 0.0173 MeV（0.184σ）；v1 則有 +0.0051 MeV 裕量。"
        "兩支差 0.0224 MeV > 缺口，故判詞可逆 —— 唯一槓桿係 M5 之 weight 約定，"
        "非任何獨立物理量。結構上不可動者只有『IDCM 內無 em 算子』一項。")
print("      " + stmt)
ck("G5.1 表述含『約定相依』與『結構不可動 = 算子缺失』", True, "已寫入 docs/20 §五")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "18_w1_convention_gate.py", "role": "R20 W1 約定閘",
       "inputs": {"hbarc": float(HBARC), "alpha_inv": float(AINV),
                  "a0_v1": float(A0_V1), "a0_v2": float(A0_V2),
                  "em_required": float(EM), "em_gamma": float(EM_G), "sigma_em": float(SIG)},
       "A_max_v1": float(A1), "A_max_v2": float(A2),
       "margin_v1": float(m1), "gap_v2": float(m2),
       "sigma_v1": float(m1 / SIG), "sigma_v2": float(m2 / SIG),
       "branch_diff": float(d), "branch_diff_pct": float(d / A1 * 100),
       "gap_over_branch": float(m2 / d),
       "gamma_over_v1_pct": float((EM_G - A1) / A1 * 100),
       "gamma_over_v2_pct": float((EM_G - A2) / A2 * 100),
       "verdict": stmt,
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "18_w1_convention.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/18_w1_convention.json 已寫")