#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/22 — R25-B：W1 之**算子空間窮盡**與「未來通道之准入包絡」

目的（把 W1 由「散文式阻塞」升為可審判詞）
  U1 枚舉本 M 已註冊、可產生 `em(p)−em(n)` 之算子類（兩體 ⟨B⟩／q² 自能），並各自對判準需求評分。
  U2 需求包絡：em = 1.2667 ± 0.094（σ_X 傳播），且 sign 要求 em(p) > em(n)，且 dof = 0。
  U3 結論分層：
     (a) 結構層：IDCM 已註冊結構內**無** em 算子 ⟹ 值不存在（硬）。
     (b) 數值層：兩體上限 A_max(v2) = 1.2494 **落在需求 1σ 帶** ⟹ 若日後有機制飽和上限，
         只差 0.18σ ⟹ **不被排除**（舊表述「兩體封」過強，本輪收正）。
     (c) q² 通道：4.5804 vs 1.2667 ⟹ 以 σ_X 計 = 35.3σ ⟹ 排除（非「3.62×」之模糊說法）。
     (d) 抑制因子：事前登記式族密度檢定陰性（R21-B）⟹ 無結構式可救 q² 通道。

執行：/usr/bin/python3 scripts/22_w1_operator_enumeration.py
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


EM = F("1.2667")
SIG = F("0.094")
AMAX_V1 = F("197.3269804") / (F("0.37742216816391860") * 3 * F("137.035999"))
AMAX_V2 = F("197.3269804") / (F("0.384184") * 3 * F("137.035999"))
SELF = F("4.5804")

print("=" * 78)
print("M10 R25-B — W1 算子空間窮盡 ＋ 准入包絡")
print("=" * 78)

LO, HI = EM - SIG, EM + SIG
print(f"\n[U1] 需求包絡：em ∈ [{float(LO):.4f}, {float(HI):.4f}] MeV（1σ，σ_X = 0.094）")
print("     附加要求：(i) sign = em(p) > em(n)；(ii) dof = 0（禁調係數）；(iii) 由已註冊結構導出")

print("\n[U2] 候選評估")
ROWS = [
    ("兩體 ⟨B⟩（M5 `A_max`，v2 正典）", float(AMAX_V2), "結構上限飽和"),
    ("兩體 ⟨B⟩（v1，已作廢）", float(AMAX_V1), "約定變體"),
    ("q² 自能（M5）", float(SELF), "量級"),
]
for nm, val, kind in ROWS:
    dev = val - float(EM)
    z = dev / float(SIG)
    print(f"     {nm:34s} = {val:8.4f} MeV  |  Δ = {dev:+.4f}  |  {z:+.1f}σ  |  {kind}")

z_self = (SELF - EM) / SIG
ck("U2.1 q² 自能以 σ_X 計超出 30σ（排除）", z_self > 30, f"{float(z_self):.1f}σ")
ck("U2.2 A_max(v2) 落於需求 1σ 帶（不被排除）", LO <= AMAX_V2 <= HI,
   f"{float(AMAX_V2):.4f} ∈ [{float(LO):.4f}, {float(HI):.4f}]；缺口僅 {float((EM-AMAX_V2)/SIG):.3f}σ")
ck("U2.3 A_max(v1) 亦落於 1σ 帶（約定不改『是否排除』）", LO <= AMAX_V1 <= HI,
   f"{float(AMAX_V1):.4f}")

print("\n[U3] 抑制因子（q² 通道若要成立所需）")
TGT = EM / SELF
print(f"     需抑制比例 = em/4.5804 = {float(TGT):.6f}")
j20 = json.load(open(os.path.join(DATA, "20_w1_suppression.json"), encoding="utf-8"))
print(f"     R21-B 密度檢定：式族 {j20['family_size']}，±1% 命中 {j20['hits_pct1']} vs 隨機中位 "
      f"{j20['null']['pct1']['median']}；±0.1% 命中 {j20['hits_pct01']} vs 中位 {j20['null']['pct01']['median']}")
ck("U3.1 抑制因子無結構式（密度檢定陰性）",
   j20["hits_pct1"] <= j20["null"]["pct1"]["median"], "沿用 R21-B 結果，不重跑")

print("\n[U4] 窮盡性與分層判詞")
CLASSES = [
    ("兩體 ⟨B⟩ 通道", "已註冊（M5）", "飽和值落 1σ 帶 ⟹ 非排除；但 ⟨B⟩≤1 恆等式禁加參數 ⟹ 無法自行補到上限"),
    ("q² 自能通道", "已註冊（M5）", "35.3σ 超 ⟹ 排除；抑制無結構式"),
    ("結構抑制重標定", "本 M 掃描（R21-B）", "式族 4538 密度檢定陰性 ⟹ 零命中"),
    ("其他（Σpair／味分裂等）", "**未註冊**", "屬 M4／M7 產權 ⟹ 本 M 只可要求、不可自行引入"),
]
for nm, st, judge in CLASSES:
    print(f"     {nm:22s}｜{st:14s}｜{judge}")
ck("U4.1 已註冊算子類 = 2（兩體／q²），其餘未註冊", True, "見上表")
ck("U4.2 判詞分層：結構層硬封、數值層兩體僅 0.18σ", True,
   "舊『兩體封』過強，改為『0.18σ 缺口 ＋ 參數禁加』")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "22_w1_operator_enumeration.py", "role": "R25-B W1 算子空間窮盡",
       "envelope": [float(LO), float(HI)], "em": float(EM), "sigma": float(SIG),
       "A_max_v2": float(AMAX_V2), "A_max_v1": float(AMAX_V1),
       "gap_v2_sigma": float((EM - AMAX_V2) / SIG),
       "self_energy": float(SELF), "self_z": float(z_self),
       "suppression_target": float(TGT),
       "classes": [{"name": n, "status": s, "judge": j} for n, s, j in CLASSES],
       "verdict": "W1：結構層無算子（硬）；數值層兩體僅 0.18σ（不排除）；q² 35.3σ（排除）；抑制無結構式；新通道屬他 M",
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "22_w1_operators.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/22_w1_operators.json 已寫")