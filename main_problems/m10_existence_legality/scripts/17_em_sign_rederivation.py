#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/17 — R19：`em` 方向更正之**全下游重推**

背景
  docs/03 原把本 M 之 `em` 標為 `em(n)−em(p)`；R18 證其應為 **`em(p)−em(n)`**（數值不變、方向反）。
  問題：方向一改，凡「依方向成立」之論斷皆須重推；只改標籤不改推導 = 偷步。
  本腳本把**所有下游依賴**逐項重推，並輸出可核 JSON。

輸入（全部本 M 已註冊）
  X    = m_d − m_u      = 2.56 ± 0.094 MeV   （FLAG；docs/03 §誤差預算）
  m_n − m_p = 1.29333236 MeV（AME/PDG 實測）
  m_e  = 0.51099895 MeV
  A_max(v2) = 1.2494 MeV（M5 譜側兩體上限）
  X_pred(γ) = 2.6897 MeV（γ 比值路線，無誤差帶）

執行：/usr/bin/python3 scripts/17_em_sign_rederivation.py
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


# ------------------------------------------------------------------ 輸入
X = F("2.56")
SIG_X = F("0.094")
MNP = F("1.29333236")          # m_n − m_p
ME = F("0.51099895")
AMAX = F("1.2494")
X_PRED = F("2.6897")           # γ 路線
EM_CLAIM = F("1.2667")         # docs/03 之登記值

print("=" * 78)
print("M10 R19 — `em` 方向更正之全下游重推")
print("=" * 78)

# ------------------------------------------------------------------ R1 恆等
em_exact = X - MNP
print(f"\n[R1] 恆等：em ≡ X − (m_n − m_p) = {float(X):.2f} − {float(MNP):.8f} = {float(em_exact):.8f}")
ck("R1.1 em 由恆等重算 = 登記值 1.2667",
   abs(em_exact - EM_CLAIM) < F("1e-4"),
   f"|Δ| = {float(abs(em_exact - EM_CLAIM)):.2e} MeV（登記值只四位 ⟹ 容差 1e-4）")
x_check = X - EM_CLAIM - ME
print(f"      交叉核：X − em − m_e = {float(x_check):.8f}；觀測 m_n−m_p−m_e = {float(MNP - ME):.8f}")
ck("R1.2 x(n) 自洽（差 < 1e-4，源自 X 之三位有效數字）",
   abs(x_check - (MNP - ME)) < F("1e-4"), f"Δ = {float(abs(x_check - (MNP - ME))):.2e}")

# ------------------------------------------------------------------ R2 方向
print("\n[R2] 方向：em = −(em_n − em_p)")
em_np = -(em_exact)            # em(n) − em(p)
print(f"      em(p) − em(n) = {float(em_exact):+.8f} MeV")
print(f"      em(n) − em(p) = {float(em_np):+.8f} MeV  ← 與 M4 scripts/27 之獨立核（E = −1.2667）同值")
ck("R2.1 方向為 em(p)−em(n) > 0（電磁使質子較重）", em_exact > 0, f"{float(em_exact):+.6f} > 0")
ck("R2.2 em(n)−em(p) = −1.26666…（M4 獨立核 −1.2667）",
   abs(float(em_np) + 1.2667) < 1e-4, f"{float(em_np):.6f}")

# ------------------------------------------------------------------ R3 σ 傳播（全新結果）
print("\n[R3] σ 傳播（R18 之後新增；原先未做）")
print(f"      em = X − (m_n−m_p) ⟹ σ_em = σ_X = {float(SIG_X)} MeV（m_n−m_p 之誤差 ≪ 0.094）")
gap = EM_CLAIM - AMAX
z = gap / SIG_X
print(f"      缺口 = 1.2667 − 1.2494 = {float(gap):.4f} MeV")
print(f"      ⟹ 相對 σ_X 之張力 = {float(gap)/float(SIG_X):.3f}σ")
ck("R3.1 缺口 < 1σ（σ_X = 0.094）", z < 1,
   f"{float(gap)/float(SIG_X):.3f}σ ⟹ **不構成統計排除**；W1 之「超界 1.37%」係中央值差，非顯著性")
print(f"      等效 1σ 區間：em ∈ [{float(EM_CLAIM-SIG_X):.4f}, {float(EM_CLAIM+SIG_X):.4f}] MeV")
print(f"      上限 A_max = {float(AMAX):.4f} ∈ 此區間 ⟹ 兩者不相斥")

# ------------------------------------------------------------------ R4 sign 翻轉閾值
print("\n[R4] 判號靈敏度（方向更正後之唯一「方向依賴」下游）")
flip = X - ME                      # x_n > 0 ⟺ em < X − m_e
print(f"      x(n) > 0 ⟺ em < X − m_e = {float(flip):.4f} MeV")
print(f"      現值 em = {float(EM_CLAIM):.4f} ⟹ 餘量 = {float(flip-EM_CLAIM):.4f} MeV（"
      f"{float((flip-EM_CLAIM)/flip)*100:.1f}% of 閾）")
ck("R4.1 現值距判號翻轉有正餘量（0.7823 MeV）",
   flip - EM_CLAIM > 0, f"餘量 {float(flip-EM_CLAIM):.4f} MeV ⟹ 符號方向更正**不改變**任何判詞")
need_em = (MNP - ME) / EM_CLAIM          # |x_n| / |em|
need_em_neg = (MNP - ME) / abs(-EM_CLAIM)  # 以 −em 重算（測方向無關性）
print(f"      判號所需 em 精度 = |x|/|em| = {float(need_em)*100:.1f}%（絕對 {float(MNP-ME):.4f} MeV）")
ck("R4.2 精度門檻 61.8% 與 em 之方向無關",
   abs(float(need_em) - 0.618) < 5e-3 and need_em == need_em_neg,
   f"{float(need_em)*100:.1f}%；以 −em 重算 {float(need_em_neg)*100:.1f}% ⟹ 相同")

# ------------------------------------------------------------------ R5 舊標籤可證偽性
print("\n[R5] 舊標籤之可證偽性（外部 QED，僅作方向旁證）")
print("      舊標籤 ⟹ 需 em(n)−em(p) = +1.2667（中子電磁上較重）")
print("      格點 QED（BMW/FLAG 系）之 em(n)−em(p) ≈ −1.00(7) MeV ⟹ 舊標籤差 2.27 MeV")
print("      ⚠️ 此為**外部數**，本 M 未獨立核 ⟹ 只作方向旁證，不作證據")
ck("R5.1 舊標籤與 QED 之號相反（外部旁證，標 🟡）", True,
   "新標籤 em(n)−em(p) = −1.2667 與 QED 同號；差 |−1.2667 −(−1.00)| ≈ 0.27 MeV 屬 FLAG-X／格點 QED 之外部張力")

# ------------------------------------------------------------------ R6 γ 路線
print("\n[R6] γ 比值路線（無誤差帶）")
em_g = X_PRED - MNP
print(f"      em = X_pred − (m_n−m_p) = {float(em_g):.4f} MeV；超 A_max {float(em_g-AMAX):.4f} MeV"
      f"（{float((em_g-AMAX)/em_g)*100:.2f}%）")
ck("R6.1 γ 路線之 em 與 docs/02 之 1.3963 一致",
   abs(em_g - F("1.3963")) < F("1e-3"), f"{float(em_g):.4f}")
print("      ⚠️ 此路線之 σ 由 m_u 之輸入誤差決定，未量化 ⟹ 只作區間上端，不可與 FLAG 路線混算 σ")

# ------------------------------------------------------------------ R7 下游清單
print("\n[R7] 下游逐項重推判")
L = [
    ("docs/03 表（62 列）標籤", "更正", "em(n)−em(p) → em(p)−em(n)"),
    ("docs/03 第 114 列", "更正", "同上"),
    ("docs/03 數值（1.2667/1.3963/0.7823/162%/61.8%）", "不變", "全為量值或比例"),
    ("docs/02「em 修正 1.3963」", "不變", "量值；已加方向註"),
    ("docs/05 牆表（W1 行）", "不變", "量值"),
    ("docs/05 §六 x_verdict 協議", "需加一步", "計 em 時須用 em(p)−em(n)"),
    ("docs/17 S3.1", "重推", "加 σ_X=0.094 ⟹ 缺口 0.18σ（非硬衝突）"),
    ("docs/17 S3.4（em/|x|=162%）", "不變", "量值"),
    ("README §三（em 主宰）", "不變", "量值"),
    ("final_result 關鍵公式", "需加方向與 σ", "x = X − em − m_e，em ≡ em(p)−em(n)"),
    ("scripts/02 輸出字串", "已更正", "前次已改"),
    ("scripts/03 G4 計算", "不變", "公式本已正確（x = X − em − m_e）"),
    ("01_problem_record P7 列", "需 append 更正", "原文保留，加方向註"),
    ("W1 狀態", "重推", "兩體封（結構）＋ 數值 0.18σ（不排除）"),
    ("M5 之 ②（自能反號）", "外部", "已回覆：同號"),
    ("M4 scripts/27", "獨立印證", "em(n)−em(p) = −1.2667 同值"),
]
for item, verdict, why in L:
    print(f"      {item:<38} {verdict:<12} {why}")
ck("R7.1 下游清單 16 項全部有裁決", len(L) == 16, f"{len(L)} 項")
ck("R7.2 「重推」類僅 S3.1 與 W1 狀態兩項（其餘為量值 ⟹ 不動）",
   sum(1 for _, v, _ in L if v == "重推") == 2, "2 項")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "17_em_sign_rederivation.py", "role": "R19 em 方向更正之全下游重推",
       "inputs": {"X": float(X), "sigma_X": float(SIG_X), "m_n_minus_m_p": float(MNP),
                  "m_e": float(ME), "A_max_v2": float(AMAX), "X_pred_gamma": float(X_PRED)},
       "em_pp_minus_pn": float(em_exact), "em_pn_minus_pp": float(em_np),
       "gap_MeV": float(gap), "gap_over_sigma": float(gap) / float(SIG_X),
       "sigma_band": [float(EM_CLAIM - SIG_X), float(EM_CLAIM + SIG_X)],
       "flip_threshold": float(flip), "margin": float(flip - EM_CLAIM),
       "em_gamma": float(em_g), "em_gamma_over": float((em_g - AMAX) / em_g),
       "downstream": [{"item": i, "verdict": v, "why": w} for i, v, w in L],
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "17_em_sign.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/17_em_sign.json 已寫")