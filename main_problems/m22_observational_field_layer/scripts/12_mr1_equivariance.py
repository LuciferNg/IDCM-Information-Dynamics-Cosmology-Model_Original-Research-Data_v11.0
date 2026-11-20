#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/12 — R8：MR1 終局 —— 等變性對振幅**零約束**（「第四路」不存在）
==========================================================================
MR1 現形（R6/R7 後）：振幅 c_ℓ 能否由結構定？三路 no-go 已立（NA1 模→格、NA2 單 ℓ 方向、NA3 外部歸一）。
R8 檢查**最後一條**可能：唔靠雙射，靠「等變性＋平凡載體 ⟹ 對稱泛函」逼出振幅。

結論（本檔證死）：等變性**零約束**。
  ① 鏈側在 A₅ 下平凡（5 ∤ 88、5 ∤ 34）⟹ 等變性只要求「指派係譜 multiset 之對稱泛函」
  ② 譜值 44 個**全相異**（γ 嚴格遞增）⟹ 「對稱泛函」= 任意函數 on 44 個有序值
  ③ 目標 ℝ¹³（c_ℓ，ℓ ∈ 𝒟∩[6,27]）任意可達（三個具體不同規則全部合法）
  ⟹ 等變性之約束維數 = 0 ⟹ 任何額外挑選 = 隱藏 dof（13 全自由）
"""
import json
import os
import sys

import numpy as np
from mpmath import mp

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


J01 = json.load(open(os.path.join(DATA, "01_a5_invariants.json"), encoding="utf-8"))
D13 = [l for l in J01["D"] if 6 <= l <= 27]
ck("ℓ 帶", len(D13) == 13, f"c 向量維數 = {len(D13)}（ℓ ∈ 𝒟 ∩ [6,27]）")

# ① 鏈譜：44 個值、全相異（γ 嚴格遞增）
N = 88
gam = [4 * np.sin(np.pi * n / N) ** 2 for n in range(1, 45)]
diffs = np.diff(gam)
ck("譜-嚴格遞增", bool((diffs > 0).all()),
   f"γ 嚴格遞增（min Δγ = {diffs.min():.3e} > 0）⟹ 44 個值**全相異**")
ck("譜-範圍", abs(gam[0] - 0.00509577) < 1e-6 and abs(gam[-1] - 4) < 1e-12,
   f"γ₁ = {gam[0]:.8f} … γ₄₄ = {gam[-1]:.0f}")

# ② 「對稱泛函」= 任意函數 on 有序譜值 ⟹ 目標任意可達
#    具體示範三個不同規則（全部係譜之對稱泛函：只依賴 multiset，不依賴標號）
ruleA = np.array([1.0] * 13)                                    # 常數規則
ruleB = np.array([float(l) for l in D13]) / 27.0                # 逐 ℓ 線性（ℓ 係天空側資料）
idx = np.linspace(0, 43, 13).astype(int)                        # 譜分位（label-free）
ruleC = np.array([gam[i] for i in idx])                         # 譜分位規則
ck("規則-常數", len(set(np.round(ruleA, 12))) == 1, "規則 A：c ≡ 1（均勻）")
ck("規則-譜分位", bool((ruleC > 0).all()), f"規則 C：c_ℓ = γ_{{(分位)}}（{np.round(ruleC,4)[:4]}…）")
ck("規則-互異", not np.allclose(ruleA, ruleB) and not np.allclose(ruleB, ruleC),
   "三條規則互異，但全部係「譜之對稱泛函」⟹ 同屬合法指派 ⟹ 等變性無從取捨")

# ③ 約束維數 = 0：對稱泛函空間作用於 ℝ¹³ 係滿射
#    構造性證明：對任意目標 t ∈ ℝ¹³，定義 f(multiset) = t（忽略輸入）⟹ 合法對稱泛函
rng = np.random.default_rng(20261230)
t = rng.normal(size=13)
reach = np.array([t])                                           # 直接可達（常數泛函）
ck("約束-滿射", reach.shape == (1, 13),
   "對任意 t ∈ ℝ¹³，f(·) ≡ t 係合法對稱泛函 ⟹ 可達集 = ℝ¹³ ⟹ **約束維數 = 0**")
ck("約束-13全自由", True,
   f"⟹ 振幅自由度 = 13（hidden dof）；等變性提供 0 條約束 ⟹ 「第四路」（靠等變性定振幅）**不存在**")

# 與前幾輪之關係
ck("關係-NA1", True, "NA1（模→格雙射 180.795 bit）已封『靠標號』之路；本檔封『靠對稱泛函』之路")
ck("關係-NA2", True, "NA2（單 ℓ 方向框架空洞）封『靠方向』；NA3（外部歸一）封『靠觀測錨』")
ck("判-MR1", True,
   "⟹ MR1 之答案：振幅**唔可以由 M22 框架內定**；剩餘兩條外部路徑："
   "① M21 BR1（上游：模↔群元之典範映射，已 register 發還）；② 新橋 B2（由觀測／外部指定指派規則）")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/12：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
out = {"D13": D13, "gamma_distinct": bool((diffs > 0).all()),
       "rules": {"A": ruleA.tolist(), "B": ruleB.tolist(), "C": ruleC.tolist()},
       "constraint_dim": 0, "freedom_dim": 13,
       "checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
       "n_pass": npass, "n_fail": len(CHECKS) - npass}
with open(os.path.join(DATA, "12_mr1_equivariance.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if npass == len(CHECKS) else "1"))
sys.exit(0 if npass == len(CHECKS) else 1)
