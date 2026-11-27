#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""攻堅唯一出路：⑤ 之 Z（正規化）究竟是否未知？
   核實一：λ = μ²/(2⟨W⟩²) 是否精確 2/123（μ² = 2、⟨W⟩² = 61.5 = 90·41/60）
   核實二：g₄ = 24λ = 16/41（精確分數）
   核實三：μ² = 2 = γ_22(/88) = 4sin²(22π/88) 精確；但 M14 之 γ 係 /90（4sin²(πn/90)，Σγ=88）
           ⟹ 鏈錯配（Dirichlet /90 vs 環 /88）—— 必須標鏈
   核實四：16/41 與 λ_SM = 0.12907 之比；除以 N_gen 後之偏差（陷阱區判定）
"""
import json

import mpmath as mp

mp.mp.dps = 40
F = mp.mpf

# 一、λ 之結構
from fractions import Fraction as Fr
mu2 = F(2)
W2 = F(90) * F(41) / F(60)          # ⟨W⟩² = Ω_Λ·90 = 61.5
lam = mu2 / (2 * W2)
lam_ex = Fr(2, 1) / (2 * Fr(123, 2))     # 精確分數：μ²=2、⟨W⟩²=123/2
g4_ex = 24 * lam_ex
print("⟨W⟩² =", mp.nstr(W2, 20), "＝ 61.5 ?", W2 == F("61.5"), "｜精確 =", Fr(90*41, 60))
print("λ = μ²/(2⟨W⟩²) =", mp.nstr(lam, 25), "｜精確 =", lam_ex, "｜＝ 2/123 ?", lam_ex == Fr(2, 123))
print("1/λ =", mp.nstr(1 / lam, 20), "（應 = 61.5）")

# 二、g₄
g4 = 24 * lam
print("g₄ = 24λ =", mp.nstr(g4, 25), "｜精確 =", g4_ex, "｜＝ 16/41 ?", g4_ex == Fr(16, 41), "（= %.12f）" % float(g4_ex))

# 三、鏈核實
g22_88 = 4 * mp.sin(mp.pi * 22 / 88) ** 2
g22_90 = 4 * mp.sin(mp.pi * 22 / 90) ** 2
print("γ_22(/88) =", mp.nstr(g22_88, 20), "（μ²=2 精確）")
print("γ_22(/90) =", mp.nstr(g22_90, 20), "（≠2 ⟹ 鏈錯配風險）")
S88 = sum(4 * mp.sin(mp.pi * mp.mpf(k) / 88) ** 2 for k in range(1, 45))
S90 = sum(4 * mp.sin(mp.pi * mp.mpf(k) / 90) ** 2 for k in range(1, 45))
print("Σγ(/88, k=1..44) =", mp.nstr(S88, 20), "｜Σγ(/90, n=1..44) =", mp.nstr(S90, 20))

# 四、與 SM 對照（陷阱區判定）
lamSM = F("0.12907")
r = g4 / lamSM
print("g₄/λ_SM =", mp.nstr(r, 12), "｜g₄/(3λ_SM) =", mp.nstr(r / 3, 12),
      "｜g₄/3 =", mp.nstr(g4 / 3, 15), " 之偏差 =", mp.nstr(100 * (g4 / 3 / lamSM - 1), 6), "%")

# 三層 λ 互換比（O39）
lam_B = (mp.mpf("3.618033988749895") ) / 4      # (φ⁴−1)/4，φ⁴ = 3φ+1 = 3.854 → 修正用精確 φ
phi = (1 + mp.sqrt(5)) / 2
lam_B = (phi ** 4 - 1) / 4
print("λ_B = (φ⁴−1)/4 =", mp.nstr(lam_B, 20), "｜1/λ_B =", mp.nstr(1 / lam_B, 20))
ratio_layers = lam_B / lam
print("λ_B/λ_C =", mp.nstr(ratio_layers, 15), "｜24λ_B =", mp.nstr(24 * lam_B, 15), "｜24λ_C = g₄ =", mp.nstr(g4, 12))
print("註：μ²/λ = 2⟨W⟩² 為場重標度不變量 ⟹ 三層 λ 不可互換（§四十八）")

out = {"W2": mp.nstr(W2, 20), "lambda": mp.nstr(lam, 25), "lambda_exact": str(lam_ex), "lambda_is_2_over_123": lam_ex == Fr(2, 123),
       "g4": mp.nstr(g4, 25), "g4_exact": str(g4_ex), "g4_is_16_over_41": g4_ex == Fr(16, 41),
       "gamma22_88": mp.nstr(g22_88, 20), "gamma22_90": mp.nstr(g22_90, 20),
       "sum88_44": mp.nstr(S88, 20), "sum90_44": mp.nstr(S90, 20),
       "ratio_g4_over_lambdaSM": mp.nstr(r, 12), "dev_if_div3_pct": mp.nstr(100 * (g4 / 3 / lamSM - 1), 6),
       "lambda_B": mp.nstr(lam_B, 20), "ratio_layers_lambdaB_over_lambdaC": mp.nstr(lam_B / lam, 15), "g4_lambdaB_layer": mp.nstr(24 * lam_B, 15),
       "verdicts": {"Z_is_unknown": False,
                    "chain_mismatch": "M14 四腿 mod 90（Dirichlet /90，Σγ=88）vs W-field λ 建於 /88（Σγ=90）",
                    "sm_correspondence": "g₄/λ_SM ≈ 3.02；除 N_gen 才落 0.78% ⟹ 陷阱區（<3%）⟹ 不升格"}}
with open("/home/wsl/idcm_ckpt/e8/mrB/op_z_normalization.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
print("→ 落 mrB/op_z_normalization.json")
