"""M5 docs/20 配套 — D4-A：contact 態依賴之**事前註冊封閉族**盲測（φ-lock 閘）
================================================================================================
背景（`docs/19` T4）：形狀比天花板 `r_K ≤ 13.4034`，數據 98.7778 ⟹ 差 7.37×。
等價講法：模型之 `|c| ∝ 1/(m_h m_l)` 需乘一個**態因子** `F(ℓ_h)`：

    |c|_K/|c|_D = (m_c/m_s)·(F_7/F_0)   ⟹   F_7/F_0 = ?
    數據 |c|_K/|c|_D = 0.0135/0.3525 = 9/235（Fraction 精確）
    m_c/m_s = 1.27/0.495 = 254/99
    ⟹ F_7/F_0 = (9/235)·(99/254) = 891/59690 = 1/66.9945…  ⟹ 目標 ≈ 1/67.0

本腳本**事前**註冊封閉候選族（6 項，全部由既有 IDCM 物件構成，不得事後增刪），逐項計 `F_7/F_0`
與目標之比，並以寫死判準判定。禁：看到 7.37 之後才加入族成員（φ-lock）；禁 φ-power 反解。

判準（寫死；先註冊後計算）：
  K1 必要性：候選須**先於**知曉目標值而可寫出（由 K_chain ⊗ K_tree 物件決定）。
  K2 命中：`|ln(F_cand/F_target)| ≤ ln(1.10)`（10% 內）。
  K3 升格另需**獨立預言**：同一公式須同時給出一個未用量（例如另一介子族之 r）之數值，
     該量不得參與本輪校準；否則記 `CANDIDATE-ONLY`（不得升格）。
  K4 拒絕：>10% 差 ⟹ `REFUTED`；若落在 ≤10% 但無獨立預言 ⟹ `HELD（疑凑數）`。
  K5 全族不命中 ⟹ 建議 D4-B（收檔 register），不得加族重試。

輸出：`data/22_registered_family.json`
"""
import json
from fractions import Fraction as F

# ---- 數據（精確 Fraction）----
C_K, C_D = F(135, 10000), F(3525, 10000)      # 0.0135 / 0.3525
M_C, M_S = F(127, 100), F(495, 1000)          # 1.27 / 0.495
TARGET = (C_K / C_D) / (M_C / M_S)            # F_7/F_0 之目標
TARGET_INV = 1 / TARGET

# ---- 既有 IDCM 物件（精確）----
LUCAS = {0: 2, 1: 1, 2: 3, 3: 4, 4: 7, 5: 11, 6: 18, 7: 29, 8: 47}
FIB = {0: 1, 1: 1, 2: 2, 3: 3, 4: 5, 5: 8, 6: 13, 7: 21, 8: 34}
B, ELL_D, ELL_K = 4, 0, 7

# ---- 事前註冊族（6 項；每項附結構來源）----
FAMILY = [
    ("F1", "樹節點重合機率：layer ℓ 有 (B−1)^ℓ 個節點 ⟹ 兩夸克同點機率 ∝ 1/(B−1)^ℓ",
     lambda l: F(1, (B - 1) ** l)),
    ("F2", "鏈-樹張量重合：K_chain 88 站 × 樹節點 ⟹ 每多一層多分母 88",
     lambda l: F(1, 88 ** l)),
    ("F3", "Lucas 譜：γ-族以 L_ℓ 標定層權重（ℓ ∈ K_φ 之 Lucas 分支）",
     lambda l: F(1, LUCAS[l]) * LUCAS[ELL_D]),
    ("F4", "Fibonacci 譜：F_ℓ 分支（Lucas 之伴侶族）",
     lambda l: F(1, FIB[l]) * FIB[ELL_D]),
    ("F5", "累積節點數：Σ_{k≤ℓ}(B−1)^k = ((B−1)^{ℓ+1}−1)/(B−2)",
     lambda l: F(B - 2, (B - 1) ** (l + 1) - 1) * ((B - 1) ** 1 - 1) / (B - 2)),
    ("F6", "鏈 Green 函數正規化：contact 支撐 ∝ 1/g(ℓ) 之鏈和（g 由 docs/01 精確式）",
     lambda l: F(12 * 88, (88 * 88 - 1) - 6 * l * (88 - l))),
]

print("=" * 104)
print("M5 docs/20（D4-A）— contact 態依賴：事前註冊封閉族盲測")
print("=" * 104)
print(f"\n[D0] 目標（由數據反解，非擬合）：")
print(f"     |c|_K/|c|_D = {C_K}/{C_D} = {C_K/C_D} = {float(C_K/C_D):.6f}")
print(f"     m_c/m_s = {M_C}/{M_S} = {M_C/M_S} = {float(M_C/M_S):.6f}")
print(f"     ⟹ F_7/F_0 目標 = {TARGET} = {float(TARGET):.6f} = 1/{float(TARGET_INV):.4f}")
print(f"     （即：K 之 contact 需較 1/(m_h m_l) 天然標度再縮小 {float(TARGET_INV):.2f}×）")

print(f"\n[K1–K4] 族成員逐項檢定（判準：10% 內 = 命中；另需獨立預言才升格）")
print(f"     {'#':4s}{'F_7/F_0':>16s}{'數值':>12s}{'目標比':>11s}{'ln 比':>9s}{'判':>18s}  結構來源")
rows = []
for tag, why, fn in FAMILY:
    val = fn(ELL_K) / fn(ELL_D)
    ratio = val / TARGET
    ok = abs(float(ratio).__class__(ratio) - 0) == 0  # placeholder, real check below
    lnr = __import__("math").log(float(ratio))
    if abs(lnr) <= __import__("math").log(1.10):
        verdict = "HELD（疑凑數）"
    else:
        verdict = "REFUTED"
    rows.append({"tag": tag, "value": str(val), "numeric": float(val), "ratio": float(ratio),
                 "ln_ratio": lnr, "verdict": verdict, "why": why})
    print(f"     {tag:4s}{str(val):>16s}{float(val):12.6f}{float(ratio):11.4f}{lnr:+9.3f}{verdict:>18s}  {why}")

print(f"\n[K3] 獨立預言檢查：族內無任何成員可提供「未用量」之預言（全部只含 ℓ 與常數）⟹ 即使 10% 內亦只可記 HELD")
print(f"\n[K5] 族內命中數（10% 內）= {sum(1 for r in rows if r['verdict'].startswith('HELD'))}／{len(rows)}")

# 誘餌檢查：1/L_4 = 1/7（5.4% 內）—— 必須明示拒絕理由
decoy = F(1, LUCAS[4])
print(f"\n[誘餌] 1/L_4 = 1/7 = {float(decoy):.6f} vs 目標 {float(TARGET):.6f} ⟹ 差 {abs(float(decoy/TARGET)-1)*100:.1f}%")
print(f"     拒絕理由：ℓ = 4 在 M5 結構中無任何角色（ℓ(d) = ℓ(s) = 7、ℓ(u) = ℓ(c) = 0）；")
print(f"     且無獨立預言 ⟹ 屬『看到 7.37 之後才有嘅形式』⟹ φ-lock 閘拒絕（`47 = L₈` 誘餌先例）")
print(f"     （`scripts/33` 之閘：φ^x 族須 x 事前指定於 K_φ，不得反解）")

out = {"target": str(TARGET), "target_inv": float(TARGET_INV), "rows": rows,
       "n_held": sum(1 for r in rows if r["verdict"].startswith("HELD")),
       "decoy_1_over_L4": {"value": float(decoy), "dev_pct": abs(float(decoy / TARGET) - 1) * 100,
                           "rejected": True,
                           "reason": "ℓ=4 無結構角色；無獨立預言；事後形式"},
       "verdict": "🔴 全族不命中或僅疑凑數 ⟹ 建議 D4-B（收檔 register）"}
with open("data/22_registered_family.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/22_registered_family.json")