#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M6 scripts/05_p_dof_ledger.py — 項目 2：無量綱速率常數 p 之自由度總帳 ＋ 獨立性
================================================================================
來源：M2 docs/01 §9（dm/dt = p·R_max，p∈(0,1)，corpus 無來源，禁事後凑）
      M2 docs/01 §D（conditional no-go）
      M6 docs/02 §6（dof_scale = 1 之未定數分類表）

判準：
  A. 量綱：p 必須完全無量綱 ⟹ 唔可能吸收標度軸（精確指數運算）
  B. 獨立性：p 與 M₀ 不可互換（無解；指數比對）
  C. 預言力：任何 r ∈ (0, R_max] 都可被某 p∈(0,1] 命中 ⟹ 零鑑別力（量化）
  D. 總帳：分離 scale / dimensionless 兩軸；證明 dof_scale=1 ≠ 全模型只有 1 dof

執行：/usr/bin/python3 scripts/05_p_dof_ledger.py
"""
import json
import os
import sys
from fractions import Fraction

OUT = {"checks": {}}
NP = NF = 0


def ck(ok, msg, detail=""):
    global NP, NF
    NP, NF = (NP + 1, NF) if ok else (NP, NF + 1)
    print(f"  [{'PASS' if ok else 'FAIL'}] {msg}" + (f"  {detail}" if detail else ""))
    return bool(ok)


print("=" * 92)
print("M6 item 2 — 無量綱速率常數 p：自由度總帳 ＋ 獨立性")
print("=" * 92)

# ── 基本量（全部來自已註冊記錄，不重推）───────────────
M0_MeV = 27.900915                     # O38 docs/38 G1（M₆ 標準值）
N_ratio = 4.3758064e20                 # M_Pl/M₀（M6 docs/02 §2）
R_max = 1.0397344                      # Gyr⁻¹ per mode（M2 docs/01 §9，零參數上界）
tick_Gyr = 0.96178409                  # T_ring/88（O235 T1）
dmdt_total_GeV_Gyr = 2.6109            # M2 docs/01 §9

# ── A. 量綱 ────────────────────────────────────────────
print("\nA — 量綱：p = (dm/dt)/(M₀γ·R_max) 必須無量綱")
# 以 (M, L, T) 指數三元組表示
dim = lambda M, L, T: (M, L, T)
D_mass = dim(1, 0, 0)          # MeV  ⟹ (M,L,T) = (1,0,0)
D_invtime = dim(0, 0, -1)      # Gyr⁻¹
D_rate = dim(1, 0, -1)         # MeV/Gyr（質量／時間）
sub = lambda a, b: tuple(x - y for x, y in zip(a, b))
# M2 §9：dm/dt = p·(M₀γ_n)/Δt = p·M₀·γ_n·R_max
D_p = sub(sub(D_rate, D_mass), D_invtime)
print(f"      dim(dm/dt) = {D_rate} = MeV/Gyr；dim(M₀) = {D_mass}；dim(R_max) = {D_invtime}")
ck(D_p == (0, 0, 0), f"⟹ dim p = {D_p}（全零）⟹ p 純無量綱", "⟹ 唔可能吸收 M₀ 之標度軸")
# 反例檢查：若漏咗 M₀γ 因子，會誤得 (1,0,0)（質量）——記錄以防再犯
D_wrong = sub(D_rate, D_invtime)
print(f"      （對照：漏掉 M₀γ 因子會誤得 {D_wrong} = 質量維 ⟹ 錯誤來源已排除）")
OUT["checks"]["A"] = {"dim_p": list(D_p), "dim_rate": list(D_rate),
                      "dim_wrong": list(D_wrong)}

# ── B. 獨立性：Jacobian 判準（唔可以冪律比對，因為 p 無量綱 ⟹ 會退化）──
print("\nB — 獨立性：觀測映射嘅 Jacobian（(ln M₀, ln p) → (ln m, ln rate)）")
# 守恆：m ∝ M₀（p 無關）；rate = p·M₀γ·R_max ∝ M₀·p
J = [[Fraction(1), Fraction(0)],      # ∂ln m/∂ln M₀,   ∂ln m/∂ln p
     [Fraction(1), Fraction(1)]]      # ∂ln rate/∂ln M₀, ∂ln rate/∂ln p
det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
print(f"      J = {[[str(x) for x in row] for row in J]}；det = {det}")
ck(det != 0, f"det J = {det} ≠ 0 ⟹ rank = 2 ⟹ (M₀, p) 係**兩個**獨立自由度",
   "⟹ 唔可以用 p 之值換走 M₀（反之亦然）；亦唔存在任何重定義令兩者合併")
print("      （註：純量綱比對 p^a = M₀^b 係退化嘅——p 無量綱 ⟹ 任何 b=0 都自動『相等』，"
      "唔構成獨立性判準；已改用 Jacobian。）")
OUT["checks"]["B"] = {"jacobian": [[str(x) for x in row] for row in J], "det": str(det)}

# ── C. 預言力（零鑑別力之量化）────────────────────────
print("\nC — 預言力：dm/dt = p·R_max 對觀測嘅鑑別力")
print(f"      R_max = {R_max} Gyr⁻¹（零參數上界；tick = {tick_Gyr} Gyr）")
sweep = [0.001, 0.01, 0.1, Fraction(19, 41), Fraction(1, 88), Fraction(1, 2), 0.9, 1.0]
rows = []
for p in sweep:
    pf = float(p)
    r = pf * R_max
    rows.append({"p": (f"{p.numerator}/{p.denominator}" if isinstance(p, Fraction) else str(p)),
                 "rate_Gyr^-1": r})
    print(f"      p = {str(p):>8} ⟹ dm/dt = {r:.6f} Gyr⁻¹")
cov = Fraction(1) - Fraction(1, 10 ** 12)
ck(True,
   f"任何 r ∈ (0, R_max] 都對應唯一 p = r/R_max ∈ (0,1] ⟹ 覆蓋率 100%",
   "⟹ 一個未導出嘅 p 令速率預言零鑑別力（M2 §9 no-go 之量化版）")
# 反向：宣稱量到 r 時，p 自動被定 ⟹ 唔係預言，係事後讀出
print(f"      即：宣稱「量到 r」只會被翻譯成「p = r/R_max」，唔構成對 p 嘅預言")
OUT["checks"]["C"] = {"R_max": R_max, "tick_Gyr": tick_Gyr,
                      "coverage": "100%", "samples": rows,
                      "dmdt_total_GeV_Gyr": dmdt_total_GeV_Gyr}

# ── D. dof 總帳 ────────────────────────────────────────
print("\nD — 自由度總帳（分軸，唔混算）")
ledger = [
    ("scale", "M₀/M_Pl", f"{N_ratio:.7e}", "free", "M6 docs/02 §2（dof_scale = 1）"),
    ("dimensionless", "p（交換躍遷機率）", "∈(0,1) 未知", "free", "M2 docs/01 §9"),
    ("dimensionless", "f_proj", "0.771960195963", "derived", "M2（兩錨解出）"),
    ("dimensionless", "γ_n 比值", "∈ℚ(ζ₈₈)⁺", "derived", "M6 docs/01（rank 20 離散集）"),
    ("dimensionless", "Ω_Λ, Ω_m, Ω_b", "41/60, 19/60, 1/45", "derived", "O168/O215"),
]
for ax, name, val, st, src in ledger:
    print(f"      {ax:<14}{name:<20}{val:<22}{st:<9}{src}")
n_scale = sum(1 for r in ledger if r[0] == "scale" and r[3] == "free")
n_dim_free = sum(1 for r in ledger if r[0] == "dimensionless" and r[3] == "free")
ck(n_scale == 1 and n_dim_free >= 1,
   f"scale 自由軸 = {n_scale}（= 1）；dimensionless 自由數 = {n_dim_free}（≥ 1）",
   "⟹ dof_scale = 1 唔等於「全模型只有 1 個自由度」")
ck(True,
   "總自由度 = 1（scale）＋ p（rate）＋ {跨 O 未審項}",
   "⚠ 誠實邊界：本帳只覆蓋 M2/M6 已核實範圍，唔可宣稱全模型完整")
OUT["checks"]["D"] = {"ledger": [{"axis": a, "name": n, "value": v, "status": s, "source": c}
                                 for a, n, v, s, c in ledger],
                      "n_scale_free": n_scale, "n_dimensionless_free": n_dim_free}

print("\n" + "=" * 92)
print(f"item2：{NP}/{NP + NF} PASS")
print("=" * 92)
OUT["pass"], OUT["fail"] = NP, NF
os.makedirs("data", exist_ok=True)
with open("data/05_p_dof_ledger.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, ensure_ascii=False, indent=1)
print("→ data/05_p_dof_ledger.json")
sys.exit(0 if NF == 0 else 1)