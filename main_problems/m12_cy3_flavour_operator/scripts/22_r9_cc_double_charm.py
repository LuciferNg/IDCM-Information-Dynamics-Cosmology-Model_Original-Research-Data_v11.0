#!/usr/bin/env python3
"""R9 — 雙粲 Ξcc 同位旋雙重態：二體對荷之 (c,c) 扇區檢定

資料源（唯讀 vendored）
  data/pdg2024/mass_width_2024.txt  （PDG 2024；**不含 Ξcc**，MC-ID 未定）
  data/22_lhcb2026_xicc_source.txt  （LHCb 2026 新測，逐字引文＋URL＋擷取日期）

Δm 定義：Δm(Ξcc⁺⁺ − Ξcc⁺) = m(ucc) − m(dcc) = **+1.77 MeV**
         （來源給 m(Ξcc⁺) − m(Ξcc⁺⁺) = −1.77 ± 0.84 ± 0.15 (+1.90/−1.30)）
         σ(stat⊕syst) = √(0.84² + 0.15²) = 0.8533；lifetime syst +1.90/−1.30 另列。

★ 計數慣例（明文，避免與 R4/R6 混淆）
Δm = Σ_(無序價對，計重數) a(f_i,f_j)。
  Ξcc⁺⁺ = {u,c,c} ⟹ 對 = (u,c)×2 + (c,c)×1
  Ξcc⁺  = {d,c,c} ⟹ 對 = (d,c)×2 + (c,c)×1
⟹ Δm = 2[a(u,c) − a(d,c)] = **−2δ_c**（δ_x ≡ a(d,x) − a(u,x)）
  （(c,c) 對相消 ⟹ 此通道**不**探測重-重相關項，只探測 δ_c 之扇區普適性）

預註冊判準（跑前凍結）
--------------------------------------------------
T1 二體對荷：預言 −2δ_c = +0.564 ± 0.568。|z| < 3 ⟹ (c,c) 扇區通；(c,c) 為**首個雙粲雙重態**。
T2 「Ξ-型偏移」假設：若 Ξc 家族常數 (+4.8262 ± 0.3254，docs/16 B1) 亦適用，
   預言 = 0.564 + 4.826 = 5.390。|z| < 3 ⟹ 該假設存活；否則判死。
T3 一體荷（R4-P1 已死之形式，作對照）：預言 −δ_u = +0.0887 ± 0.1124。
T4 內部一致性：由 Ξcc 反推 δ_c^eff = −Δm/2；與 Σc 錨之 δ_c 比較（z）。
T5 lifetime 系統誤差（+1.90/−1.30）不得忽略：報 T1/T2 之 z 區間。
禁：不得由本通道回擬 δ_c；不得因 lifetime syst 而挑選有利端。
"""
import pathlib
import hashlib
from decimal import Decimal
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40

# ---- PDG 錨（與 scripts/19,20 同源同值）
import re
SRC = pathlib.Path("data/pdg2024/mass_width_2024.txt")
MS = {}
for line in SRC.read_text().splitlines():
    if line.startswith("*"):
        continue
    m = re.match(r"^\s*(-?\d+)\s+(.+)$", line)
    if not m:
        continue
    nums = re.findall(r"[+-]?\d\.\d+E[+-]\d+", m.group(2))
    if len(nums) < 2:
        continue
    MS[m.group(1)] = (F(Decimal(nums[0])) * 1000, abs(F(Decimal(nums[1]))) * 1000)

Q = {  # 質量（MeV, Fraction）
    "n": MS["2112"][0], "p": MS["2212"][0], "Sm": MS["3112"][0], "S0": MS["3212"][0],
    "Xm": MS["3312"][0], "X0": MS["3322"][0], "Sc0": MS["4112"][0], "Scp": MS["4212"][0],
    "Sb-": MS["5112"][0], "Sb+": MS["5222"][0],
}
SG = {"n": MS["2112"][1], "p": MS["2212"][1], "Sm": MS["3112"][1], "S0": MS["3212"][1],
      "Xm": MS["3312"][1], "X0": MS["3322"][1], "Sc0": MS["4112"][1], "Scp": MS["4212"][1],
      "Sb-": MS["5112"][1], "Sb+": MS["5222"][1]}

d_s = (Q["Xm"] - Q["X0"]) / 2
d_d = (Q["Sm"] - Q["S0"]) - d_s
d_u = (Q["n"] - Q["p"]) - d_d
d_c = (Q["Sc0"] - Q["Scp"]) - d_d

# σ（線性組合傳播）
def sd(*terms):
    """terms: (label, coeff)"""
    return mp.sqrt(sum(mp.mpf(float(c * SG[k])) ** 2 for k, c in terms))


sd_ds = sd(("Xm", F(1, 2)), ("X0", F(-1, 2)))
sd_dd = sd(("Sm", 1), ("S0", -1), ("Xm", F(-1, 2)), ("X0", F(1, 2)))
sd_du = sd(("n", 1), ("p", -1), ("Sm", -1), ("S0", 1), ("Xm", F(1, 2)), ("X0", F(-1, 2)))
sd_dc = sd(("Sc0", 1), ("Scp", -1), ("Sm", -1), ("S0", 1), ("Xm", F(1, 2)), ("X0", F(-1, 2)))

# ---- LHCb 2026
DM_OBS = F(Decimal("1.77"))          # Δm(Ξcc++ − Ξcc+) = +1.77
SG_STAT = mp.mpf("0.84")
SG_SYST = mp.mpf("0.15")
SG_LIFE_UP = mp.mpf("1.90")
SG_LIFE_DN = mp.mpf("1.30")
SG_MEAS = mp.sqrt(SG_STAT ** 2 + SG_SYST ** 2)

print("=" * 94)
print("R9 — 雙粲 Ξcc 同位旋雙重態：二體對荷之 (c,c) 扇區")
print("=" * 94)
print("  PDG 2024 sha256(前32) =", hashlib.sha256(SRC.read_bytes()).hexdigest()[:32])
print("  LHCb 來源 sha256(前32) =", hashlib.sha256(pathlib.Path("data/22_lhcb2026_xicc_source.txt").read_bytes()).hexdigest()[:32])
print("\n[錨] δ 值（Fraction 精確；σ 線性傳播）")
for nm, v, s in [("δ_u", d_u, sd_du), ("δ_d", d_d, sd_dd), ("δ_s", d_s, sd_ds), ("δ_c", d_c, sd_dc)]:
    print("   %-4s = %+9.4f ± %7.4f MeV" % (nm, float(v), float(s)))
print("\n[觀測] Δm(Ξcc⁺⁺ − Ξcc⁺) = %+0.4f ± %0.4f MeV（stat⊕syst）；lifetime syst +%.2f/−%.2f"
      % (float(DM_OBS), float(SG_MEAS), float(SG_LIFE_UP), float(SG_LIFE_DN)))

# ---- T1
print("\n" + "=" * 94)
pred1 = -2 * d_c
sg1 = 2 * sd_dc
res1 = mp.mpf(float(DM_OBS)) - mp.mpf(float(pred1))
tot1 = mp.sqrt(sg1 ** 2 + SG_MEAS ** 2)
z1 = res1 / tot1
print("[T1] 二體對荷：預言 −2δ_c = %+8.4f ± %6.4f MeV   （2[a(u,c)−a(d,c)]；(c,c) 對相消）"
      % (float(pred1), float(sg1)))
print("     實測 %+8.4f ± %6.4f ；殘差 %+8.4f ；合成 σ %6.4f ；z = %+6.2f  ⟹ %s"
      % (float(DM_OBS), float(SG_MEAS), float(res1), float(tot1), z1,
         "通 🟢（(c,c) 扇區延續）" if abs(z1) < 3 else "判死 🔴"))
for nm, extra in [("lifetime +1.90", SG_LIFE_UP), ("lifetime −1.30", SG_LIFE_DN)]:
    t = mp.sqrt(tot1 ** 2 + extra ** 2)
    print("     含 %s 之合成 σ = %6.4f ⟹ |z| = %5.2f" % (nm, float(t), abs(res1 / t)))

# ---- T2
print("\n" + "=" * 94)
MU_XI = mp.mpf("4.8262")
SG_MU = mp.mpf("0.3254")
pred2 = mp.mpf(float(pred1)) + MU_XI
tot2 = mp.sqrt(sg1 ** 2 + SG_MU ** 2 + SG_MEAS ** 2)
z2 = (mp.mpf(float(DM_OBS)) - pred2) / tot2
print("[T2] 「Ξ-型偏移」假設（Ξc 家族常數 +4.8262 ± 0.3254，docs/16 B1）")
print("     預言 = %+8.4f ± %6.4f MeV ；實測 %+8.4f ⟹ z = %+6.2f  ⟹ %s"
      % (float(pred2), float(tot2), float(DM_OBS), z2,
         "存活 🟡" if abs(z2) < 3 else "判死 🔴（(c,c) 不帶 Ξ-型偏移）"))

# ---- T3
print("\n" + "=" * 94)
pred3 = -d_u
sg3 = sd_du
tot3 = mp.sqrt(sg3 ** 2 + SG_MEAS ** 2)
z3 = (mp.mpf(float(DM_OBS)) - mp.mpf(float(pred3))) / tot3
print("[T3] 一體荷（R4-P1 已死形式，對照）：預言 −δ_u = %+8.4f ± %6.4f ⟹ z = %+6.2f  %s"
      % (float(pred3), float(sg3), z3, "通 🟢" if abs(z3) < 3 else "判死 🔴"))

# ---- T4
print("\n" + "=" * 94)
dc_eff = -DM_OBS / 2
sg_eff = mp.mpf(float(SG_MEAS)) / 2
zc = (mp.mpf(float(dc_eff)) - mp.mpf(float(d_c))) / mp.sqrt(sg_eff ** 2 + sd_dc ** 2)
print("[T4] 內部一致性：由 Ξcc 反推 δ_c^eff = −Δm/2 = %+8.4f ± %6.4f" % (float(dc_eff), float(sg_eff)))
print("     與 Σc 錨 δ_c = %+8.4f ± %6.4f 比較 ⟹ z = %+6.2f  %s"
      % (float(d_c), float(sd_dc), zc, "一致 🟢" if abs(zc) < 3 else "不一致 🔴"))
print("     註：Δm/2 = %+0.4f MeV 與 |δ_c| = %0.4f MeV 同量級，此通道對 δ_c 之靈敏度甚高（σ_eff %0.3f）"
      % (float(DM_OBS / 2), float(abs(d_c)), float(sg_eff)))

print("\n" + "=" * 94)
print("[匯總 — 只報事實]")
print("   T1 二體對荷 (−2δ_c = %+0.4f)：z = %+0.2f ⟹ %s" % (float(pred1), z1, "通" if abs(z1) < 3 else "死"))
print("   T2 Ξ-型偏移 (+%0.4f)：z = %+0.2f ⟹ %s" % (float(MU_XI), z2, "存活" if abs(z2) < 3 else "判死"))
print("   T3 一體荷（對照）：z = %+0.2f ⟹ %s" % (z3, "通" if abs(z3) < 3 else "死"))
print("   T4 δ_c 扇區普適性：z = %+0.2f ⟹ %s" % (zc, "一致" if abs(zc) < 3 else "不一致"))
print("=" * 94)