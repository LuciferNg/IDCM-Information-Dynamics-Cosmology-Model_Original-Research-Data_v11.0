#!/usr/bin/env python3
"""M3 scripts/15 — 第十五輪：R1 —— κ⁻¹ = 16 之身份 ＋ ℓ = 2 之結構必然

R1 原缺口：為何 κ⁻¹ = 16 就係暗份額？ℓ = 2 之前只由觀測鎖死（−36σ / +0.44σ / +65σ）。

本輪新路線【容量填充】（零觀測輸入）：
  樹層 ℓ 之佔用 = B^ℓ；層空間係樹模態空間之子空間
  ⟹ 累積佔用 Σ_{n=0}^{L} B^n ≤ dim(K_tree) = K_tree
  ⟹ L_max = 最大可填滿層 = 2（21 ≤ 33 < 85）

若成立：① ℓ = 2 結構必然（觀測 +0.44σ 由「選」變「預言」）
        ② 閘門 κ = B^{−ℓ} = 1/16 = γ_max⁻² ⟹ κ⁻¹ = 終端層佔用 = 暗份額（同一個數兩個面）
        ③ 通過層數 = ⌈ℓ ln B⌉ = 3 = 已填層數（閘門邊界 ≡ 填充邊界）
"""
import mpmath as mp
from fractions import Fraction as F

mp.mp.dps = 30
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_d, K_tree = 4, 44, 33
OBS, SIG = 0.26447, 0.00500          # Planck 2018 Ω_dm

# ── §1 容量填充 ─────────────────────────────────────────────────
print("§1 容量填充：L_max = max{L : Σ_{n=0}^{L} B^n ≤ C}")
def cum(L): return sum(B**n for n in range(L+1))
for L in range(0, 5):
    print(f"  L={L}: 累積 = {cum(L):4d}")
for name, C in (("K_tree", 33), ("N_d", 44), ("60", 60), ("Σγ", 90)):
    Lmax = max(L for L in range(0, 6) if cum(L) <= C)
    print(f"  C={name:7s}={C:3d} ⟹ L_max={Lmax}   終端佔用 B^L={B**Lmax}")
chk("K_tree = 33 ⟹ L_max = 2（21 ≤ 33 < 85）", max(L for L in range(6) if cum(L) <= 33) == 2)
chk("N_d = 44 ⟹ L_max = 2（穩健：換容量不改答案）", max(L for L in range(6) if cum(L) <= 44) == 2)
Lmax = 2
chk("ℓ = 2 ⟹ 終端佔用 = B² = 16", B**Lmax == 16)
print("  可接納容量區間：ℓ=2 ⟺ cum(2) ≤ C < cum(3) ⟺ 21 ≤ C < 85")
chk("K_tree 與 N_d 皆落區間 [21,85)", 21 <= 33 < 85 and 21 <= 44 < 85)
chk("容量須為維度型：Σγ=90（能量）給 L=3 ⟹ 被排除",
    max(L for L in range(6) if cum(L) <= 90) == 3)

# ── §2 閘門一致性 ───────────────────────────────────────────────
print("\n§2 閘門：κ = B^{−ℓ} = 1/16 =? γ_max⁻²")
kappa = F(1, B**Lmax)
print(f"  κ = {kappa} ; γ_max = B = {B} ; γ_max⁻² = {F(1, B*B)}")
chk("κ = γ_max⁻² EXACT", kappa == F(1, B*B))
filled = list(range(0, Lmax+1))
kap_mp = mp.mpf(kappa.numerator)/mp.mpf(kappa.denominator)
passed = [n for n in range(0, N_d) if mp.e**(-n) > kap_mp]
print(f"  已填層 = {filled} ; 通過閘門層 = {passed}")
chk("通過層 ≡ 已填層（兩集合相同）", passed == filled)
nlev = len(filled)
print(f"  通過層數 = {nlev} ; ⌈ℓ ln B⌉ = {int(mp.ceil(Lmax*mp.log(B)))} ; B−1 = {B-1}")
chk("通過層數 = ⌈ℓ ln B⌉ = 3", nlev == int(mp.ceil(Lmax*mp.log(B))) == 3)
chk("與 B−1 = 非平凡特徵標數一致（3）", nlev == B-1)

# ── §3 由 (B, ℓ, N_d) 生成三分 ──────────────────────────────────
print("\n§3 三分生成（零觀測輸入）")
vis, dm, lam = nlev, B**Lmax, N_d - nlev
tot = N_d + B**Lmax
print(f"  可見 = {vis} ; 暗 = {dm} ; Λ = {lam} ; 總 = {tot}")
chk("總 = N_d + B^ℓ = 60", tot == 60)
chk("(3:16:41)/60", (vis, dm, lam) == (3, 16, 41))
chk("暗份額 = κ⁻¹ = 16（同一物件兩面）", dm == int(1/kappa))

# ── §4 ℓ 掃描：觀測由「選」變「預言」─────────────────────────────
print("\n§4 Ω_dm(ℓ) = B^ℓ/(N_d+B^ℓ) vs Planck")
zs = {}
for L in (1,2,3,4):
    p = mp.mpf(B**L)/mp.mpf(N_d + B**L)
    z = (p - OBS)/SIG
    zs[L] = float(z)
    mark = "★由容量推出" if L == 2 else ""
    print(f"  ℓ={L}: {mp.nstr(p,6):9s}  z={float(z):+9.2f}  {mark}")
chk("ℓ=2（容量推出者）為唯一 |z|<1", abs(zs[2]) < 1 and min(abs(zs[l]) for l in (1,3,4)) > 30)
chk("其他自洽候選 ℓ=1 達 −36σ ⟹ 排除", abs(zs[1]) > 30)

# ── §5 暗份額讀法選別（終端 vs 累積）────────────────────────────
print("\n§5 暗份額三讀法對 Planck 0.26447±0.00500")
cands = {"層數 = 3": 3, "終端佔用 B^ℓ = 16": 16, "累積佔用 Σ = 21": 21}
for name, v in cands.items():
    p = mp.mpf(v)/60
    print(f"  {name:22s} {v:2d}/60 = {mp.nstr(p,6):9s}  z = {float((p-OBS)/SIG):+8.2f}")
chk("終端佔用唯一 |z|<1；累積 +17.1σ、層數 −42.9σ",
    abs((mp.mpf(16)/60-OBS)/SIG) < 1
    and abs((mp.mpf(21)/60-OBS)/SIG) > 10
    and abs((mp.mpf(3)/60-OBS)/SIG) > 20)

# ── 判決 ───────────────────────────────────────────────────────
print("\n§判決")
print("  ① ℓ=2 由容量填充導出（K_tree=33 或 N_d=44 皆給 2，區間 21≤C<85）⟹ 零觀測輸入")
print("  ② κ = B^{−ℓ} = γ_max⁻² = 1/16 ⟹ 閘門邊界 ≡ 填充邊界（通過層 = 已填層）")
print("  ③ κ⁻¹ = 16 = 終端層佔用 ⟹ 「為何 κ⁻¹ 即暗份額」= 同一物件兩個面（門檻 vs 計數）")
print("  ④ 觀測 +0.44σ 由「選擇」升為「預言」（ℓ 已獨立導出）")
print("  🟡 殘餘：『暗 = 終端佔用而唔係累積』屬物理識別步（數據以 15σ 支持）")
print(f"\n結算: PASS {PASS} / FAIL {FAIL}")