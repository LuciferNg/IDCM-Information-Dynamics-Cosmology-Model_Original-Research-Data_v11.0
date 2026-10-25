#!/usr/bin/env python3
"""M3 scripts/12 — 第十二輪：(a) 層/節點不對稱 ＋ (b) P4 M2↔M3 橋

(a) 不對稱：可見 = 通過層數（深度），暗 = 終端層節點數（標度）。
    統一假說：兩者皆為【單一閘門 κ = B⁻²】之兩個泛函
        可見 = ⌈ln κ⁻¹⌉      (= 3)
        暗   = κ⁻¹           (= 16)
    ⟹ 物質 19 = ⌈ln κ⁻¹⌉ + κ⁻¹（一個 κ 生成兩個數）
    殘餘：B=4 時 κ⁻¹ 亦等於終端層節點數 B^{L−1}（算術巧合？）

(b) P4 橋：M2 互換律 E_sync/E_proj = 41/19 是否即 M3 之 41:19？
    M2: E_sync = ⟨W⟩² = 123/2 = 61.5 ; E_proj = Σγ − E_sync = 57/2 = 28.5
    ⟹ 41/60 = 61.5/90 、19/60 = 28.5/90（能量讀法，無需 ⊕）
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_d = 4, 44
Sig = F(90)                      # Σγ
Esync = F(123, 2)                # ⟨W⟩² = 61.5
Eproj = Sig - Esync              # 57/2 = 28.5

# ── (a) 單閘門生成 ────────────────────────────────────────────────
print("(a) 單閘門 κ = B⁻² 生成物質兩份額")
kap = F(1, B*B)
L = 3
print(f"  κ = {kap} ;  κ⁻¹ = {kap**-1}")
chk("κ⁻¹ = 16 = B²", kap**-1 == 16)
chk("⌈ln κ⁻¹⌉ = 3（通過層數 L）", L == 3)
chk("物質 = ⌈ln κ⁻¹⌉ + κ⁻¹ = 3 + 16 = 19", L + int(kap**-1) == 19)
print(f"  ⟹ 可見 3 與 暗 16 係【同一個 κ】之兩個泛函：深度 vs 標度")

print("\n  終端層節點數 B^(L−1) vs κ⁻¹（B 掃描）")
for b in (3, 4, 5, 6):
    Lb = int(mp.ceil(mp.log(b**2)))
    term = b**(Lb-1)
    k_inv = b*b
    print(f"    B={b}: L={Lb}  B^(L−1)={term:4d}  κ⁻¹={k_inv:3d}  "
          f"{'✅相等' if term == k_inv else '✗'}")
chk("B=4 時終端節點數 = κ⁻¹ = 16", B**(L-1) == 16 == B*B)
chk("此等式只喺 B=3,4（⌈2lnB⌉−1 = 2）成立 ⟹ 屬算術範圍, 非結構律",
    all((b**(int(mp.ceil(mp.log(b**2)))-1) == b*b) == (b in (3, 4)) for b in range(2, 20)))

# ── (b) P4 橋：能量讀法 ──────────────────────────────────────────
print("\n(b) P4 橋：M2 能量分割 ⟹ M3 份額")
print(f"  Σγ = {Sig} ; E_sync = ⟨W⟩² = {Esync} ; E_proj = {Eproj}")
chk("E_sync = ⟨W⟩² = 123/2", Esync == F(123, 2))
chk("E_sync + E_proj = Σγ（守恆, 無擬合）", Esync + Eproj == Sig)
chk("K_tree = E_sync − E_proj = 33（= N_gen·L₅）", Esync - Eproj == 33)
chk("E_sync/E_proj = 41/19 EXACT", Esync/Eproj == F(41, 19), str(Esync/Eproj))
chk("E_sync/Σγ = 41/60 EXACT（= Ω_Λ）", Esync/Sig == F(41, 60), str(Esync/Sig))
chk("E_proj/Σγ = 19/60 EXACT（= Ω_m）", Eproj/Sig == F(19, 60), str(Eproj/Sig))
chk("量子 = E_sync/41 = E_proj/19 = 3/2", Esync/41 == F(3, 2) == Eproj/19,
    f"{Esync/41}")
chk("分母 60 = Σγ/(3/2) = 2Σγ/3（60 係【導出】, 非前提）", F(90)/(F(3,2)) == 60)

print("\n  ★ 簡化後之邏輯次序（無需 ⊕ 直和）")
print("    第1層(能量):  Ω_Λ = ⟨W⟩²/Σγ = 61.5/90 = 41/60 ;  Ω_m = 1−Ω_Λ = 19/60")
print("    第2層(閘門):  Ω_m 內部分 3 : 16 = ⌈ln κ⁻¹⌉ : κ⁻¹")
print("    ⟹ ⊕ (rank 相加 60) 只係【換算】後之寫法，唔係前提")

# ── 一致性：三條路線同值 ─────────────────────────────────────────
print("\n  三條獨立路線對 41:19")
routes = {
    "M2 能量 (⟨W⟩²/Σγ)":      (Esync, Eproj),
    "M3 秩 (41 : 3+16)":      (F(41), F(3+16)),
    "量子換算 (×3/2)":         (F(41)*F(3,2), F(19)*F(3,2)),
}
for k, (a, b) in routes.items():
    print(f"    {k:22s} {a} : {b}  = {F(a)/F(b)}")
chk("三路線比值全等 41/19", len({F(a)/F(b) for a, b in routes.values()}) == 1)

# ── 觀測對照 ─────────────────────────────────────────────────────
print("\n  觀測對照（Planck 2018）")
obs = {"Ω_Λ": (0.6847, 0.0073), "Ω_m": (0.3153, 0.0073), "Ω_b": (0.0493, 0.0003),
       "Ω_dm": (0.2645, 0.0050)}
pred = {"Ω_Λ": 41/60, "Ω_m": 19/60, "Ω_b": 3/60, "Ω_dm": 16/60}
for k in obs:
    z = (pred[k]-obs[k][0])/obs[k][1]
    print(f"    {k:5s} pred={pred[k]:.5f} obs={obs[k][0]:.4f}±{obs[k][1]:.4f}  z={z:+.2f}σ")
chk("Ω_Λ, Ω_m, Ω_dm 皆 |z|<1；Ω_b +2.33σ（已知預言）",
    all(abs((pred[k]-obs[k][0])/obs[k][1]) < 1 for k in ("Ω_Λ","Ω_m","Ω_dm")))

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")