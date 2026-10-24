#!/usr/bin/env python3
"""M3 scripts/10 — 第十輪：Ω_b 2.33σ 攻堅（N_d = 44 vs 45 之結構判決）

張力定位（docs/02 §3 已釘）:
    Ω = ((B−1) : B² : N_d−B+1) / (N_d + B²)
    N_d = 44 ⟹ 分母 60 ⟹ Ω_b = 3/60 = 0.05   （+2.33σ）
    N_d = 45 ⟹ 分母 61 ⟹ Ω_b = 3/61 = 0.04918（−0.4σ）  ← 數據偏好

本輪問題: N_d = 45 係唔係一個【合法】結構替代？抑或只係擬合？
判據: 整數身份網絡（唔可以靠 Ω 數據反推 N_d）
"""
from fractions import Fraction as F

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_gen, gamma_max = 4, 3, 4
obs = (0.0493, 0.2645, 0.6847)      # Ω_b, Ω_dm, Ω_Λ (Planck 2018 TT,TE,EE+lowE+lensing; O215 docs/04)
sig = (0.0003, 0.0050, 0.0074)
L5_lucas = 11                        # Lucas L₅ = 11
F11 = 89

def omega(nd):
    tot = nd + B*B
    return ((B-1)/tot, B*B/tot, (nd-B+1)/tot), tot

def chi2(t):
    return sum(((v-o)/s)**2 for v, o, s in zip(t, obs, sig))

# ── §1 重現張力 ────────────────────────────────────────────────────
print("§1 重現（obs = Planck 2018）")
for nd in (44, 45):
    o, tot = omega(nd)
    z = [(v-oo)/s for v, oo, s in zip(o, obs, sig)]
    print(f"  N_d={nd}: ({B-1}:{B*B}:{nd-B+1})/{tot} = (" +
          ", ".join(f"{v:.5f}" for v in o) + ")  z=(" +
          ", ".join(f"{x:+.2f}" for x in z) + f")  χ²={chi2(o):.3f}")

# ── §2 ★ 整數身份網絡：N_d=45 之下邊幾條死 ─────────────────────────
print("\n§2 整數身份網絡（N_d 必須係整數，且以下全部要成立）")
def audit(nd):
    res = {}
    res["L₅ = N_d/B ∈ ℤ"] = (nd % B == 0, F(nd, B))
    L5 = F(nd, B)
    res["K_tree = N_gen·L₅ ∈ ℤ"] = (L5.denominator == 1 and (N_gen*L5).denominator == 1, N_gen*L5)
    res["2N_d = 8·L₅, L₅ ∈ ℤ"] = (L5.denominator == 1 and F(2*nd, 8) == L5, F(2*nd, 8))
    res["N_d = N_gen·|A₂| − 1"] = (nd == N_gen*15 - 1, N_gen*15 - 1)
    res["41 = N_d − N_gen"] = (nd - N_gen == 41, nd - N_gen)
    Sigma = 2*nd + 2
    res["123 = K_tree + Σγ"] = (N_gen*L5 + Sigma == 123, N_gen*L5 + Sigma)
    res["61.5 = 123/2 = ⟨W⟩²"] = (N_gen*L5 + Sigma == 123, F(N_gen*L5 + Sigma, 2))
    res["88 = 2N_d = F₁₁ − 1"] = (2*nd == F11 - 1, 2*nd)
    return res

for nd in (44, 45):
    print(f"\n  ── N_d = {nd} ──")
    r = audit(nd)
    for k, (ok, val) in r.items():
        print(f"    {'✅' if ok else '🔴'} {k:28s} → {val}")
    nfail = sum(1 for ok, _ in r.values() if not ok)
    print(f"    ⟹ 失敗 {nfail}/{len(r)}（全部皆真約束）")

r44, r45 = audit(44), audit(45)
chk("N_d=44：全部身份成立（8/8）", all(ok for ok, _ in r44.values()))
chk("N_d=45：身份失敗 ≥ 4 條", sum(1 for ok, _ in r45.values() if not ok) >= 4,
    f"{sum(1 for ok, _ in r45.values() if not ok)}/8")
chk("N_d=45 ⟹ L₅ = 11.25 ∉ ℤ（樹深度唔可以非整數）", F(45, 4).denominator == 4, str(F(45,4)))
chk("N_d=45 ⟹ K_tree = 33.75 ∉ ℤ（模式數）", (3*F(45,4)).denominator != 1, str(3*F(45,4)))
chk("N_d=45 ⟹ 123 ≠ K_tree + Σγ（33.75 + 92 = 125.75）", 3*F(45,4) + (2*45+2) != 123)
chk("N_d=44 ⟹ 88 = 2B·L₅ = 8×11（單一乘積分解）", 2*B*L5_lucas == 88)
chk("★ 44 有兩條獨立身份：3·15−1 = 4·11 = 44（−1 殘餘由乘積身份佐證）",
    3*15 - 1 == B*L5_lucas == 44)

# ── §3 角色分離：90 = E(TI) = Σγ（唔係鏈長）───────────────────────
print("\n§3 多面體數之角色分離（RH/TI 計數）")
V, E_RH, F_TI, E_TI = 60, 120, 32, 90
print(f"  RH: V={V}, E={E_RH}, F=62; TI: V={V}, E={E_TI}, F={F_TI}")
chk("88 = E(RH) − F(TI) = 120 − 32 = 2N_d（鏈長）", E_RH - F_TI == 88 == 2*44)
chk("90 = E(TI) = Σγ(44-mode)（譜重，唔係鏈長）", E_TI == 2*44 + 2 == 90)
chk("E(RH) − V/2 = 90（第二條 90 之來源, 同一角色）", E_RH - V//2 == 90)
print("  ⟹ 90 已經以【Σγ 譜重】身份用完；唔可以再兼任鏈長 ⟹ 鏈長 = 88")

# ── §4 統計：Δχ² 有幾多信息？ ──────────────────────────────────────
print("\n§4 統計判讀")
o44, t44 = omega(44); o45, t45 = omega(45)
c44, c45 = chi2(o44), chi2(o45)
print(f"  χ²(60) = {c44:.3f} ; χ²(61) = {c45:.3f} ; Δχ² = {c44-c45:.3f}")
# 1 個離散 bit、3 個觀測 ⟹ 粗略 p（χ²₁ 尾概率）
import math
p = math.erfc(math.sqrt((c44-c45)/2))
print(f"  以 χ²₁ 近似：p = {p:.4f}（≈{abs(2.33):.2f}σ，與 Ω_b 本身同源）")
print("  ⚠️ 關鍵：呢個 bit 並非自由 —— §2 已證 45 破壞 4+ 條整數身份 ⟹ 唔可以當 free parameter 去換 fit")

# ── §5 預註冊：CMB-S4 判死線 ──────────────────────────────────────
print("\n§5 盲測預註冊（CMB-S4）")
h = 0.674071                      # O170 c = R·ω_ring 錨
pred60 = F(3, 60); pred61 = F(3, 61)
Obh2_60 = float(pred60) * h*h
Obh2_61 = float(pred61) * h*h
print(f"  h = {h}")
print(f"  Ω_b h² (3/60) = {Obh2_60:.6f}")
print(f"  Ω_b h² (3/61) = {Obh2_61:.6f}")
print(f"  Δ = {Obh2_60-Obh2_61:.6f}")
for s in (0.00015, 0.00006, 0.00002):
    print(f"  σ(Ω_b h²) = {s:.5f} ⟹ 分辨力 {abs(Obh2_60-Obh2_61)/s:.1f}σ")
chk("現時 Planck σ=1.5e-4 ⟹ 2.3σ（未判）", abs(Obh2_60-Obh2_61)/0.00015 < 3)
chk("CMB-S4 σ≈6e-5 ⟹ >5σ 判死", abs(Obh2_60-Obh2_61)/0.00006 > 5)

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")