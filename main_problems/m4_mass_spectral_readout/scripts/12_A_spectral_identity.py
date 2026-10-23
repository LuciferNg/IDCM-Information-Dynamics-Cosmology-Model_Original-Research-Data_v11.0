#!/usr/bin/env python3
"""
M4 R4 — A 嘅譜身份（A = 1.64337 eV/γ²，唯一未閉合標度）(2026-10-24)
=====================================================================
上游（O257 R7 定案 / P6 docs）：m_ν = A·γ_n²，A = v²/M_comp×1e9 eV。
  A 係「譜讀數單位」；M_comp 依賴 = R7 唯一殘留 🔴（O66 1/76 已否證）。
歸一化鐵律：A 有兩個版本 ——
  v = 246.0     → A = 1.64044 eV/γ²（P6 docs/02,07 舊版）
  v_SM=246.2196 → A = 1.64337 eV/γ²（P6 docs/12 vSM_update，現行 IDCM 錨）
本輪一律用 v_SM 版本（跨版本引用先核歸一化）。

問：A 有冇零參數譜表達（唔經 M_comp）？
Protocol（事前定義，先寫低族再掃）：
  族 F = { anchor × γ_i^k × S }
    anchor ∈ 8 個已註冊標度：1 eV, m_e, Λ_IDCM, M₀, v_SM, √(M_Pl·Λ), M_Pl, ħω_ring(以 eV)
    k ∈ [-45,45]（γ₁ 級梯指數範圍，覆蓋 10^±115）
    S ∈ 18 個結構常數（B=4, N_d=44, L=88, F/L 數列, N_gen=3 …）
  判準：best dev < 0.1% 才進入候選；再用 density（同族對隨機目標）判斷可辨識性
  禁：見 A 反揀 anchor/k/S（族必須先列印，後掃）
輸出：data/12_A_identity.json
"""
import json
import numpy as np
import mpmath as mp

mp.mp.dps = 40
L = 88
g1 = mp.mpf(4)*mp.sin(mp.pi/L)**2

# ── 標度（eV）──
eV = mp.mpf(1)
m_e = mp.mpf('510998.950')
v_SM = mp.mpf('246.2196')*mp.mpf('1e9')          # eV
Lambda_IDCM = mp.mpf('45.28487')*mp.mpf('1e6')   # eV（O142 v16 機制值 45.28487 MeV）
M0 = mp.mpf('27.900915')*mp.mpf('1e6')           # eV
M_Pl = mp.mpf('1.2209e19')*mp.mpf('1e9')         # eV
hbar_omega_ring = mp.mpf('1.0')                  # placeholder（未註冊 eV 值 → 唔用）
anchors = {
    "1eV": eV,
    "m_e": m_e,
    "Λ_IDCM": Lambda_IDCM,
    "M₀": M0,
    "v_SM": v_SM,
    "√(M_Pl·Λ_IDCM)": mp.sqrt(M_Pl*Lambda_IDCM),
    "M_Pl": M_Pl,
    "Λ_IDCM²/M₀": Lambda_IDCM**2/M0,
}
S_pool = {
    "1": 1, "B=4": 4, "N_d=44": 44, "L=88": 88, "N_gen=3": 3, "F5=5": 5, "F7=13": 13,
    "F8=21": 21, "L5=11": 11, "L7=29": 29, "L9=76": 76, "N_gen²=9": 9, "L6=18": 18,
    "B·N_gen=12": 12, "L+N_d=132": 132, "N_d/2=22": 22, "L/2=44": 44, "√5": float(np.sqrt(5)),
}
ks = list(range(-45, 46))

A_old = mp.mpf('1.64044')
A_target = mp.mpf('1.64337')

print("="*80)
print("M4 R4: A 譜身份 — 零參數族掃描（先列族，後掃）")
print("="*80)
print(f"γ₁ = {mp.nstr(g1, 20)}")
print(f"目標 A = {A_target} eV/γ²  (v_SM=246.2196 版本; 舊版 1.64044 @ v=246.0)")
print("\n預先列出族（anchor × γ₁^k × S）:")
print(f"  anchors ({len(anchors)}): " + ", ".join(f"{k}={mp.nstr(v,6)}eV" for k, v in anchors.items()))
print(f"  k ∈ [{ks[0]}, {ks[-1]}]  ({len(ks)} 個)")
print(f"  S ({len(S_pool)}): " + ", ".join(S_pool))
print(f"  ⟹ 族大小 = {len(anchors)*len(ks)*len(S_pool)} 組合")

fam = []
for an, av in anchors.items():
    for k in ks:
        base = av*(g1**k)
        for sn, S in S_pool.items():
            fam.append((an, k, sn, base*mp.mpf(str(S))))
print(f"  族實際條目 = {len(fam)}")

def scan(target, thr=0.1):
    hits = [(an, k, sn, val) for an, k, sn, val in fam if abs(val/target - 1)*100 < thr]
    hits.sort(key=lambda z: abs(z[3]/target - 1))
    return hits

print("\n" + "-"*80)
print(f"[S1] A 目標 {A_target} eV，thr 0.1%")
print("-"*80)
h = scan(A_target, 0.1)
print(f"  hits = {len(h)}")
for an, k, sn, val in h[:10]:
    print(f"    {an:<16} γ₁^{k:<3} × {sn:<10} = {mp.nstr(val, 8):>14} eV  dev {float(abs(val/A_target-1)*100):.5f}%")
if not h:
    # 放寬看最接近
    best = min(fam, key=lambda z: abs(z[3]/A_target-1))
    print(f"  最近: {best[0]} γ₁^{best[1]} ×{best[2]} = {mp.nstr(best[3],8)} eV dev {float(abs(best[3]/A_target-1)*100):.3f}%")

# ── density：同一族對隨機目標 ──
print("\n" + "-"*80)
print("[S2] Density：同族對隨機目標（域 10^-3 … 10^9 eV，log-均勻）")
print("-"*80)
rng = np.random.default_rng(20261028)
NMC = 20000
thr_list = (0.1, 0.01)
cnt = {t: 0 for t in thr_list}
# 族值取 log（固定）
logfam = np.array([float(mp.log10(v)) for _, _, _, v in fam])
for _ in range(NMC):
    t = 10**rng.uniform(-3, 9)
    lt = np.log10(t)
    devmin = np.min(np.abs(10**(logfam-lt) - 1))*100
    for thr in thr_list:
        if devmin < thr:
            cnt[thr] += 1
for thr in thr_list:
    print(f"  thr {thr}%: 隨機目標命中率 = {cnt[thr]/NMC*100:.2f}%")
print("  ⟹ 命中率 >5% = 族密度足以解釋任何單點命中（唔可辨識）")

# ── M_comp 反解 + 結構比對 ──
print("\n" + "-"*80)
print("[S3] A 反解 M_comp（A = v²/M_comp × 1e9）→ 標度身份比對")
print("-"*80)
M_comp_GeV = (mp.mpf('246.2196')**2)/(A_target/mp.mpf('1e9'))
print(f"  M_comp = v²·1e9/A = {mp.nstr(M_comp_GeV, 10)} GeV")
refs = {"M_Pl·γ₁⁰": mp.mpf('1.2209e19'),
        "M_Pl·γ₁¹": mp.mpf('1.2209e19')*g1,
        "M_Pl·γ₁²": mp.mpf('1.2209e19')*g1**2,
        "M_Pl·γ₁³": mp.mpf('1.2209e19')*g1**3,
        "M_Pl·γ₁⁴": mp.mpf('1.2209e19')*g1**4}
for rn, rv in refs.items():
    q = M_comp_GeV/rv
    nint = mp.nint(q)
    if nint == 0:
        print(f"  M_comp/{rn} = {mp.nstr(q, 8)}  (太小，無整數身份)")
        continue
    print(f"  M_comp/{rn} = {mp.nstr(q, 8)}  (最近整數 {mp.nstr(nint,6)}, dev {float(abs(q/nint-1)*100):.3f}%)")
print("  ⟹ 若比值只係「近整數」而無結構身份（如 N_d、L、B^k），仍屬未閉合")

out = {
    "A_target_eV": str(A_target), "A_old_v246": str(A_old),
    "family_size": len(fam),
    "hits_0.1pct": [[an, k, sn, mp.nstr(val, 12)] for an, k, sn, val in h],
    "density_pct": {str(t): cnt[t]/NMC*100 for t in thr_list},
    "M_comp_GeV": mp.nstr(M_comp_GeV, 12),
}
with open("data/12_A_identity.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\n→ data/12_A_identity.json")
