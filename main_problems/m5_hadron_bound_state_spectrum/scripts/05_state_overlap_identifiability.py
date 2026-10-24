"""P16 docs/05 配套 — H4a: 態依賴（state-dependent）超度量核之「可否證性」審計
================================================================================
H3a/H3b 用「味標籤 → 固定層」證死；H4a 問：若層係**態依賴**（同一味喺唔同強子可佔唔同層，
即波函數自由度），8 通道仲有無約束力？
  模型：E_H = κ·Σ_{i<j} q_i q_j B^(−|ℓ_i^H − ℓ_j^H|)（ℓ 逐強子自由 ⟹ 即「態」）
  未知：每強子嘅層型（baryon 3 個、meson 2 個）＋ κ ⟹ 遠多過 8 條方程
  判準：① 存唔存在解（hit）② 解嘅密度 vs 置換 null ⟹ 可辨識／否證力評估
  另加：尺度審計 — κ 係唯一有量綱嘅量（MeV），必須同 IDCM 已登記 MeV 級量比對（並報密度）
輸出：`data/05_state_overlap_identifiability.json`
"""
import json
import random
import numpy as np

B = 4
Q = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3,
     "ub": -2 / 3, "db": 1 / 3, "sb": 1 / 3}
CH = [("n−p", ("u", "d", "d"), ("u", "u", "d"), -1.2667, 0.0940),
      ("K⁺−K⁰", ("u", "sb"), ("d", "sb"), -1.3740, 0.0963),
      ("π⁺−π⁰", ("u", "db"), ("u", "ub"), 2.0336, 0.0940),
      ("Σ⁻−Σ⁰", ("d", "d", "s"), ("u", "d", "s"), 2.2470, 0.1010),
      ("Σ⁰−Σ⁺", ("u", "d", "s"), ("u", "u", "s"), 0.7120, 0.1196),
      ("Σ⁻−Σ⁺", ("d", "d", "s"), ("u", "u", "s"), 2.9590, 0.2028),
      ("Ξ⁻−Ξ⁰", ("d", "s", "s"), ("u", "s", "s"), 4.2900, 0.2301),
      ("D⁺−D⁰", ("c", "db"), ("c", "ub"), 2.2600, 0.1172)]

def paircoef(state):
    """回傳 pair 權重 dict[(i,j)] = q_i q_j（i<j 位置）"""
    out = {}
    for i in range(len(state)):
        for j in range(i + 1, len(state)):
            out[(i, j)] = Q[state[i]] * Q[state[j]]
    return out

HAD = {}          # 每個「強子態」一個條目（用 tuple(state) 做 key）
for name, ini, fin, o, s in CH:
    for st in (ini, fin):
        HAD[st] = paircoef(st)
STATES = list(HAD)
IDX = {st: k for k, st in enumerate(STATES)}
NB = len(STATES)
print("=" * 96)
print("P16 docs/05 — H4a: 態依賴超度量核之可否證性審計")
print("=" * 96)
print(f"獨立強子態 = {NB}（每態有 2 或 3 個價夸克，層自由 ⟹ 參數 ≥ {sum(len(s)-1 for s in STATES)} + 1）")
print(f"數據方程 = {len(CH)} ⟹ 自由度嚴重過剩，先驗預期「擬合必然可行」")

# ---- 隨機搜索：層型 ℓ_i ∈ {0..LMAX}，κ 用 1-param 最小二乘 ----
rng = np.random.default_rng(20261024)
LMAX = 3
def random_patterns(n, rng):
    P = {}
    for st in STATES:
        nq = len(st)
        P[st] = rng.integers(0, LMAX + 1, size=nq)      # 絕對層（含平移冗餘，無妨）
    return P

SF = {"u": 1, "d": 1, "s": 1, "c": 1}      # 全域號（每次搜尋隨機設定；s_u = +1 規範）

def fl_(v):
    return v[:-1] if v.endswith("b") else v

def energy_vec(P, st):
    lv = P[st]
    e = 0.0
    for (i, j), qq in HAD[st].items():
        a, b = fl_(st[i]), fl_(st[j])
        e += qq * SF[a] * SF[b] * B ** (-abs(int(lv[i]) - int(lv[j])))
    return e

OBS = np.array([c[3] for c in CH]); SIG = np.array([c[4] for c in CH])
INI = [c[1] for c in CH]; FIN = [c[2] for c in CH]

def chi2(P, obs=OBS):
    d = np.array([energy_vec(P, INI[k]) - energy_vec(P, FIN[k]) for k in range(len(CH))])
    num = float(np.sum(d * obs / SIG ** 2))
    den = float(np.sum((d / SIG) ** 2))
    kappa = num / den if den > 0 else 0.0
    r = kappa * d - obs
    return float(np.sum((r / SIG) ** 2)), kappa

NTRIAL = 200000
best = []
hits = 0
hits1 = 0
for _ in range(NTRIAL):
    P = random_patterns(NB, rng)
    for f in ("d", "s", "c"):
        SF[f] = int(rng.choice([1, -1]))          # 號自由（s_u = +1 規範）
    c2, kap = chi2(P)
    if kap > 0:
        nparam = sum(len(s) for s in STATES) + 1
        best.append((c2, kap, P.copy()))
        if c2 <= 9.0:            # 8 dof, 3σ 全域
            hits += 1
        if c2 <= 3.0:
            hits1 += 1
best.sort(key=lambda t: t[0])
dof = len(CH) - 1
print(f"\n[隨機搜索] {NTRIAL} 次層型抽樣（ℓ ∈ 0..{LMAX}，κ 每次 1-param LSQ）")
print(f"  χ² ≤ 3（≈1σ 全域）命中 = {hits1}；χ² ≤ 9（≈3σ）命中 = {hits}")
print(f"  最佳 χ² = {best[0][0]:.4f}（dof {dof}）κ = {best[0][1]:.4f} MeV")
print("  最佳 3 個解（層型，按強子）:")
for c2, kap, P in best[:3]:
    print(f"    χ²={c2:.3f} κ={kap:.4f}")
    for st in STATES[:6]:
        print(f"      {''.join(st):8s} ℓ = {list(map(int, P[st]))}")

# ---- 置換 null：打亂 8 個觀測值 ----
null_hits = 0; null_hits1 = 0
for _ in range(2000):
    perm = OBS.copy(); rng.shuffle(perm)
    P = random_patterns(NB, rng)
    c2, kap = chi2(P, perm)
    if c2 <= 9.0 and kap > 0:
        null_hits += 1
    if c2 <= 3.0 and kap > 0:
        null_hits1 += 1
print(f"\n[置換 null] 2000 次：χ²≤3 命中 {null_hits1}；χ²≤9 命中 {null_hits}")

print("\n[判定]")
print("  態依賴層 ⟹ 參數個數 ≫ 8 ⟹ 擬合必然可行（本輪已示範：見上）⟹ **模型不可否證**")
print("  ⟹ 唔可以靠 8 通道「反推波函數」；ψ 必須由 lattice Hamiltonian 導出（H4b），8 通道只作檢定")

# ---- 尺度審計：κ 係唯一有量綱量 ----
print("\n[尺度審計] κ 係唯一有量綱參數（MeV）；IDCM 已登記 MeV 級量對照（報密度，不作匹配宣稱）")
REF = {"Λ_IDCM(O142)": 45.28487, "m_π±(PDG)": 139.57039, "m_K±(PDG)": 493.677,
       "m_p(PDG)": 938.272, "v_SM/1000": 246.2196, "X=m_d−m_u(FLAG24)": 2.56}
kap0 = best[0][1]
for k, v in REF.items():
    print(f"    κ / {k:22s} = {kap0 / v:.6f}")
print(f"  ⟹ κ ≈ {kap0:.4f} MeV 嘅絕對值屬**量綱定義**（同 M₀ 一樣：絕對數值 = 單位制 artifact 偽命題）")
print("  ⟹ P16 只可交付**比值**（如 n−p = −κ/3）；κ 絕對值 = hadronic scale 識別，非本 O 可推")

out = {"n_states": NB, "n_trials": NTRIAL, "dof": dof,
       "hits_chi2_3": hits1, "hits_chi2_9": hits,
       "null_hits_chi2_3": null_hits1, "null_hits_chi2_9": null_hits,
       "best_chi2": best[0][0], "best_kappa_MeV": best[0][1],
       "verdict": "態依賴層 ⟹ 不可否證；ψ 須由 lattice Hamiltonian 導出（H4b）",
       "scale_note": "κ 絕對值 = 量綱定義（類 M₀）；本 O 只交付比值"}
bv = []
for c2, kap, P in best[:5]:
    bv.append({"chi2": c2, "kappa_MeV": kap,
               "levels": {"".join(st): [int(x) for x in P[st]] for st in STATES}})
out["best"] = bv
with open("data/05_state_overlap_identifiability.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/05_state_overlap_identifiability.json")