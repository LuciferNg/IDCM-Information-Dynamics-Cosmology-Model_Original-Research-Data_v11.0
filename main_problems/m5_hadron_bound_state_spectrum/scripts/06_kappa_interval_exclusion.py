"""P16 docs/06 配套 — H4b: 態＋號依賴核之精確否證（κ-區間交集）
================================================================================
H4a 隨機搜索 0 命中只算探測；本輪決定性計算：
  E_H/κ ∈ S_H = { Σ_{i<j} q_i q_j s_i s_j 4^(−|ℓ_i−ℓ_j|) }，ℓ 逐態自由（= 態依賴），s_f = ±1（s_u=+1 規範）
  通道 k：obs_k = κ(E_A − E_B) ⟹ κ ∈ ∪_{a∈S_A,b∈S_B} [obs_k±σ_k]/(a−b)
  8 通道之 κ 准許集**無公共交集** ⟹ 核族（態依賴＋Z₂ 號，任意層）被否證（精確，非抽樣）
實現：Fraction 生成值集 → float 廣播算差值 → 去重 → 區間交集
輸出：`data/06_kappa_interval_exclusion.json`
"""
import json
from fractions import Fraction as F
from itertools import product
import numpy as np

B = 4
Q = {"u": F(2, 3), "d": F(-1, 3), "s": F(-1, 3), "c": F(2, 3),
     "ub": F(-2, 3), "db": F(1, 3), "sb": F(1, 3)}
CH = [("n−p", ("u", "d", "d"), ("u", "u", "d"), -1.2667, 0.0940),
      ("K⁺−K⁰", ("u", "sb"), ("d", "sb"), -1.3740, 0.0963),
      ("π⁺−π⁰", ("u", "db"), ("u", "ub"), 2.0336, 0.0940),
      ("Σ⁻−Σ⁰", ("d", "d", "s"), ("u", "d", "s"), 2.2470, 0.1010),
      ("Σ⁰−Σ⁺", ("u", "d", "s"), ("u", "u", "s"), 0.7120, 0.1196),
      ("Σ⁻−Σ⁺", ("d", "d", "s"), ("u", "u", "s"), 2.9590, 0.2028),
      ("Ξ⁻−Ξ⁰", ("d", "s", "s"), ("u", "s", "s"), 4.2900, 0.2301),
      ("D⁺−D⁰", ("c", "db"), ("c", "ub"), 2.2600, 0.1172)]
LMAX = 12
STATES = sorted({s for c in CH for s in (c[1], c[2])})

def fl(v):
    return v[:-1] if v.endswith("b") else v

def values(state, sgn):
    n = len(state)
    out = set()
    for lev in product(range(LMAX + 1), repeat=n):
        e = F(0)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = fl(state[i]), fl(state[j])
                e += Q[state[i]] * Q[state[j]] * sgn[a] * sgn[b] * F(1, B ** abs(lev[i] - lev[j]))
        out.add(e)
    return np.array(sorted(float(x) for x in out))

print("=" * 96)
print("P16 docs/06 — H4b: 態＋號依賴核精確否證（κ 區間交集）")
print("=" * 96)
print(f"獨立強子態 = {len(STATES)}；ℓ ∈ 0..{LMAX}；號組合 8（s_u = +1 規範）")

def interval_intersect(cur, iv):
    out = []
    for a1, b1 in cur:
        for a2, b2 in iv:
            lo, hi = max(a1, a2), min(b1, b2)
            if hi > lo:
                out.append((lo, hi))
    out.sort()
    merged = []
    for lo, hi in out:
        if merged and lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged[:2000]

def solve(obs_list, V, verbose=False):
    """obs_list = [(name, ini, fin, obs, sig)]；回傳 κ 區間 list"""
    cur = [(0.0, float("inf"))]
    for name, ini, fin, obs, sig in obs_list:
        d = (V[ini][:, None] - V[fin][None, :]).ravel()
        d = d[np.abs(d) > 1e-12]
        d = np.unique(np.round(d, 12))
        lo = (obs - sig) / d; hi = (obs + sig) / d
        iv = np.stack([np.minimum(lo, hi), np.maximum(lo, hi)], axis=1)
        iv = iv[iv[:, 1] > 0]
        iv = iv[np.argsort(iv[:, 0])]
        cur = interval_intersect(cur, [(float(x), float(y)) for x, y in iv])
        if not cur:
            break
    return cur

# ---------- 自測（pipeline 校驗）：用已知配置造 obs，必須搵返 κ_true ----------
rng = np.random.default_rng(7)
sgn_t = {"u": 1, "d": 1, "s": -1, "c": 1}
Vt = {st: values(st, sgn_t) for st in STATES}
lev_t = {st: tuple(int(x) for x in rng.integers(0, 4, size=len(st))) for st in STATES}

def e_state(st):
    n = len(st); e = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = fl(st[i]), fl(st[j])
            q = float(Q[st[i]] * Q[st[j]] * F(1, 1))
            e += q * sgn_t[a] * sgn_t[b] * B ** (-abs(lev_t[st][i] - lev_t[st][j]))
    return e
K_TRUE = 7.0
obs_t = [(n, i, f, K_TRUE * (e_state(i) - e_state(f)), 0.05) for n, i, f, o, s in CH]
cur_t = solve(obs_t, Vt)
print(f"[自測] 合成數據（κ_true = {K_TRUE}）⟹ κ 區間數 = {len(cur_t)}"
      + (f"，首個區間 = ({cur_t[0][0]:.4f}, {cur_t[0][1]:.4f})"
         + ("  ✅ 包含 κ_true" if cur_t[0][0] <= K_TRUE <= cur_t[0][1] else "  🔴 唔包含 ⟹ pipeline 有 bug")
         if cur_t else "  🔴 空 ⟹ pipeline 有 bug，否定結論作廢"))

res = []
for sd in (1, -1):
    for ss in (1, -1):
        for sc in (1, -1):
            sgn = {"u": 1, "d": sd, "s": ss, "c": sc}
            V = {st: values(st, sgn) for st in STATES}
            cur = solve(CH, V)
            res.append({"sgn": sgn, "intervals": len(cur),
                        "range": [[float(a), float(b)] for a, b in cur[:3]],
                        "sizes": {st: int(len(V[st])) for st in STATES}})
            print(f"  s=({sd:+d},{ss:+d},{sc:+d}): κ 區間數={len(cur)}"
                  + (f" 例 {[ (round(a,4), round(b,4)) for a,b in cur[:2]]}" if cur else "  ⟹ ∅"))

alive = [r for r in res if r["intervals"] > 0]
print(f"\n[判定] 存活號組合 = {len(alive)} / {len(res)}")
print("  → 空集 ⟹ 「態依賴層 + Z₂ 號」核族在 8 通道下**無解**（精確枚舉，非隨機）")
print("  → 有解 ⟹ 族過寬、不可否證，須改由 lattice Hamiltonian 導出 ψ")

out = {"LMAX": LMAX, "n_states": len(STATES), "alive": len(alive),
       "selftest_intervals": len(cur_t),
       "detail": [{"s": [int(r["sgn"][k]) for k in "udsc"], "intervals": int(r["intervals"]),
                   "range": [[float(a), float(b)] for a, b in r["range"]]} for r in res],
       "verdict": ("核族（態依賴＋號）8 通道無解" if not alive else "族過寬：存在解，不可否證")}
with open("data/06_kappa_interval_exclusion.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/06_kappa_interval_exclusion.json")