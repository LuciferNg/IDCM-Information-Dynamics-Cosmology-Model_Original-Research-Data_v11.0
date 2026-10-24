"""P16 docs/04 配套 — H3b: 含號（Z₂）超度量核之可辨識性檢定 + 置換 null
================================================================================
H3a 否證正核類（K 通道符號矛盾）⟹ 出路 = 核帶 Z₂ 號：x_AB = s_A s_B B^(−|ℓ_A−ℓ_B|)，s_f = ±1。
本輪：含號核類**能否**用 8 通道定住 (ℓ_f, s_f)？（只檢定，不擬合）
  參數：ℓ_d, ℓ_s, ℓ_c ∈ 0..33（ℓ_u = 0 平移規範）；s_d, s_s, s_c = ±1（s_u = +1 規範）
        ⟹ 39 304 × 8 = 314 432 組合；κ 由 n−p 結構關係 −κ/3 = −1.2667 釘死（非擬合）
  判準：8 通道全部 |z| ≤ {1,2,3} 才算命中；命中數必須同**置換 null**（200 次重排觀測）比
向量化：每通道 E/κ = c + Σ_p a_p x_p（Fraction 精確算係數），x 由 (ℓ, s) 生成 ⟹ numpy 全枚舉
輸出：`data/04_signed_kernel_identifiability.json`
"""
import json
import random
from fractions import Fraction as F
import numpy as np

B = 4
QC = {"u": F(2, 3), "d": F(-1, 3), "s": F(-1, 3), "c": F(2, 3),
      "ub": F(-2, 3), "db": F(1, 3), "sb": F(1, 3)}
CH = [("n−p", ("u", "d", "d"), ("u", "u", "d"), -1.2667, 0.0940),
      ("K⁺−K⁰", ("u", "sb"), ("d", "sb"), -1.3740, 0.0963),
      ("π⁺−π⁰", ("u", "db"), ("u", "ub"), 2.0336, 0.0940),
      ("Σ⁻−Σ⁰", ("d", "d", "s"), ("u", "d", "s"), 2.2470, 0.1010),
      ("Σ⁰−Σ⁺", ("u", "d", "s"), ("u", "u", "s"), 0.7120, 0.1196),
      ("Σ⁻−Σ⁺", ("d", "d", "s"), ("u", "u", "s"), 2.9590, 0.2028),
      ("Ξ⁻−Ξ⁰", ("d", "s", "s"), ("u", "s", "s"), 4.2900, 0.2301),
      ("D⁺−D⁰", ("c", "db"), ("c", "ub"), 2.2600, 0.1172)]
PAIRS = [("d", "u"), ("s", "u"), ("d", "s"), ("c", "u"), ("c", "d"), ("c", "s")]   # sorted 形式（配合 sorted(())）
KAPPA = 3.0 * 1.2667

def fl(v):
    return v[:-1] if v.endswith("b") else v

def coefs(state):
    c = F(0); a = {p: F(0) for p in PAIRS}
    for i in range(len(state)):
        for j in range(i + 1, len(state)):
            x, y = fl(state[i]), fl(state[j])
            q = QC[state[i]] * QC[state[j]]
            if x == y:
                c += q
            else:
                a[tuple(sorted((x, y)))] += q
    return c, a

CO = []
for name, i, f, o, s in CH:
    c0, a0 = coefs(i); c1, a1 = coefs(f)
    CO.append((name, float(c0 - c1), np.array([float(a0[p] - a1[p]) for p in PAIRS]), o, s))
print("=" * 96)
print("P16 docs/04 — H3b: 含號超度量核可辨識性（314 432 組合 ＋ 置換 null 200 次）")
print("=" * 96)
print(f"κ = {KAPPA:.4f} MeV（由 n−p 結構關係 −κ/3 = −1.2667 定死）")

# ---- 生成組合矩陣：E/κ = c + Σ a_p s_A s_B 4^(−d) ----
grid = []
for ld in range(34):
    for ls in range(34):
        for lc in range(34):
            d = {("d", "u"): abs(ld), ("s", "u"): abs(ls), ("d", "s"): abs(ls - ld),
                 ("c", "u"): abs(lc), ("c", "d"): abs(lc - ld), ("c", "s"): abs(lc - ls)}
            grid.append((ld, ls, lc, d))
print(f"層指派組合 = {len(grid)}；× 8 號組合 = {len(grid)*8}")

lev_arr = np.array([g[:3] for g in grid])
dx = {p: np.array([g[3][p] for g in grid], dtype=float) for p in PAIRS}
X = {p: (B ** (-dx[p])) for p in PAIRS}

def eval_all(sgn):
    """sgn = (s_d, s_s, s_c)；回傳 (n_grid, 8) 殘差矩陣"""
    sd, ss, sc = sgn
    S = {("d", "u"): sd * 1, ("s", "u"): ss * 1, ("d", "s"): sd * ss,
         ("c", "u"): sc * 1, ("c", "d"): sc * sd, ("c", "s"): sc * ss}
    Xs = {p: X[p] * S[p] for p in PAIRS}
    cols = []
    for name, c, a, o, s in CO:
        e = c + sum(a[k] * Xs[PAIRS[k]] for k in range(len(PAIRS)))
        cols.append(KAPPA * e)
    return np.stack(cols, axis=1)

OBS = np.array([c[3] for c in CO]); SIG = np.array([c[4] for c in CO])

def count_hits(obs):
    hits = {1: 0, 2: 0, 3: 0}; best = []
    for sgn in [(sd, ss, sc) for sd in (1, -1) for ss in (1, -1) for sc in (1, -1)]:
        M = eval_all(sgn)
        z = np.abs((M - obs) / SIG).max(axis=1)
        for t in (1, 2, 3):
            hits[t] += int((z <= t).sum())
        k = int(np.argmin(z))            # 不論門檻，一律記錄最佳（用於報告最小偏差）
        best.append((float(z[k]),) + tuple(int(v) for v in lev_arr[k]) + sgn)
    best.sort()
    return hits, best

hits, best = count_hits(OBS)
print(f"\n[真數據] 命中（8 通道全部 |z| ≤ 門檻）：{hits}")
print("  最佳 10 個（z, ℓ_d, ℓ_s, ℓ_c, s_d, s_s, s_c）：")
for b in best[:10]:
    print(f"    z={b[0]:6.2f}  ℓ=({b[1]},{b[2]},{b[3]})  s=({b[4]:+d},{b[5]:+d},{b[6]:+d})")

rng = random.Random(20261024)
null = {1: 0, 2: 0, 3: 0}
NRUN = 200
for _ in range(NRUN):
    perm = list(OBS); rng.shuffle(perm)
    h, _ = count_hits(np.array(perm))
    for t in null:
        null[t] += h[t]
print(f"\n[置換 null] {NRUN} 次：總假命中 {null} ⟹ 每次平均 "
      f"{null[1]/NRUN:.2f} / {null[2]/NRUN:.2f} / {null[3]/NRUN:.2f}")

print("\n[判定]")
verdict = {}
for t in (1, 2, 3):
    nm = null[t] / NRUN
    if hits[t] == 0 and nm < 0.05:
        tag, ok = "🔴 核類否證（零命中；null 亦零 ⟹ 檢定有力）", False
    elif hits[t] > 0 and nm >= 1:
        tag, ok = "🟡 族過密（真命中不顯著於 null）", False
    elif hits[t] == 1 and nm < 0.2:
        tag, ok = "✅ 可辨識候選（唯一命中且 null 極低；須樣本外盲測）", True
    else:
        tag, ok = f"🟡 混合（真 {hits[t]} vs null {nm:.2f}）", False
    verdict[t] = {"true_hits": hits[t], "null_mean": nm, "identifiable": bool(ok), "tag": tag}
    print(f"  |z|≤{t}: 真 {hits[t]} vs null {nm:.2f} ⟹ {tag}")

out = {"kappa_MeV": KAPPA, "combos": len(grid) * 8,
       "hits_true": {str(k): v for k, v in hits.items()},
       "null_total": {str(k): v for k, v in null.items()},
       "null_mean": {str(k): v / NRUN for k, v in null.items()},
       "verdict": {str(k): v for k, v in verdict.items()},
       "best": [{"z": b[0], "ld": b[1], "ls": b[2], "lc": b[3],
                 "sd": b[4], "ss": b[5], "sc": b[6]} for b in best[:20]]}
with open("data/04_signed_kernel_identifiability.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/04_signed_kernel_identifiability.json")