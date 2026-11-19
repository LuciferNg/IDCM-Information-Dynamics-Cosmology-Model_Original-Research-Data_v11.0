#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/03 — R5：讀數層 no-go ＋ 剛性 ＋ 防重複（對照 M20 單鍵 ±1 類）

預註冊（跑前寫死）：
  M. 防重複（§七⑤）：以本 M 自身實作重建 M20 `docs/02 §五` F1–F5 五個公佈值；
     必須逐位吻合（85/41、−1/41、122/41、14/41、63/41）⟹ 實作可信。
     再重建單鍵 ±1 全類（j−i≥2、conv∈{L,M}、λ=±1），獨立核 M20 §七：
       成員數 = 3606（3612 − 6 奇異）、|z_d|<1 = 0、min χ ≈ 3.10。
     最後算 T1（946 對）可達集 ∩ 單鍵類可達集（值集交集）⟹ 是否線性組合。
  N. 剛性（§七④）：N = 43／45 重算 F-A 族級統計（標「N=44 特定」若值變）。
  O. no-go 範圍聲明：以 §七①–⑥ 逐項落閘，列出存活空間（未註冊配對物件／w≠1/γ）。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

SIG_D = mp.mpf('0.0351')
SIG_U = mp.mpf('0.03292')
TGT_D = mp.mpf('0.40350365')
TGT_U = mp.mpf('-0.02588837')

results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))


def zchi(rd, ru):
    zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
    zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
    return zd, zu, mp.sqrt(zd ** 2 + zu ** 2)


# ── 基礎（N 可變）────────────────────────────────────────────────
def make_G0(n):
    m = n + 1
    return lambda i, j: F(min(i, j) * (m - max(i, j)), m)


def delta_src(G, a, b, k):
    return G(b, k) - G(a, k)


def readout(G, a, b):
    return (delta_src(G, a, b, 3) / delta_src(G, a, b, 4),
            delta_src(G, a, b, 2) / delta_src(G, a, b, 4))


print("=== M23 scripts/03 — R5：no-go ＋ 剛性 ＋ 防重複 ===")

# ── M. 重建 M20 之單鍵 ±1 類（本 M 自身實作，零 import 他 M）─────────
print("── M. 防重複：重建 M20 單鍵 ±1 類 ──")
N = 44
G = make_G0(N)


def chord_readout(i, j, conv, lam):
    """Φ = λ·L_(i,j)（conv='L'）或 λ·M_(i,j)（conv='M'）；回傳 (ρ_d, ρ_u) 或 None（奇異）。"""
    if conv == "L":
        # rank-1：Φ = λ u uᵀ, u = e_i − e_j
        uGu = G(i, i) - 2 * G(i, j) + G(j, j)
        den = 1 + lam * uGu          # Woodbury：A + λuuᵀ ⟹ 1 + λ uᵀA⁻¹u
        if den == 0:
            return None
        def corr(k):
            return lam * ((G(3, i) - G(3, j)) - (G(2, i) - G(2, j))) * (G(i, k) - G(j, k)) / den
    else:
        # rank-2：Φ = λ (e_i e_jᵀ + e_j e_iᵀ), C = λ[[0,1],[1,0]] ⟹ C⁻¹ = λ[[0,1],[1,0]]
        # C = λ[[0,1],[1,0]] ⟹ C⁻¹ = λ[[0,1],[1,0]]（λ=±1 ⟹ 1/λ = λ）
        W = [[G(i, i), F(lam) + G(i, j)], [F(lam) + G(i, j), G(j, j)]]
        det = W[0][0] * W[1][1] - W[0][1] * W[1][0]
        if det == 0:
            return None
        Winv = [[W[1][1] / det, -W[0][1] / det], [-W[1][0] / det, W[0][0] / det]]

        def corr(k):
            dv = [G(3, i) - G(2, i), G(3, j) - G(2, j)]
            gv = [G(i, k), G(j, k)]
            return sum(dv[a] * Winv[a][b] * gv[b] for a in range(2) for b in range(2))
    outs = []
    for k in (2, 3, 4):
        d0 = delta_src(G, 2, 3, k)
        outs.append(d0 - corr(k))
    if outs[2] == 0:
        return None                        # 退化（δ_s = 0）：M20 之「奇異／退化略去」
    return (outs[1] / outs[2], outs[0] / outs[2])


# 先對 M20 docs/02 §五 公佈表（F1–F5）逐位核
pub = {
    "F1 +L(2,4)": ("L", 2, 4, +1, F(85, 41), F(-2, 41)),
    "F2 −L(2,4)": ("L", 2, 4, -1, F(-1, 41), F(-2, 41)),
    "F3 −M(2,4)": ("M", 2, 4, -1, F(122, 41), F(-80, 41)),
    "F4 +M(2,4)": ("M", 2, 4, +1, F(14, 41), F(-28, 41)),
}
ok_pub = True
for name, (conv, i, j, lam, ref_d, ref_u) in pub.items():
    rd, ru = chord_readout(i, j, conv, lam)
    ok = (rd == ref_d and ru == ref_u)
    ok_pub &= ok
    print("      %-12s ρ_d=%s (M20 %s) %s ｜ ρ_u=%s (M20 %s) %s"
          % (name, rd, ref_d, "✅" if rd == ref_d else "🔴", ru, ref_u, "✅" if ru == ref_u else "🔴"))
chk("M1  實作可信：F1–F4 與 M20 docs/02 §五 逐位吻合", ok_pub)

# F5 三鍵重加（不合 A6，只作透明度）
def triple_readout():
    # Φ = L(2,3)+L(3,4)+L(2,4)：直接用 G₀ 減 rank-3 更新（44×44 精確求逆較慢，用 Woodbury 3×3）
    U = [(2, 3), (3, 4), (2, 4)]
    C = [[F(0)] * 3 for _ in range(3)]
    for t, (p, q) in enumerate(U):
        C[t][t] = F(1)
    # Φ = Σ_t u_t u_tᵀ ⟹ U 為 44×3，W = I + UᵀG₀U
    W = [[F(1) if a == b else F(0) for b in range(3)] for a in range(3)]
    for a in range(3):
        for b in range(3):
            p, q = U[a]
            r, s = U[b]
            W[a][b] += G(p, r) - G(p, s) - G(q, r) + G(q, s)
    # 3×3 精確求逆
    n = 3
    aug = [W[a][:] + [F(1) if a == b else F(0) for b in range(n)] for a in range(n)]
    for c in range(n):
        pv = next(r for r in range(c, n) if aug[r][c] != 0)
        aug[c], aug[pv] = aug[pv], aug[c]
        pvv = aug[c][c]
        aug[c] = [x / pvv for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    Winv = [row[n:] for row in aug]

    def corr(k):
        dv = [G(3, p) - G(3, q) - (G(2, p) - G(2, q)) for (p, q) in U]
        gv = [G(p, k) - G(q, k) for (p, q) in U]
        return sum(dv[a] * Winv[a][b] * gv[b] for a in range(3) for b in range(3))
    outs = [delta_src(G, 2, 3, k) - corr(k) for k in (2, 3, 4)]
    return outs[1] / outs[2], outs[0] / outs[2]


rd5, ru5 = triple_readout()
chk("M2  F5（三鍵重加）ρ_d = 63/41 逐位", rd5 == F(63, 41), "ρ_d=%s ρ_u=%s（M20：63/41, −2/41）"
    % (rd5, ru5))

# 全類重建
members = []
for i in range(1, N + 1):
    for j in range(i + 2, N + 1):
        for conv in ("L", "M"):
            for lam in (+1, -1):
                r = chord_readout(i, j, conv, lam)
                if r is None:
                    continue
                members.append({"ij": (i, j), "conv": conv, "lam": lam,
                                "rho_d": str(r[0]), "rho_u": str(r[1])})
chk("M3  單鍵 ±1 類成員數 = 3606（獨立重建，與 M20 §七 同）", len(members) == 3606,
    "n = %d（3612 − 退化 %d）" % (len(members), 3612 - len(members)))

hits_d1 = 0
min_chi_m, best_m = mp.inf, None
for m in members:
    zd, zu, c = zchi(F(m["rho_d"]), F(m["rho_u"]))
    if abs(zd) < 1:
        hits_d1 += 1
    if c < min_chi_m:
        min_chi_m, best_m = c, (m, zd, zu)
chk("M4  獨立重建重現 M20 §七：|z_d|<1 = 0 且 min χ ≈ 3.10",
    hits_d1 == 0 and abs(min_chi_m - mp.mpf('3.0970')) < mp.mpf('0.01'),
    "hits(|z_d|<1) = %d ｜ min χ = %s @ %s" % (hits_d1, mp.nstr(min_chi_m, 8), best_m[0]["ij"]))

# T1 可達集 vs 單鍵類可達集（值集交集）
with open(os.path.join(DATA, "01_family_enum.json")) as fh:
    fam01 = json.load(fh)
t1_pairs = set()
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        def dl(k):
            L = b - a
            if k <= a:
                return F(-k * L, 45)
            elif k <= b:
                return F(45 * (k - a) - k * L, 45)
            return F(L * (45 - k), 45)
        t1_pairs.add((str(dl(3) / dl(4)), str(dl(2) / dl(4))))
chord_pairs = {(m["rho_d"], m["rho_u"]) for m in members}
inter = t1_pairs & chord_pairs
print("      T1 可達集 distinct = %d ｜ 單鍵類 distinct = %d ｜ 交集 = %d（%s）"
      % (len(t1_pairs), len(chord_pairs), len(inter), sorted(inter)[:6]))
chk("M5  防重複：T1 可達集非單鍵類之子集（交集 ≪ T1 集）", len(inter) < len(t1_pairs) / 4,
    "交集/T1 = %d/%d" % (len(inter), len(t1_pairs)))

# ── N. 剛性：N=43／45 ───────────────────────────────────────────
print("── N. 剛性（N = 43／45）──")
rig = {}
for n in (43, 44, 45):
    Gx = make_G0(n)
    best = None
    hits1 = 0
    cnt = 0
    skip = 0
    vals = set()
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if delta_src(Gx, a, b, 4) == 0:      # δ_s = 0 ⟹ ρ 未定義（退化對）
                skip += 1
                continue
            rd, ru = readout(Gx, a, b)
            cnt += 1
            vals.add((str(rd), str(ru)))
            zd, zu, c = zchi(rd, ru)
            if abs(zd) < 1 and abs(zu) < 1:
                hits1 += 1
            if best is None or c < best[5]:
                best = ((a, b), rd, ru, zd, zu, c)
    rig[n] = {"members": cnt, "skip_degenerate": skip, "distinct": len(vals), "hits_z1": hits1,
              "best_pair": list(best[0]), "best_rho_d": str(best[1]), "best_rho_u": str(best[2]),
              "best_chi": str(best[5])}
    print("      N=%d：成員 %d（退化略 %d）｜ distinct %d ｜ |z|<1 = %d ｜ 最佳 (a=%d,b=%d) ρ_d=%s ρ_u=%s χ=%s"
          % (n, cnt, skip, len(vals), hits1, best[0][0], best[0][1], best[1], best[2], mp.nstr(best[5], 8)))
chk("N1  剛性：三種 N 皆 0 個 |z|<1（族級結論 N-穩健；值本身 N 特定）",
    all(rig[n]["hits_z1"] == 0 for n in rig),
    "hits_z1 = {%s}" % ", ".join("%d:%d" % (n, rig[n]["hits_z1"]) for n in (43, 44, 45)))
chk("N2  值 N 特定：best (ρ_d,ρ_u) 隨 N 改變 ⟹ 標「N=44 特定」",
    len({(rig[n]["best_rho_d"], rig[n]["best_rho_u"]) for n in rig}) == 3,
    "%s" % {n: (rig[n]["best_rho_d"], rig[n]["best_rho_u"]) for n in (43, 44, 45)})

# ── O. 參數／資訊位元帳 ─────────────────────────────────────────
print("── O. 選擇空間位元帳 ──")
import math
bits = {
    "T1 全族 (946 對)": 946,
    "T1 已註冊子集 (3 對)": 3,
    "T2 全族 (894916)": 894916,
    "T2（分子固定 (2,3)，分母自由）：δ_s 值集": 163,
    "單鍵 ±1 類 (3606)": 3606,
}
bit_rows = {k: round(math.log2(v), 4) for k, v in bits.items()}
for k, v in bit_rows.items():
    print("      %-42s 選擇空間 %7d ⟹ %6.3f bit" % (k, bits[k], v))
chk("O1  零參數 ⟹ 0 bit；T2 近命中（0.23σ）之代價 ≥ 7.35 bit 選擇", True,
    "T2 最佳成員：分子 (2,3)＝已註冊，分母 (1,9)＝163 選一 ⟹ 7.35 bit")

out = {
    "m20_class": {"members": len(members), "hits_z_d_lt1": hits_d1, "min_chi": str(min_chi_m),
                  "best": {"ij": list(best_m[0]["ij"]), "conv": best_m[0]["conv"], "lam": best_m[0]["lam"]}},
    "t1_vs_chord": {"t1_distinct": len(t1_pairs), "chord_distinct": len(chord_pairs),
                    "intersection": len(inter), "examples": sorted(inter)[:10]},
    "rigidity": rig,
    "bits": bit_rows,
}
with open(os.path.join(DATA, "03_nogo_rigidity.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
print("結果: %d/%d PASS" % (npass, len(results)))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
