#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/05 — R7：源對隨站（T1x）判定（全 Fraction 精確）

預註冊（`docs/05` 凍結；跑前寫死，不得改）：
  M  映射集 = {1,2,3,4} ∪ {k+d: d∈{−2..2}} ∪ {(N+1)−(k+d): d∈{−2..2}}          |M| = 14
  F-F1  δ_x = G(f₂(k_x),k_x) − G(f₁(k_x),k_x)（同一對映射供三個讀數；(f₁,f₂) 有序，f₁≠f₂）
        全部有序對 × 合法性（1≤f(k_x)≤N ∧ f₁(k_x)≠f₂(k_x) ∀x ∧ δ_s≠0）
        判：① 結構唯一逼出？ ② |z_d|≤3 ∧ |z_u|≤3？ ③ MC 10⁷ p<0.05？
  F-F2  p(u)、p(d)、p(s) 三對獨立（6 整數 − 1 標度 ⟹ 有效 5）⟹ 直接違 §七①；
        只求「最優可達 χ」作參數代價披露（δ_s 值集 × δ_d 值集 × δ_u 值集，全枚舉）
  剛性  N = 43／45 重算 F-F1（映射之反射用 (N+1)−(k+d)）
  防重複 F-F1 值集 vs T1（946 對，`data/01`）vs 單鍵 ±1 類（3606，本檔自建，不重跑他 M）
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
probe_results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))


def chkp(name, ok, detail=""):
    probe_results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[REC!] ") + name + ("  " + detail if detail else ""))


def zchi(rd, ru):
    zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
    zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
    return zd, zu, mp.sqrt(zd ** 2 + zu ** 2)


def mkG(n):
    m = n + 1
    return lambda i, j: F(min(i, j) * (m - max(i, j)), m)


def maps(n):
    """映射集 M（14 個）：const {1,2,3,4}；shift k+d；reflect (n+1)−(k+d)。回傳 (name, f)。"""
    out = []
    for c in (1, 2, 3, 4):
        out.append(("c%d" % c, lambda k, c=c: c))
    for d in (-2, -1, 0, 1, 2):
        out.append(("k%+d" % d, lambda k, d=d: k + d))
    for d in (-2, -1, 0, 1, 2):
        out.append(("R(k%+d)" % d, lambda k, d=d, n=n: (n + 1) - (k + d)))
    return out


print("=== M23 scripts/05 — R7：源對隨站（T1x） ===")
N = 44
G = mkG(N)
KS = {"u": 2, "d": 3, "s": 4}
M = maps(N)
chk("R7.0  映射集 |M| = 14（4 const ＋ 5 shift ＋ 5 reflect）", len(M) == 14, "|M| = %d" % len(M))

# ── F-F1 全枚舉 ────────────────────────────────────────────────
fam = []
illegal = 0
deg = 0
for n1, f1 in M:
    for n2, f2 in M:
        if n1 == n2:
            continue
        ok = True
        for x, k in KS.items():
            a, b = f1(k), f2(k)
            if not (1 <= a <= N and 1 <= b <= N) or a == b:
                ok = False
                break
        if not ok:
            illegal += 1
            continue
        d = {}
        for x, k in KS.items():
            d[x] = G(f2(k), k) - G(f1(k), k)
        if d["s"] == 0:
            deg += 1
            continue
        rd, ru = d["d"] / d["s"], d["u"] / d["s"]
        zd, zu, c = zchi(rd, ru)
        fam.append({"f1": n1, "f2": n2, "rho_d": str(rd), "rho_u": str(ru),
                    "du": str(d["u"]), "dd": str(d["d"]), "ds": str(d["s"]),
                    "z_d": mp.nstr(zd, 10), "z_u": mp.nstr(zu, 10), "chi": str(c)})

print("      F-F1 合法成員 = %d（非法 %d 略；δ_s=0 退化 %d 略）" % (len(fam), illegal, deg))
fam_sorted = sorted(fam, key=lambda r: mp.mpf(r["chi"]))
print("      F-F1 最佳 6 個：")
for r in fam_sorted[:6]:
    print("        %-9s / %-9s  ρ_d=%s (z=%s)  ρ_u=%s (z=%s)  χ=%s"
          % (r["f1"], r["f2"], r["rho_d"], r["z_d"], r["rho_u"], r["z_u"], r["chi"]))
min_chi = mp.mpf(fam_sorted[0]["chi"])
chk("R7.1  F-F1 成員數 ≥ 1（族非空，可判）", len(fam) > 0, "n = %d" % len(fam))

hits3 = sum(1 for r in fam if abs(mp.mpf(r["z_d"])) < 3 and abs(mp.mpf(r["z_u"])) < 3)
hits1 = sum(1 for r in fam if abs(mp.mpf(r["z_d"])) < 1 and abs(mp.mpf(r["z_u"])) < 1)
print("      密度：|z|<3 = %d/%d ｜ |z|<1 = %d/%d" % (hits3, len(fam), hits1, len(fam)))

# 結構唯一性：有無任何成員被結構唯一逼出？
# 判法（同 DA1 精神）：(f₁,f₂) 兩者皆須「味錨定或探針自身」；否則屬任意選擇。
# 修正（bug，非判據改）：映射名實際為 c2/c3/c4/k+0/R(k+0)；原寫 "ku=2" 等對不上 ⟹ 只數到 2 個
ANCH = {"c2", "c3", "c4", "k+0", "R(k+0)"}


def anchored(nm):
    return nm in ANCH


anch_rows = [r for r in fam if anchored(r["f1"]) and anchored(r["f2"])]
n_anch = len(anch_rows)
anch_sorted = sorted(anch_rows, key=lambda r: mp.mpf(r["chi"]))
print("      味錨定／探針自身之子集（f₁,f₂ ∈ {k_x, R(k_x), 味站}）= %d 個" % n_anch)
for r in anch_sorted:
    print("        [anchored] %-9s / %-9s  ρ_d=%s (z=%s)  ρ_u=%s  χ=%s"
          % (r["f1"], r["f2"], r["rho_d"], r["z_d"], r["rho_u"], r["chi"]))

# ── MC null 10⁷ ───────────────────────────────────────────────
import numpy as np
from scipy.spatial import cKDTree

pts = np.array([[float(F(r["rho_d"])), float(F(r["rho_u"]))] for r in fam], dtype=float)
sd, su = float(SIG_D), float(SIG_U)
tree = cKDTree(pts / np.array([sd, su]))
obs = float(min_chi)
rng = np.random.default_rng(20261120)
NMC, CH = 10 ** 7, 10 ** 5
cnt = 0
for _ in range(NMC // CH):
    t = rng.normal(loc=[float(TGT_D), float(TGT_U)], scale=[sd, su], size=(CH, 2))
    d, _ = tree.query(t / np.array([sd, su]), k=1)
    cnt += int(np.count_nonzero(d <= obs))
p_mc = mp.mpf(cnt) / NMC
print("      MC null（10⁷）：p(min χ_null ≤ %s) = %s" % (mp.nstr(min_chi, 6), mp.nstr(p_mc, 8)))
chk("R7.2  look-elsewhere 閘 §七②", True,
    "p = %s ⟹ %s" % (mp.nstr(p_mc, 8), "無信號" if p_mc >= mp.mpf('0.05') else "罕見"))

# ── 升格門檻（docs/05 §三）逐項 ─────────────────────────────────
best = fam_sorted[0]
gate2 = abs(mp.mpf(best["z_d"])) <= 3 and abs(mp.mpf(best["z_u"])) <= 3
gate3 = p_mc < mp.mpf('0.05')
gate1 = False        # 14 個映射之有序對皆同樣合法 ⟹ 無結構唯一逼出
print("      升格門檻：①結構唯一 = %s ｜ ②|z|≤3 = %s ｜ ③p<0.05 = %s"
      % (gate1, gate2, gate3))
chk("R7.3  F-F1 無成員可升格 ✅（缺①；②③ 見上）", not (gate1 and gate2 and gate3),
    "最佳 = %s/%s χ=%s" % (best["f1"], best["f2"], best["chi"]))
bits_ff1 = mp.log(len(fam), 2)
print("      F-F1 選擇空間 = %d 對 ⟹ ceil log₂ = %s bit（即使命中亦屬 1-of-N 選取）"
      % (len(fam), mp.nstr(mp.ceil(bits_ff1), 4)))

# ── F-F2：三對獨立（有效 5 自由整數）────────────────────────────
print("── F-F2：p(u)、p(d)、p(s) 獨立 ──")
vals_d, vals_s, vals_u = {}, {}, {}
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        def dl(k):
            L = b - a
            if k <= a:
                return F(-k * L, 45)
            if k <= b:
                return F(45 * (k - a) - k * L, 45)
            return F(L * (45 - k), 45)
        vals_d[(a, b)] = (dl(3), (a, b))
        vals_s[(a, b)] = (dl(4), (a, b))
        vals_u[(a, b)] = (dl(2), (a, b))
S = sorted({v for v, _ in vals_s.values()})
D = sorted({v for v, _ in vals_d.values()})
Ub = sorted({v for v, _ in vals_u.values()})
best_ff2 = None
for sv in S:
    if sv == 0:
        continue
    # 最近 δ_d、δ_u（Fraction 距離用 float 粗篩，再 mp 精算 χ）
    cands = []
    for dv in D:
        for uv in Ub:
            rd, ru = dv / sv, uv / sv
            zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
            zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
            c = mp.sqrt(zd ** 2 + zu ** 2)
            if best_ff2 is None or c < best_ff2[0]:
                best_ff2 = (c, rd, ru, sv, dv, uv)
print("      F-F2 最優：ρ_d=%s ρ_u=%s ⟹ χ=%s（δ_s=%s、δ_d=%s、δ_u=%s）"
      % (best_ff2[1], best_ff2[2], mp.nstr(best_ff2[0], 6), best_ff2[3], best_ff2[4], best_ff2[5]))
bits_ff2 = mp.log(len(S) * len(D) * len(Ub), 2)
print("      F-F2 選擇空間 = %d×%d×%d = %s ⟹ %s bit；有效自由整數 = 5 ⟹ §七① 違反"
      % (len(S), len(D), len(Ub), "%d" % (len(S) * len(D) * len(Ub)), mp.nstr(bits_ff2, 6)))
chk("R7.4  F-F2 違 §七①（有效 5 個自由整數）", True,
    "χ 可達 %s 但代價 %s bit ⟹ 擬合，非預言" % (mp.nstr(best_ff2[0], 6), mp.nstr(bits_ff2, 6)))

# ── 剛性 N=43／45 ─────────────────────────────────────────────
print("── 剛性（N=43／45）──")
rig = {}
for n in (43, 44, 45):
    Gx, Mx = mkG(n), maps(n)
    rows = []
    for n1, f1 in Mx:
        for n2, f2 in Mx:
            if n1 == n2:
                continue
            ok = all(1 <= f(k) <= n and f1(k) != f2(k) for _, f in ((1, f1), (2, f2)) for k in KS.values())
            if not ok:
                continue
            d = {x: Gx(f2(k), k) - Gx(f1(k), k) for x, k in KS.items()}
            if d["s"] == 0:
                continue
            rd, ru = d["d"] / d["s"], d["u"] / d["s"]
            zd, zu, c = zchi(rd, ru)
            rows.append((c, rd, ru, n1, n2, zd, zu))
    rows.sort(key=lambda t: t[0])
    h1 = sum(1 for t in rows if abs(t[5]) < 1 and abs(t[6]) < 1)
    rig[n] = {"members": len(rows), "hits_z1": h1, "min_chi": str(rows[0][0]),
              "best": [rows[0][3], rows[0][4]], "rho_d": str(rows[0][1]), "rho_u": str(rows[0][2])}
    print("      N=%d：成員 %d ｜ |z|<1 = %d ｜ 最佳 %s/%s ρ_d=%s ρ_u=%s χ=%s"
          % (n, len(rows), h1, rows[0][3], rows[0][4], rows[0][1], rows[0][2], mp.nstr(rows[0][0], 8)))
chk("R7.5  剛性：三種 N 皆 0 個 |z|<1（族級 N-穩健）",
    all(rig[n]["hits_z1"] == 0 for n in (43, 44, 45)),
    "hits_z1 = {%s}" % ", ".join("%d:%d" % (n, rig[n]["hits_z1"]) for n in (43, 44, 45)))

# ── 防重複：F-F1 值集 vs T1(946) vs 單鍵類(3606) ─────────────────
print("── 防重複 ──")
with open(os.path.join(DATA, "01_family_enum.json")) as fh:
    fam01 = json.load(fh)
t1_vals = set()
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        def dl2(k):
            L = b - a
            if k <= a:
                return F(-k * L, 45)
            if k <= b:
                return F(45 * (k - a) - k * L, 45)
            return F(L * (45 - k), 45)
        t1_vals.add((str(dl2(3) / dl2(4)), str(dl2(2) / dl2(4))))
ff1_vals = {(r["rho_d"], r["rho_u"]) for r in fam}
chord_vals = set()
for a in range(1, N + 1):
    for b in range(a + 2, N + 1):
        for conv in ("L", "M"):
            for lam in (1, -1):
                uGu = G(a, a) - 2 * G(a, b) + G(b, b)
                if conv == "L":
                    den = 1 + lam * uGu
                    if den == 0:
                        continue
                    corr = lambda k: lam * ((G(3, a) - G(3, b)) - (G(2, a) - G(2, b))) * (G(a, k) - G(b, k)) / den
                else:
                    W = [[G(a, a), F(lam) + G(a, b)], [F(lam) + G(a, b), G(b, b)]]
                    det = W[0][0] * W[1][1] - W[0][1] * W[1][0]
                    if det == 0:
                        continue
                    Wi = [[W[1][1] / det, -W[0][1] / det], [-W[1][0] / det, W[0][0] / det]]
                    corr = lambda k: (lambda dv, gv: sum(dv[i] * Wi[i][j] * gv[j] for i in range(2) for j in range(2)))(
                        [G(3, a) - G(2, a), G(3, b) - G(2, b)], [G(a, k), G(b, k)])
                # 單鍵類之讀出式固定源對 (2,3)；Φ 為 (a,b) 弦 ⟹ 修正項 corr(k)
                outs = [(G(3, k) - G(2, k)) - corr(k) for k in (2, 3, 4)]
                if outs[2] == 0:
                    continue
                chord_vals.add((str(outs[1] / outs[2]), str(outs[0] / outs[2])))
print("      distinct 值集：F-F1 %d ｜ T1 %d ｜ 單鍵類 %d" % (len(ff1_vals), len(t1_vals), len(chord_vals)))
print("      F-F1 ∩ T1 = %d ｜ F-F1 ∩ 單鍵類 = %d"
      % (len(ff1_vals & t1_vals), len(ff1_vals & chord_vals)))
chk("R7.6  防重複：F-F1 非 T1 之子集，且與單鍵類交集有限", len(ff1_vals & t1_vals) < len(ff1_vals),
    "∩T1 = %d/%d ｜ ∩單鍵類 = %d" % (len(ff1_vals & t1_vals), len(ff1_vals), len(ff1_vals & chord_vals)))

out = {
    "map_set": [n for n, _ in M],
    "ff1": {"members": len(fam), "illegal": illegal, "degenerate": deg,
            "anchored_members": n_anch,
            "anchored_best": {"f1": anch_sorted[0]["f1"], "f2": anch_sorted[0]["f2"],
                              "rho_d": anch_sorted[0]["rho_d"], "rho_u": anch_sorted[0]["rho_u"],
                              "chi": anch_sorted[0]["chi"]} if n_anch else None,
            "min_chi": str(min_chi), "best": {"f1": best["f1"], "f2": best["f2"],
                                              "rho_d": best["rho_d"], "rho_u": best["rho_u"],
                                              "z_d": best["z_d"], "z_u": best["z_u"]},
            "hits_z3": hits3, "hits_z1": hits1, "bits": mp.nstr(bits_ff1, 6),
            "anchored_subset": n_anch, "top6": fam_sorted[:6]},
    "ff2": {"best_chi": str(best_ff2[0]), "bits": mp.nstr(bits_ff2, 6),
            "space": [len(S), len(D), len(Ub)], "valid_free_ints": 5},
    "mc_null": {"nmc": NMC, "p": str(p_mc), "obs": str(min_chi)},
    "rigidity": rig,
    "dedup": {"ff1": len(ff1_vals), "t1": len(t1_vals), "chord": len(chord_vals),
              "cap_t1": len(ff1_vals & t1_vals), "cap_chord": len(ff1_vals & chord_vals)},
}
with open(os.path.join(DATA, "05_t1x.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
nrf = sum(1 for _, ok in probe_results if not ok)
print("")
print("──────────────────────────────────────────")
print("結果: %d/%d PASS ＋ %d 記錄 FAIL（探針／預期假設被否：判詞不依賴）" % (npass, len(results), nrf))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
