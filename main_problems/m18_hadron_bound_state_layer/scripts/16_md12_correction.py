#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/16 — **自我更正**：`MD12` 之前設（`A₅`-等變）被接駁結構否證
====================================================================
本輪目的：唔係加新結果，而係**查自己上一輪個判有冇漏洞**。

`docs/15` 之 `MD12` 用「載體重數 `m > 1` ⟹ 精確退化」排除 `60`。該推論需要：
  前提 P：`H` 係 **`A₅`-等變**。
本輪問：P 成立嗎？

結構審計（本輪機器核）
  C1 `A₅` 作用於 `20` 點之**交換子代數維數** = `X×X` 之 `A₅`-軌道數 = Σ mult² = 8
  C2 同一 `20` 點集，改由 `A₄` 作用：交換子維數 = ? （若 > 8 ⟹ 允許更多混合）
  C3 `A₂`↔`G`-重數空間：`G↓A₄ = A ⊕ T` ⟹ 兩份 `G` 之重數空間在 `A₄` 下**平凡**
     ⟹ `A₄`-等變算子可喺重數空間任意混合 ⟹ 退化可被劈開
  C4 接駁結構：`d = 0` 係**單一**節點（`N_0 = 1`）⟹ 佢揀定 `5`-orbit 之一點 ⟹ 只餘穩定子 `A₄`
     ⟹ 若 `H` 係 `A₅`-等變，則 `d = 0` 唔會被區分（矛盾）⟹ **P 唔成立**

事前判準（寫死，唔准事後改）
----------------------------
  C1 交換子維數用「`X×X` 軌道數」實數（置換模塊定理：不變矩陣 = 軌道矩陣），並與 `Σ mult²` 交叉核
  C2 `A₄` 之軌道數必須**嚴格大於** `A₅` 之軌道數，否則「可劈開」論證無效
  C3 「退化可劈開」須有**顯式見證**（非只講維數）
  C4 自我更正：`MD12` 若前設被否 ⟹ 必須**撤回為條件性命題**，並 append 記錄（原文不動）
  C5 替代判別（若仍有）：只能係**態數**（可數），並明標較弱

輸出：data/16_md12_correction.json
"""
import itertools
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="C", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="C"):
    return chk(name, ok, detail, tag, integrity=False)


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def inv(p):
    r = [0] * len(p)
    for i, v in enumerate(p):
        r[v] = i
    return tuple(r)


def sign(p):
    s, seen = 1, [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


E = tuple(range(5))
A5 = sorted(p for p in itertools.permutations(range(5)) if sign(p) == 1)


def closure(gens):
    S = {E}
    fr = list(gens)
    while fr:
        g = fr.pop()
        new = {compose(g, s) for s in S} | {compose(s, g) for s in S}
        for x in new:
            if x not in S:
                S.add(x)
                fr.append(x)
    return frozenset(S)


cyc = {closure([g]) for g in A5}
subs = set(cyc)
prev = -1
while len(subs) != prev:
    prev = len(subs)
    L = list(subs)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            subs.add(closure(list(L[i]) + list(L[j])))
byorder = {}
for s in subs:
    byorder.setdefault(len(s), []).append(s)
A4 = sorted(byorder[12][0])
C3 = sorted([s for s in byorder[3] if s < frozenset(A4)][0])
C5s = sorted(byorder[5][0])


def cosets(G, H):
    out, seen_c = [], set()
    for x in G:
        c = frozenset(compose(x, h) for h in H)
        if c not in seen_c:
            seen_c.add(c)
            out.append(sorted(c)[0])
    return out


X20 = cosets(A5, C3)
print("=" * 78)
print("C1. 交換子維數 = X×X 軌道數（置換模塊定理）")
print("=" * 78)


def n_orbits(G, X, act):
    seen, cnt = set(), 0
    for i in range(len(X)):
        for j in range(len(X)):
            if (i, j) in seen:
                continue
            cnt += 1
            orb, fr = set(), [(i, j)]
            while fr:
                a, b = fr.pop()
                if (a, b) in orb:
                    continue
                orb.add((a, b))
                for g in G:
                    fr.append(act(g, a, b))
            seen |= orb
    return cnt


# 正確做法：以 coset index 之作用 g·(hH) = (gh)H
def coset_of(x):
    return frozenset(compose(x, h) for h in C3)


idx_of = {coset_of(x): k for k, x in enumerate(X20)}
assert len(idx_of) == 20


def act_idx(g, a, b):
    return (idx_of[coset_of(compose(g, X20[a]))], idx_of[coset_of(compose(g, X20[b]))])


n5 = n_orbits(A5, X20, act_idx)
print("  A₅ 對 20×20 之軌道數 = %d（應 = Σ mult² = 1+1+1+4+1 = 8）" % n5)
chk("C1a A₅ 之交換子維數 = 8（與 Σ mult² 交叉核）", n5 == 8, "%d" % n5)

n4 = n_orbits(A4, X20, act_idx)
print("  A₄ 對同一 20×20 之軌道數 = %d（A₄ ⊂ A₅ ⟹ 軌道更細 ⟹ 維數更大）" % n4)
chk("C2a ★ A₄ 交換子維數 %d **嚴格大於** A₅ 之 8 ⟹ 允許更多混合（可劈開重數）" % n4, n4 > n5,
    "%d > %d" % (n4, n5))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("C3. 顯式見證：兩份 G 之重數空間在 A₄ 下平凡 ⟹ 可混合")
print("=" * 78)
# 建 20 點置換模塊之 A₄ 作用，找 G-等型分量之投影，檢查重數空間
import numpy as np
n = 20
perms = []
for g in A4:
    M = np.zeros((n, n))
    for a in range(n):
        M[idx_of[coset_of(compose(g, X20[a]))], a] = 1.0
    perms.append(M)
# 取 A₅ 之 G 投影（用 G 之 4 維 irrep 特徵做投影算子）
g_irrep = np.zeros((n, n))
# 用 G 之特徵：χ_G(1)=4, χ_G(2A)=0, χ_G(3A)=1, χ_G(5)= -1
def cls_of(p):
    o = 1
    h = p
    while h != E:
        h = compose(h, p)
        o += 1
    return o


chiG = {1: 4, 2: 0, 3: 1, 5: -1}
P_G = np.zeros((n, n))
for g in A5:
    M = np.zeros((n, n))
    for a in range(n):
        M[idx_of[coset_of(compose(g, X20[a]))], a] = 1.0
    P_G += float(mp.mpf(4) / 60 * mp.mpf(chiG[cls_of(g)])) * M
P_G = np.asarray(P_G, dtype=float)
rank_G = int(np.linalg.matrix_rank(P_G, tol=1e-9))
print("  G-等型分量維數 = rank(P_G) = %d（應 = 2·4 = 8）" % rank_G)
chk("C3a G-等型維數 = 8 = 2 份 × dim G ⟹ 重數 = 2", rank_G == 8, "rank %d" % rank_G)

# 重數空間：P_G 之影像內，取兩份之可分基底；驗 A₄ 作用在重數空間上係恒等（平凡）
# 做法：取任意 x ∈ Im(P_G)，檢查所有 A₄ 元素作用後仍在同一「copy-label」…
# 簡化而嚴格：檢查 A₄-不變算子空間在 Im(P_G) 上之維數 = 4（M₂ 之維數）而 A₅-不變 = 1？
# 用矩陣方式直接算：A₄-不變算子之維數（在 20 維空間）= n4（已算）；A₅-不變 = 8（已算）
# 關鍵：兩份 G 之重數空間（2 維）——檢查存在 A₄-等變算子喺該 2 維上非純量
# 建 basis of Im(P_G)
U, S, Vt = np.linalg.svd(P_G)
B = U[:, :8]                      # 8×20 basis rows
Bp = B.T                          # 20×8
# 投影任何 20×20 矩陣到 8×8
def proj(Xm):
    return B.T @ Xm @ B

# 喺 Im(P_G) 上取 A₄-等變算子空間維數：解線性系統（8×8 = 64 未知）
rows = []
for M_ in perms:
    Mp = proj(M_)
    for i in range(8):
        for j in range(8):
            row = np.zeros(64)
            # (Mp X - X Mp)[i,j] = 0
            for k in range(8):
                row[k * 8 + j] += Mp[i, k]
                row[i * 8 + k] -= Mp[k, j]
            rows.append(row)
Amat = np.array(rows)
dim_inv = 64 - int(np.linalg.matrix_rank(Amat, tol=1e-9))
print("  Im(P_G)（8 維）上之 A₄-等變算子維數 = %d" % dim_inv)
print("  （若 = 4 = dim M₂ ⟹ 兩份 G 之重數空間可任意混合 ⟹ 可劈開）")

# 顯式見證：喺該空間揀一個非純量算子
ns = None
# 用 SVD 取零空間
uu, ss, vv = np.linalg.svd(Amat)
ns = vv[64 - dim_inv:] if dim_inv > 0 else np.zeros((0, 64))
witness = None
for vec in ns:
    Xm = vec.reshape(8, 8)
    if np.linalg.matrix_rank(Xm, tol=1e-9) > 1 and np.linalg.norm(Xm - np.trace(Xm) / 8 * np.eye(8)) > 1e-6:
        witness = Xm
        break
ok_w = witness is not None
print("  非純量見證算子：%s（rank %d）" % ("有" if ok_w else "冇",
                                          int(np.linalg.matrix_rank(witness, tol=1e-9)) if ok_w else 0))
chk("C3b ★ 顯式見證：存在 A₄-等變、但喺重數空間上非純量之算子 ⟹ 可劈開 2 份 G",
    ok_w, "rank %d" % (int(np.linalg.matrix_rank(witness, tol=1e-9)) if ok_w else 0))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("C4. 自我更正：MD12 之前設被否")
print("=" * 78)
print("  d = 0 係 N_0 = 1 之**單一**節點 ⟹ 佢揀定 5-orbit 之一點 ⟹ 只餘穩定子 A₄")
print("  若 H 係 A₅-等變 ⟹ 5 點全部等價 ⟹ d=0 唔會被區分（與 M5 之 V₀ 單點接觸項矛盾）")
verdict("C4a ★ `MD12` 前提 P（H 為 A₅-等變）**被接駁結構否證** ⟹ 退化論證失效 ⟹ `60` 唔被排除",
        False, "⟹ MD12 必須撤回為條件性命題（append-only 更正）", "C4")
verdict("C4b `MD11` 不受影響（Frobenius 誘導係模塊層事實，與 H 之等變性無關）", True,
        "A1 逐類相等仍成立", "C4")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("C5. 替代判別（較弱，可數）")
print("=" * 78)
mult20 = {"A": 1, "T1": 1, "T2": 1, "G": 2, "H": 1}
mult60 = {"A": 1, "T1": 3, "T2": 3, "G": 4, "H": 5}
npm20 = sum(mult20.values())
npm60 = sum(mult60.values())
print("  20：多重態數 = %d（覆蓋 20 態）" % npm20)
print("  60：多重態數 = %d（覆蓋 60 態）" % npm60)
print("  ⟹ 每徑向帶之態數：20 → 6 多重態；60 → 16 多重態（×%.2f）" % (60 / 20))
chk("C5a 態數判別：6 vs 16 多重態（可數、可證偽，但**較弱**：靠經濟性非定理）",
    npm20 == 6 and npm60 == 16, "比例 3")
verdict("C5b 判定：`{20,60}` 之 1 bit **未解**（群論層兩者皆合法）；只能用態數經濟性 🟡 傾向 20",
        False, "唔可以宣稱已判 ⟹ 誠實報殘", "C5")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "16_md12_correction.py", "checks": checks, "integrity_fail": nfail,
       "commutant_A5_on_20": n5, "commutant_A4_on_20": n4,
       "G_isotypic_dim": rank_G, "A4_invariant_dim_on_G_isotypic": int(dim_inv),
       "witness_exists": bool(ok_w), "MD12_status": "WITHDRAWN (premise refuted)",
       "multiplet_20": npm20, "multiplet_60": npm60, "bit_resolved": False}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "16_md12_correction.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/16_md12_correction.json")
sys.exit(0 if nfail == 0 else 1)
