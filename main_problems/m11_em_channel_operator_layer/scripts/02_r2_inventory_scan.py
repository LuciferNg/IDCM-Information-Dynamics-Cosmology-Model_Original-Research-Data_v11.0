#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M11 scripts/02 — R2：零參數 inventory ℐ 之預註冊掃描 ＋ null ＋ 剛性 ＋ 盲測

判準（事前寫死於 docs/03，本腳本只執行）：
  C3  零參數：候選式之常數只准取 inventory 集 INV（IDCM 已註冊整數），不得擬合
  C5  錨分離：in-sample = 5 錨（π 0／N／Σ／Ξ／K）；out-of-sample = 3 條分割通道
  C6  |z|<1 PASS／1–2 TENSION／≥3 KILL
  C7  null：10⁷ 隨機 4-向量；覆蓋率 > 5% ⟹ 無鑑別力
  C8  剛性：N_d 88 → 86 / 90；候選判詞須不變
  C10 誤差傳播：全部用 data/01 之 σ⟨O⟩

禁：不得見結果加候選；不得以任何方式擬合係數；不得改寫 ANCHORS。
"""

import json
import math
import os
import sys
from fractions import Fraction as F

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J1 = os.path.join(HERE, "data", "01_r1_criteria.json")
OUT = os.path.join(HERE, "data", "02_r2_inventory_scan.json")

# ── school-book constants from R1 (single source) ───────────────────────────
D1 = json.load(open(J1, encoding="utf-8"))
X = D1["C2"]["X"]
ROWS = D1["C2"]["rows"]
C9 = D1["C9"]["channels"]

ANCH = ["n-p", "K+-K0", "S--S0", "Xi--Xi0"]      # 4 個非零錨（π 由 annihilation 律自動 0）
PI = "pi+-pi0"
BLIND = ["S0-S+", "S--S+", "D+-D0"]
A_O = {k: ROWS[k]["O"] for k in ANCH}
A_S = {k: ROWS[k]["sigma_O"] for k in ANCH}

# ── 通道結構特徵（final state 量子數；非擬合，來自夸克內容）──────────────────
#  nq  = 夸克數（2/3）｜ns = |奇異數|｜Q = |末態總電荷|｜ann = 價夸克湮滅通道
#  nsub = 被替換夸克數（0 表湮滅）
CH = {
    "n-p":     dict(nq=3, ns=0, Q=0, ann=0, nsub=1),
    "K+-K0":   dict(nq=2, ns=1, Q=1, ann=0, nsub=1),
    "pi+-pi0": dict(nq=2, ns=0, Q=1, ann=1, nsub=0),
    "S--S0":   dict(nq=3, ns=1, Q=1, ann=0, nsub=1),
    "S0-S+":   dict(nq=3, ns=1, Q=0, ann=0, nsub=1),
    "S--S+":   dict(nq=3, ns=1, Q=1, ann=0, nsub=2),
    "Xi--Xi0": dict(nq=3, ns=2, Q=1, ann=0, nsub=1),
    "D+-D0":   dict(nq=2, ns=0, Q=1, ann=0, nsub=1),
}

INV = [1, 2, 3, 4, 5, 7, 9, 12, 18, 29, 41, 44, 45, 47, 55, 88, 89, 90]
RANGE = 2.6          # null 之均勻取樣範圍（錨值全落 [0.0, 2.5]）

FEATS = {
    "nq":   lambda c: c["nq"],
    "ns":   lambda c: c["ns"],
    "ns1":  lambda c: c["ns"] + 1,
    "nqns": lambda c: c["nq"] * c["ns"],
    "nsp":  lambda c: c["nq"] - 1,
    "Q":    lambda c: c["Q"],
    "nsub": lambda c: c["nsub"],
}


def build_family(inv):
    """預註冊候選族：常數只來自 inv；形式固定 8 族。"""
    fam = []
    for kn, fn in FEATS.items():
        for k in inv:
            fam.append((f"G1|{k}/({k}-{kn})", lambda c, f=fn, k=k: k / (k - f(c))))
            fam.append((f"G2|({k}+{kn})/{k}", lambda c, f=fn, k=k: (k + f(c)) / k))
            fam.append((f"G3|1+{kn}/{k}", lambda c, f=fn, k=k: 1.0 + f(c) / k))
            fam.append((f"G6|{kn}/{k}", lambda c, f=fn, k=k: f(c) / k))
        for k in inv:
            fam.append((f"G4|({k}-nq)/({k}-nq-ns)", lambda c, k=k: (k - c["nq"]) / (k - c["nq"] - c["ns"])))
            fam.append((f"G5a|{k}/({k}-nq-ns)", lambda c, k=k: k / (k - c["nq"] - c["ns"])))
            fam.append((f"G5b|{k}/({k}-nq*ns)", lambda c, k=k: k / (k - c["nq"] * c["ns"])))
            fam.append((f"G7|{k}/({k}-nsub*ns)", lambda c, k=k: k / (k - c["nsub"] * c["ns"])))
            fam.append((f"G8|({k}+nsub*ns)/{k}", lambda c, k=k: (k + c["nsub"] * c["ns"]) / k))
    # G9｜γ 譜階梯比（純 inventory：K_chain 譜比，不含通道→模映射）
    G = gamma_ladder(44)
    fam.append(("G9|g[nq+ns]/g[nq]", lambda c: G[c["nq"] + c["ns"]] / G[c["nq"]]))
    fam.append(("G9|g[nq+ns+1]/g[nq+1]", lambda c: G[c["nq"] + c["ns"] + 1] / G[c["nq"] + 1]))
    fam.append(("G9|g[nsub*ns+nq]/g[nq]", lambda c: G[c["nsub"] * c["ns"] + c["nq"]] / G[c["nq"]]))
    fam.append(("G9|g[nq+ns]/g[nq+ns-1]", lambda c: G[c["nq"] + c["ns"]] / G[c["nq"] + c["ns"] - 1]))
    return fam


def gamma_ladder(nmax):
    """γ_n = 2 − 2T_n(c)，c = 1 − γ₁/2（GAMMA1 由 R4 譜代數，M4 docs/10）。無自由參數。"""
    G1 = 0.005095770779492918
    c = 1.0 - G1 / 2.0
    th = math.acos(c)
    return [None] + [2.0 - 2.0 * math.cos(n * th) for n in range(1, nmax + 1)]


def val(fn, chname):
    c = CH[chname]
    try:
        v = fn(c)
    except ZeroDivisionError:
        return None
    if not np.isfinite(v):
        return None
    return v * (1 - c["ann"])          # 湮滅律：價夸克湮滅 ⟹ ⟨O⟩ = 0（Dashen）


def main():
    log = []
    def emit(s=""):
        log.append(s); print(s)

    emit("M11 R2 — 零參數 inventory 掃描（in-sample 4 錨｜blind 3 通道）")
    emit("=" * 100)
    emit(f"錨值來源 data/01（實驗 total − lattice QED 法）；X = {X:.2f}")
    for k in ANCH:
        emit(f"  {k:9s} ⟨O⟩ = {A_O[k]:.6f} ± {A_S[k]:.6f}")
    emit(f"  {PI:9s} ⟨O⟩ = 0（Dashen；由 annihilation 律自動滿足，不具檢定力）")
    emit(f"inventory INV = {INV}")
    emit(f"特徵 = {list(FEATS)}")

    fam = build_family(INV)
    emit(f"候選總數 = {len(fam)}（9 形式族 G1–G9；常數只取 INV；零自由係數）")

    # ── 評估全部候選 ──
    recs = []
    for name, fn in fam:
        vs = {ch: val(fn, ch) for ch in CH}
        if any(v is None for v in vs.values()):
            continue
        z = {k: (vs[k] - A_O[k]) / A_S[k] for k in ANCH}
        maxz = max(abs(v) for v in z.values())
        recs.append(dict(name=name, vals=vs, z=z, maxz=maxz))
    emit(f"有效候選（分母非零） = {len(recs)}")

    # ── 判詞分級（C6）──
    for r in recs:
        m = r["maxz"]
        r["grade"] = "PASS" if m < 1 else ("TENSION" if m < 2 else ("STRAIN" if m < 3 else "KILL"))
    cnt = {g: sum(1 for r in recs if r["grade"] == g) for g in ("PASS", "TENSION", "STRAIN", "KILL")}
    emit("")
    emit(f"C6 分級：{cnt}")
    passes = [r for r in recs if r["grade"] == "PASS"]
    passes.sort(key=lambda r: r["maxz"])
    emit(f"PASS 候選 = {len(passes)}")
    for r in passes[:12]:
        emit(f"  {r['maxz']:6.3f}  {r['name']:28s}  "
             + " ".join(f"{k}={r['vals'][k]:.3f}" for k in ANCH))

    # ── 最接近者（即使 KILL 也列，供結構診斷）──
    best = sorted(recs, key=lambda r: r["maxz"])[:10]
    emit("")
    emit("最接近之 10 個候選（KILL 亦列）：")
    for r in best:
        killer = max(r["z"], key=lambda k: abs(r["z"][k]))
        emit(f"  {r['maxz']:7.2f}  {r['name']:26s} " +
             " ".join(f"{k}={r['vals'][k]:6.3f}" for k in ANCH) +
             f"   | 最惡錨={killer} z={r['z'][killer]:+.2f}")
    kill_anchor = {}
    for r in recs:
        kk = max(r["z"], key=lambda k: abs(r["z"][k]))
        kill_anchor[kk] = kill_anchor.get(kk, 0) + 1
    emit(f"  失敗錨分佈（首惡）：{kill_anchor}")

    # ── C7 null（10⁷ 隨機 4-向量；尺度化 Chebyshev KD-tree，精確判箱）──
    uniq = sorted({tuple(round(r["vals"][k], 9) for k in ANCH) for r in recs})
    U = np.array(uniq, float)
    S = np.array([A_S[k] for k in ANCH], float)
    tree = cKDTree(U / S)            # 每維除以 σ ⟹ 箱變單位立方
    M = 10 ** 7
    rng = np.random.default_rng(20261028)
    hits = 0
    CHUNK = 250_000
    done = 0
    while done < M:
        n = min(CHUNK, M - done)
        T = rng.uniform(0.0, RANGE, size=(n, 4)) / S      # 尺度化
        d, _ = tree.query(T, k=1, p=np.inf)
        hits += int(np.count_nonzero(d <= 1.0))
        done += n
    cov = hits / M
    analytic = len(uniq) * float(np.prod(2 * S / RANGE))
    emit("")
    emit(f"C7 null：候選族（去重 {len(uniq)} 個不同 4-向量）× {M:.0e} 隨機目標")
    emit(f"  隨機命中 = {hits}  ⟹ 覆蓋率 = {cov:.3e}（判準 > 5% 即無鑑別力）")
    emit(f"  獨立箱體積和 = {analytic:.3e}（**上界，未計重疊**；候選在部分維度高度重疊 ⟹ 大於 MC 真值屬正常）")
    emit(f"  ⟹ 檢定力：α = {cov:.3%}（隨機目標被族命中之機率）；實際 PASS = {len(passes)} 個")

    # ── C8 剛性（N_d 88 → 86/90）──
    emit("")
    rigid = {}
    for alt in (86, 90):
        inv2 = [alt if k == 88 else k for k in INV]
        fam2 = build_family(inv2)
        npass = 0
        for name, fn in fam2:
            vs = {ch: val(fn, ch) for ch in CH}
            if any(v is None for v in vs.values()):
                continue
            m = max(abs((vs[k] - A_O[k]) / A_S[k]) for k in ANCH)
            if m < 1:
                npass += 1
        rigid[alt] = npass
        emit(f"C8 N_d→{alt}：PASS 候選 = {npass}（原 88：{len(passes)}）")

    # ── C5 盲測：PASS 候選之 3 條 out-of-sample 分割 ──
    emit("")
    emit("C5 盲測（out-of-sample）：Σ⁰−Σ⁺、Σ⁻−Σ⁺、D⁺−D⁰ 之 (⟨O⟩, em) 分割")
    k0 = D1["C9"]["k0_MeV"]
    blind_out = []
    for r in passes[:20]:
        line = []
        ok = True
        for ch in BLIND:
            c = C9[ch]
            O = r["vals"][ch]
            coeff = ROWS[ch]["coeff"]
            em = ROWS[ch]["dm"] - coeff * X * O
            t = em / k0
            lo, hi = min(c["box_min"], c["box_max"]), max(c["box_min"], c["box_max"])
            feas = lo - 1e-9 <= t <= hi + 1e-9
            ok = ok and feas
            line.append(f"{ch}: ⟨O⟩={O:.3f} em={em:+.3f} {'OK' if feas else 'OUT'}")
        blind_out.append(dict(name=r["name"], maxz=r["maxz"], ok=bool(ok), detail=line))
        emit(f"  [{r['maxz']:.2f}] {'PASS' if ok else 'FAIL'}  {r['name']}")
        for s in line:
            emit(f"        {s}")

    # ── JSON ──
    res = dict(
        anchors={k: dict(O=A_O[k], sigma=A_S[k]) for k in ANCH},
        pi_rule="annihilation ⟹ ⟨O⟩=0（Dashen），對所有候選自動成立 ⟹ 不具檢定力",
        inventory_INV=INV, features=list(FEATS), n_candidates=len(recs),
        grades=cnt,
        pass_candidates=[dict(name=r["name"], maxz=r["maxz"], vals=r["vals"],
                              z=r["z"]) for r in passes],
        null=dict(M=M, uniq_vectors=len(uniq), hits=hits, coverage=cov, analytic=analytic),
        rigidity={str(k): v for k, v in rigid.items()},
        blind=blind_out,
    )
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    emit("")
    emit(f"[JSON] {OUT}")
    emit(f"[LINES] {len(log)}")


if __name__ == "__main__":
    main()