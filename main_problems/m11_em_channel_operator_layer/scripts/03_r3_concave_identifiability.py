#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M11 scripts/03 — R3：凹形族之鑑別力崩塌（ℐ₃，**後驗設計，已標記**）

依據 docs/03 §五.4：R2 之 9 族全為線性／凸；數據在 ns 上係凹（+0.783 → +0.263）。
本輪把凹形補上，問一個更關鍵嘅問題：**補得上之後仲有冇鑑別力？**

判準（事前寫死）：
  R3-1 後驗標記：本族由 R2 結果反推 ⟹ 全族標 posthoc，任何命中不得升格為結論
  R3-2 鑑別力門檻：α（隨機目標被族命中率）必須 ≤ 1%；α > 5% ⟹ 族無鑑別力（即「任何凹形擬合皆不可辨識」）
  R3-3 盲測：若 PASS 集合非空，其 3 條盲測預言之散佈須窄（否則無預言力）
  R3-4 剛性：N_d 88 → 86/90

形式（零自由係數；常數只取 INV）：
  G10｜⟨O⟩ = A − B/(C + f)，f ∈ {ns, nq, nsub}, A,B,C ∈ INV
  G11｜⟨O⟩ = A/(C + f)          （反向檢查：單調遞減族是否同樣無力）
全部乘 (1−ann)（湮滅律，同 R2）。
"""

import json
import math
import os

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D1 = json.load(open(os.path.join(HERE, "data", "01_r1_criteria.json"), encoding="utf-8"))
OUT = os.path.join(HERE, "data", "03_r3_concave.json")

ROWS = D1["C2"]["rows"]
ANCH = ["n-p", "K+-K0", "S--S0", "Xi--Xi0"]
BLIND = ["S0-S+", "S--S+", "D+-D0"]
A_O = {k: ROWS[k]["O"] for k in ANCH}
A_S = {k: ROWS[k]["sigma_O"] for k in ANCH}

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
RANGE = 2.6


def family(inv):
    fam = []
    for fn_name in ("ns", "nq", "nsub"):
        for A in inv:
            for B in inv:
                for C in inv:
                    fam.append((f"G10|{A}-{B}/({C}+{fn_name})",
                                lambda c, A=A, B=B, C=C, f=fn_name: A - B / (C + c[f])))
        for A in inv:                      # G11 不含 B ⟹ 獨立迴圈（避免 18× 重複）
            for C in inv:
                fam.append((f"G11|{A}/({C}+{fn_name})",
                            lambda c, A=A, C=C, f=fn_name: A / (C + c[f])))
    return fam


def ev(fn, ch):
    c = CH[ch]
    try:
        v = fn(c)
    except ZeroDivisionError:
        return None
    if not np.isfinite(v):
        return None
    return v * (1 - c["ann"])


def main():
    log = []
    def emit(s=""):
        log.append(s); print(s)

    emit("M11 R3 — 凹形族 ℐ₃：補上形狀缺口之後，鑑別力剩幾多？（後驗設計）")
    emit("=" * 100)
    emit("R3-1｜本族由 R2 形狀診斷反推 ⟹ 全族 posthoc；任何命中不得升格")
    fam = family(INV)
    emit(f"候選 = {len(fam)}（G10＋G11 × A,B,C ∈ INV × f ∈ {{ns,nq,nsub}}）")

    recs = []
    for name, fn in fam:
        vs = {ch: ev(fn, ch) for ch in CH}
        if any(v is None for v in vs.values()):
            continue
        z = {k: (vs[k] - A_O[k]) / A_S[k] for k in ANCH}
        recs.append(dict(name=name, vals=vs, z=z, maxz=max(abs(v) for v in z.values())))
    emit(f"有效候選 = {len(recs)}")

    for r in recs:
        m = r["maxz"]
        r["grade"] = "PASS" if m < 1 else ("TENSION" if m < 2 else ("STRAIN" if m < 3 else "KILL"))
    cnt = {g: sum(1 for r in recs if r["grade"] == g) for g in ("PASS", "TENSION", "STRAIN", "KILL")}
    passes = sorted([r for r in recs if r["grade"] == "PASS"], key=lambda r: r["maxz"])
    emit(f"C6 分級：{cnt}｜PASS = {len(passes)}")
    for r in passes[:8]:
        emit(f"  {r['maxz']:6.3f}  {r['name']:26s} " +
             " ".join(f"{k}={r['vals'][k]:6.3f}" for k in ANCH))

    best = sorted(recs, key=lambda r: r["maxz"])[:10]
    emit("")
    emit("最接近之 10 個（KILL 亦列）：")
    for r in best:
        kk = max(r["z"], key=lambda k: abs(r["z"][k]))
        emit(f"  {r['maxz']:7.2f}  {r['name']:26s} " +
             " ".join(f"{k}={r['vals'][k]:6.3f}" for k in ANCH) +
             f"  | 首惡={kk} z={r['z'][kk]:+.2f}")
    kill_a = {}
    for r in recs:
        kk = max(r["z"], key=lambda k: abs(r["z"][k]))
        kill_a[kk] = kill_a.get(kk, 0) + 1
    emit(f"  失敗錨分佈（首惡）：{kill_a}")

    # α（10⁷）
    uniq = sorted({tuple(round(r["vals"][k], 9) for k in ANCH) for r in recs})
    U = np.array(uniq, float)
    S = np.array([A_S[k] for k in ANCH], float)
    tree = cKDTree(U / S)
    M = 10 ** 7
    rng = np.random.default_rng(20261029)
    hits, done, CHUNK = 0, 0, 250_000
    while done < M:
        n = min(CHUNK, M - done)
        T = rng.uniform(0.0, RANGE, size=(n, 4)) / S
        d, _ = tree.query(T, k=1, p=np.inf)
        hits += int(np.count_nonzero(d <= 1.0))
        done += n
    cov = hits / M
    emit("")
    emit(f"R3-2｜α = {cov:.4f}（{hits}/{M}）｜判準：α ≤ 1% 才可用；>5% ⟹ 無鑑別力")
    emit(f"  去重 4-向量 = {len(uniq)}")
    verdict = ("無鑑別力（α > 5%）" if cov > 0.05 else
               ("勉強可用（1% < α ≤ 5%）" if cov > 0.01 else "可用（α ≤ 1%）"))
    emit(f"  ⟹ {verdict}")

    # 盲測散佈
    emit("")
    emit("R3-3｜PASS 集合之盲測預言散佈（若闊 ⟹ 無預言力）")
    blind_stats = {}
    if passes:
        for ch in BLIND:
            vals = sorted(r["vals"][ch] for r in passes)
            lo, hi = vals[0], vals[-1]
            blind_stats[ch] = dict(lo=lo, hi=hi, n=len(vals),
                                   med=vals[len(vals) // 2])
            emit(f"  {ch:9s} PASS 預言 ⟨O⟩ ∈ [{lo:.3f}, {hi:.3f}]（{len(vals)} 個，中位 {vals[len(vals)//2]:.3f}）")
        # coherent check：PASS 集合中有幾多個同時令 3 條盲測都箱可行
        k0 = D1["C9"]["k0_MeV"]
        C9 = D1["C9"]["channels"]
        coherent = 0
        for r in passes:
            ok = True
            for ch in BLIND:
                c = C9[ch]
                em = ROWS[ch]["dm"] - ROWS[ch]["coeff"] * D1["C2"]["X"] * r["vals"][ch]
                t = em / k0
                lo, hi = min(c["box_min"], c["box_max"]), max(c["box_min"], c["box_max"])
                if not (lo - 1e-9 <= t <= hi + 1e-9):
                    ok = False
                    break
            coherent += int(ok)
        emit(f"  同時通過 3 條盲測箱可行性者 = {coherent}/{len(passes)}")
    else:
        emit("  PASS 集合為空 ⟹ 無盲測可言")
        coherent = 0

    # 剛性
    emit("")
    rigid = {}
    for alt in (86, 90):
        inv2 = [alt if k == 88 else k for k in INV]
        npass = 0
        for name, fn in family(inv2):
            vs = {ch: ev(fn, ch) for ch in CH}
            if any(v is None for v in vs.values()):
                continue
            m = max(abs((vs[k] - A_O[k]) / A_S[k]) for k in ANCH)
            if m < 1:
                npass += 1
        rigid[alt] = npass
        emit(f"R3-4｜N_d→{alt}：PASS = {npass}（原 88：{len(passes)}）")

    res = dict(family_n=len(fam), valid=len(recs), grades=cnt,
               n_pass=len(passes), alpha=cov, uniq=len(uniq), verdict=verdict,
               pass_examples=[dict(name=r["name"], maxz=r["maxz"], vals=r["vals"])
                              for r in passes[:20]],
               best10=[dict(name=r["name"], maxz=r["maxz"], vals=r["vals"], z=r["z"])
                       for r in best],
               blind_spread=blind_stats, coherent_blind=coherent,
               rigidity={str(k): v for k, v in rigid.items()},
               posthoc=True)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    emit("")
    emit(f"[JSON] {OUT}")
    emit(f"[LINES] {len(log)}")


if __name__ == "__main__":
    main()