#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M11 scripts/04 — R4：必要性定理 ＋ 錨飽和度 ＋ 盲測箱檢定（I2 最小實現之前置分析）

問題：`⟨O⟩_ch` 要由結構推出，最少需要哪些接口（I1 電荷讀數／I2 iso-breaking 頂點／I3 通道狀態）？
做法：把「接口之有無」翻譯成對 `⟨O⟩` 之**可加性／狀態依賴**限制，用 4 個錨（in-sample）
      逐個檢定，再用 3 條分割通道之**箱約束**（out-of-sample）作第二道閘。

判準（事前寫死）：
  R4-1 類別模型零擬合：每類一共同值 = 加權平均（唯一 MLE）
  R4-2 錨層：χ²/dof ⟹ z ≥ 3 KILL；2–3 TENSION；<2 SURVIVE
  R4-3 盲測層：由模型推 3 條通道之 (⟨O⟩, em)，以 data/01 C9 之箱（t = em/k₀）判可行；
        任一通道不可行 ⟹ 該模型被盲測否證（OUT）
  R4-4 飽和度：類別數 vs 錨數；0 dof ⟹ 錨耗盡，該模型只係「未被否證」
禁：不得為救任何模型增加類別以外之參數；不得事後改判準。
"""

import json
import os

import numpy as np
from scipy.stats import chi2 as chi2dist, norm

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D1 = json.load(open(os.path.join(HERE, "data", "01_r1_criteria.json"), encoding="utf-8"))
OUT = os.path.join(HERE, "data", "04_r4_necessity.json")

ROWS, C9 = D1["C2"]["rows"], D1["C9"]["channels"]
X = D1["C2"]["X"]
K0 = D1["C9"]["k0_MeV"]

A = [
    ("n-p",     3, 0, ROWS["n-p"]["O"],     ROWS["n-p"]["sigma_O"]),
    ("S--S0",   3, 1, ROWS["S--S0"]["O"],   ROWS["S--S0"]["sigma_O"]),
    ("Xi--Xi0", 3, 2, ROWS["Xi--Xi0"]["O"], ROWS["Xi--Xi0"]["sigma_O"]),
    ("K+-K0",   2, 1, ROWS["K+-K0"]["O"],   ROWS["K+-K0"]["sigma_O"]),
]
BLIND = ["S0-S+", "S--S+", "D+-D0"]


def wmean(pts):
    w = np.array([1.0 / p[1] ** 2 for p in pts])
    v = np.array([p[0] for p in pts])
    return float((w * v).sum() / w.sum()), float(1.0 / np.sqrt(w.sum()))


def build(keyfn, label, note):
    groups = {}
    for n, nq, ns, o, s in A:
        groups.setdefault(keyfn(nq, ns), []).append((o, s))
    chi2v, dof, params = 0.0, 0, {}
    for k, pts in groups.items():
        m, sm = wmean(pts)
        params[k] = dict(value=m, sigma=sm, n=len(pts))
        for o, s in pts:
            chi2v += ((o - m) / s) ** 2
        dof += len(pts) - 1
    p = float(chi2dist.sf(chi2v, dof)) if dof > 0 else 0.0
    z = float(norm.isf(max(p, 1e-300) / 2)) if dof > 0 else float("inf")
    return dict(label=label, note=note, chi2=chi2v, dof=dof, p=p, z=z,
                params=params, n_classes=len(groups))


def main():
    log = []
    def emit(s=""):
        log.append(s); print(s)

    emit("M11 R4 — 必要性定理：⟨O⟩_ch 最少需要哪些接口？（零擬合類別模型 ＋ 盲測箱）")
    emit("=" * 104)
    emit("錨（in-sample）：")
    for n, nq, ns, o, s in A:
        emit(f"  {n:9s} nq={nq} ns={ns}  ⟨O⟩ = {o:.6f} ± {s:.6f}")
    emit(f"盲測（out-of-sample）：{BLIND}；k₀ = {K0:.4f} MeV；X = {X}")

    MODELS = [
        (lambda nq, ns: "all",  "A valence-additive（無 I3）",
         "⟨O⟩ 只依被替換夸克 ⟹ 與通道無關"),
        (lambda nq, ns: f"nq{nq}", "B nq-only", "只記介子／重子（1 bit 狀態）"),
        (lambda nq, ns: f"ns{ns}", "C ns-only", "只記奇異數"),
        (lambda nq, ns: f"({nq},{ns})", "D (nq, ns)", "聯合類別（最小可存活者）"),
    ]

    out = []
    emit("")
    emit("R4-1／R4-2｜錨層（每類一共同值；z ≥ 3 KILL）")
    for keyfn, label, note in MODELS:
        m = build(keyfn, label, note)
        gr = "KILL" if (m["dof"] > 0 and m["z"] >= 3) else (
             "TENSION" if (m["dof"] > 0 and m["z"] >= 2) else (
             "SURVIVE" if m["dof"] > 0 else "SATURATED"))
        m["grade_anchor"] = gr
        emit(f"  {label:34s} k={m['n_classes']}  χ²={m['chi2']:7.2f}  dof={m['dof']}  "
             f"p={m['p']:.2e}  z={m['z']:6.2f}  ⟹ {gr}")
        m["verdicts"] = {}
        for k, v in m["params"].items():
            emit(f"      {k:8s} = {v['value']:.6f} ± {v['sigma']:.6f} (n={v['n']})")

        # R4-3 盲測箱
        emit(f"      R4-3 盲測：")
        nfail = 0
        for ch in BLIND:
            nq_ch = 3 if ch in ("S0-S+", "S--S+") else 2
            ns_ch = 1 if ch in ("S0-S+", "S--S+") else 0
            key = keyfn(nq_ch, ns_ch)
            coeff, dm = ROWS[ch]["coeff"], ROWS[ch]["dm"]
            if key not in m["params"]:
                emit(f"        {ch:8s} 類別 {key} 無錨 ⟹ 無預言（真盲測位）")
                m["verdicts"][ch] = dict(category=key, O=None, em=None, t=None, feasible=None)
                continue
            O = m["params"][key]["value"]
            em = dm - coeff * X * O
            t = em / K0
            c = C9[ch]
            lo, hi = min(c["box_min"], c["box_max"]), max(c["box_min"], c["box_max"])
            feas = bool(lo - 1e-9 <= t <= hi + 1e-9)
            nfail += int(not feas)
            emit(f"        {ch:8s} 類別 {key} ⟹ ⟨O⟩={O:.4f}（強制） em={em:+8.4f} MeV  "
                 f"t={t:+.4f}  箱[{lo:+.4f},{hi:+.4f}] ⟹ {'可行' if feas else '❌ 不可行'}")
            m["verdicts"][ch] = dict(category=key, O=O, em=em, t=t, box=[lo, hi], feasible=feas)
        m["blind_fails"] = nfail
        m["grade"] = "OUT（盲測否證）" if nfail else (
            "KILL（錨）" if gr == "KILL" else ("SATURATED" if gr == "SATURATED" else "SURVIVE"))
        emit(f"      ⟹ 綜合：{m['grade']}（盲測失敗 {nfail}/3）")
        out.append(m)

    emit("")
    emit("═" * 104)
    emit("判詞鏈（全部由上面計算得出，非預設）")
    surv = [m for m in out if m["grade"] == "SURVIVE"]
    sat = [m for m in out if m["grade"] == "SATURATED"]
    for m in out:
        emit(f"  {m['label']:34s} ⟹ {m['grade']}")
    emit("")
    if sat:
        for m in sat:
            emit(f"  ★ 唯一未被否證者 = {m['label']}，但 dof = {m['dof']}（{m['n_classes']} 類 vs 4 錨）")
            emit(f"    ⟹ 錨已耗盡：模型中每個參數都被一個錨釘死 ⟹ 無殘餘自由度可作預言")
            for ch, v in m["verdicts"].items():
                if v["O"] is not None:
                    emit(f"    盲測 {ch:8s}：⟨O⟩ = {v['O']:.6f}（強制） em = {v['em']:+.4f} MeV "
                         f"t = {v['t']:+.4f} ⟹ {'可行' if v['feasible'] else '不可行'}")
                else:
                    emit(f"    盲測 {ch:8s}：類別 {v['category']} 新 ⟹ **無預言**")
    else:
        emit("  ★ 無模型存活 ⟹ 連類別層亦不成立（需要更強之否定陳述）")

    emit("")
    emit("接口必要性（本輪結果 ＋ 已立事實）：")
    emit("  L1 I2 必需：無 iso-breaking 頂點 ⟹ ∂M/∂(m_d−m_u) ≡ 0（M5 docs/43④／45③；m_k=M₀γ_k flavour-blind）⟹ ⟨O⟩≡0，與 4 錨矛盾")
    mL2 = [m for m in out if m["label"].startswith("A")][0]
    emit(f"  L2 I3 必需：A（通道無關）⟹ 錨 z={mL2['z']:.2f}（{mL2['grade_anchor']}）"
         + (f"＋盲測失敗 {mL2['blind_fails']}/3" if mL2["blind_fails"] else ""))
    mB = [m for m in out if m["label"].startswith("B")][0]
    mC = [m for m in out if m["label"].startswith("C")][0]
    emit(f"  L3 ns 依賴必需：B（只記 nq）⟹ 錨 z={mB['z']:.2f}（{mB['grade_anchor']}）"
         + (f"＋盲測失敗 {mB['blind_fails']}/3" if mB["blind_fails"] else ""))
    emit(f"  L4 nq 依賴：C（只記 ns）⟹ 錨 z={mC['z']:.2f}（{mC['grade_anchor']}）"
         + (f"；但盲測失敗 {mC['blind_fails']}/3 ⟹ {mC['grade']}" if mC["blind_fails"] else "；盲測全過"))

    res = dict(anchors=[dict(name=n, nq=nq, ns=ns, O=o, sigma=s) for n, nq, ns, o, s in A],
               models=out, k0=K0, X=X)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    emit("")
    emit(f"[JSON] {OUT}")
    emit(f"[LINES] {len(log)}")


if __name__ == "__main__":
    main()