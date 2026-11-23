#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — M6 文檔交叉核（idcm-audit-verify §2.4）
================================================================================
唯讀；期望值**只取自 L2**（data/l2_independent_verify.json）。禁抄 docs 數字。
★ 千分位整數（如 4,862,007,014,624,589,824）唔可作數值 claim —— token 會被逗號切碎 ⟹
  該項改以「容量 2511.082350 MeV」作代表（同一結構量，L2 key capacity_MeV）。
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG, RES, NFAIL = [], {}, 0

CLAIMS = [
    {"doc": "docs/14_scale_axis_theorem.md", "key": "sum_gamma", "label": "Σγ = 90（環內精確）",
     "anchor": "Fraction(90", "span": 40},
    {"doc": "docs/14_scale_axis_theorem.md", "key": "TA_v2", "label": "T_A v2 = 18.408939（附註 B3′）",
     "anchor": "18.408939", "span": 30},
    {"doc": "docs/13_h2_verdict.md", "key": "zero_root", "label": "零交叉連續根 18.59330797",
     "anchor": "18.59330797", "span": 40},
    {"doc": "docs/13_h2_verdict.md", "key": "g18", "label": "g(18) = +0.17329545",
     "anchor": "g(18)=+0.17329545", "span": 40},
    {"doc": "docs/07_hidden_layer_mechanism.md", "key": "nq_vsm", "label": "n_q(v_SM) = 99",
     "anchor": "8824.7859", "span": 30},
    {"doc": "docs/07_hidden_layer_mechanism.md", "key": "capacity_MeV", "label": "容量 2511.082350 MeV",
     "anchor": "2511.082350", "span": 30},
    {"doc": "docs/07_hidden_layer_mechanism.md", "key": "lam_max_1", "label": "λ_max(1) = 4 = 4k",
     "anchor": "λ_max(1) = γ₄₄ = 4 = 4k", "span": 40},
    {"doc": "docs/15_m_series_handoff_adjudication.md", "key": "kappa_lo", "label": "κ/M₀ 下界 2.457611",
     "anchor": "2.457611", "span": 40},
    {"doc": "docs/15_m_series_handoff_adjudication.md", "key": "kappa_hi", "label": "κ/M₀ 上界 2.529465",
     "anchor": "2.529465", "span": 40},
    {"doc": "docs/15_m_series_handoff_adjudication.md", "key": "Tt_v2", "label": "T_t v2 = 0.06347820493",
     "anchor": "0.06347820493", "span": 30},
    {"doc": "docs/15_m_series_handoff_adjudication.md", "key": "TA_v1", "label": "T_A v1 = 18.738758425",
     "anchor": "18.738758425", "span": 30},
    {"doc": "final_result.md", "key": "H0_z", "label": "H₀ z = 0.087σ",
     "anchor": "0.087σ", "span": 40},
    {"doc": "01_problem_record.md", "key": "n_strict", "label": "支撐定理 43 嚴格",
     "anchor": "43/43", "span": 60},
    {"doc": "docs/08_framework_self_audit.md", "key": "rank_Q", "label": "rank_ℚ{γ_n} = 20",
     "anchor": "rank = 20", "span": 40},
    {"doc": "docs/09_upgrade_routes.md", "key": "n_orbits", "label": "7 個 Galois 軌道",
     "anchor": "7 個軌道", "span": 40},
    {"doc": "docs/10_forced_single_scale.md", "key": "n_orbits", "label": "7 個軌道（T2′ 上界）",
     "anchor": "7 個軌道", "span": 40},
    {"doc": "docs/10_forced_single_scale.md", "key": "max_orbit", "label": "最大軌道 20",
     "anchor": "20/10/5/5/2/1/1", "span": 20},
    {"doc": "docs/12_economy_and_homogeneity.md", "key": "rank_Q", "label": "rank 20（參數經濟）",
     "anchor": "rank_ℚ{γ_n} = 20", "span": 40},
    {"doc": "final_result.md", "key": "sum_gamma_sq", "label": "Σγ² = 272",
     "anchor": "Σγ² = 272", "span": 30},
    {"doc": "docs/11_beta_no_go.md", "key": "delta_alpha_inv", "label": "Δα⁻¹ = 9 = (B−1)²",
     "anchor": "Δα⁻¹ = N_gen²", "span": 40},
    {"doc": "README.md", "key": "sum_gamma", "label": "Σγ = 90（上游依賴列）",
     "anchor": "Σγ=90", "span": 20},
    {"doc": "docs/16_full_chain_audit.md", "key": "rank_Q", "label": "rank_ℚ = 20（本輪精確算出）",
     "anchor": "rank_ℚ{γ_n} = 20 本輪首次", "span": 40},
]

TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def tokens(text):
    out = []
    for m in TOK.finditer(text.replace("\u2212", "-")):
        mant, exp = m.group(1), m.group(2)
        try:
            v = float(mant) * (10.0 ** int(exp)) if exp else float(mant)
        except (ValueError, OverflowError):
            continue
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        out.append((v, sig))
    return out


def match(toks, target, tol_override=None):
    t = abs(target)
    best, best_d = None, float("inf")
    for v, sig in toks:
        if v == 0:
            continue
        rel = abs(abs(v) - t) / t
        tol = tol_override if tol_override else (1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3))
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


def window_around(text, anchor, span=70):
    flat = re.sub(r"[\s\u00a0]+", "", text)
    a = re.sub(r"[\s\u00a0]+", "", str(anchor))
    i = flat.find(a)
    if i < 0:
        return ""
    return flat[max(0, i - 12): i + len(a) + span]


def main():
    hdr("§0 讀 L2（期望值唯一來源）")
    if not L2PATH.exists():
        print("✗ 缺 data/l2_independent_verify.json —— 先跑 scripts/l2_independent_verify.py")
        return 2
    l2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    der = l2.get("_derived", {})
    chk("E1 L2 有 _derived 且非空", bool(der), f"keys = {sorted(der)[:6]}...")

    hdr("§1 逐 claim 交叉核（唯讀）")
    cache, by_doc = {}, {}
    for c in CLAIMS:
        doc, key, label = c["doc"], c["key"], c.get("label", c["key"])
        p = ROOT / doc
        if not p.exists():
            chk(f"C.{doc} 存在", False, "檔案不存在")
            continue
        cache.setdefault(doc, p.read_text(encoding="utf-8"))
        if key not in der:
            chk(f"C.{doc}.{label}（L2 key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1
        exp = der[key].get("float")
        if exp is None:
            chk(f"C.{doc}.{label}", False, "L2 值超出 float 範圍")
            continue
        txt = cache[doc]
        if c.get("anchor"):
            txt = window_around(cache[doc], c["anchor"], c.get("span", 70))
            if not txt:
                chk(f"C.{doc.split('/')[-1][:18]}.{label}", False, f"✗ 錨點「{c['anchor']}」不存在")
                continue
        ok, got, rel = match(tokens(txt), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{doc.split('/')[-1][:18]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok
            else f"✗ L2 {exp:.6g}，文檔最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率")
    for d, (h, t) in by_doc.items():
        print(f"    {d:52s} {h}/{t}")
        chk(f"F.{d} 全命中（{h}/{t}）", h == t)

    hdr("§3 匯總")
    n = len(RES)
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    _pas = n - NFAIL
    print(f"{_pas}/{n} PASS（L3）")
    LOG.append(f"{_pas}/{n} PASS（L3）")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "by_doc": by_doc}
    RES["pass"] = n - NFAIL
    RES["fail"] = NFAIL
    (DATA / "l3_doc_crosscheck.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\n寫出: data/l3_doc_crosscheck.json\ndata/run_l3.log")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())