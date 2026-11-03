#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — M11 文檔交叉核（唯讀；idcm-audit-verify §2.4；模板之填版）
================================================================================
R1 唯讀：只讀 doc，只寫 data/l3_doc_crosscheck.json ＋ data/run_l3.log。
R2 期望值**只來自 L2**（`data/l2_independent_verify.json` 之 `_derived`），零抄 docs。
R3 容差依文檔有效位；R4 負值取絕對值；R5 必要時 `tol=` 收緊。
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG, RES, NFAIL = [], {}, 0

CLAIMS = [
    # docs/01（R1 判準化）—— 錨值與 k₀
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "k0_from_a0", "label": "k₀ = αħc/a₀ = 3.815262"},
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "hc_over_a0", "label": "ħc/a₀ = 522.8283"},
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "O_n-p", "label": "⟨O⟩_n−p = 1.1263"},
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "O_K+-K0", "label": "⟨O⟩_K = 2.4213"},
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "O_S--S0", "label": "⟨O⟩_Σ = 1.9090"},
    {"doc": "docs/01_r1_criteria_and_target.md", "key": "O_Xi--Xi0", "label": "⟨O⟩_Ξ = 2.1719"},
    # docs/06（k₀ 溯源）
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "k0_from_a0", "label": "k₀ 3.815262"},
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "hc_over_a0", "label": "522.828"},
    # docs/06（k₀^crit 四模型）—— 本輪由 L2 補算，獨立核
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "k0_crit_D", "label": "k₀^crit(D) = 2.5425"},
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "k0_crit_A", "label": "k₀^crit(A) = 4.0730"},
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "k0_crit_B", "label": "k₀^crit(B) = 4.1358"},
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "k0_crit_C", "label": "k₀^crit(C) = 5.3358"},
    {"doc": "docs/06_r5_k0_identity_and_sensitivity.md", "key": "a_crit_all", "label": "a^crit(all) = 0.26987"},
    {"doc": "docs/05_r4_necessity_and_saturation.md", "key": "a_crit_all", "label": "docs/05 引 a^crit 0.26987"},
    {"doc": "01_problem_record.md", "key": "k0_crit_C", "label": "R5 行 k₀^crit C 5.3358"},
    # docs/07（M5 docs/52 對接）—— T₂ 兩支
    {"doc": "docs/07_m5_docs52_xref.md", "key": "T2_aeff_over_k0", "label": "T₂ = A_eff/k₀ = 1.166158", "tol": 1e-5},
    {"doc": "docs/07_m5_docs52_xref.md", "key": "T2_v2_branch", "label": "T₂(v2 支) = 1.187051", "tol": 1e-5},
    # docs/09（em 錨溯源）
    {"doc": "docs/09_em_anchor_provenance.md", "key": "em_eps_flag45", "label": "em(ε=0.79) = 2.2773"},
    {"doc": "docs/09_em_anchor_provenance.md", "key": "sig_em_flag45", "label": "σ_em = 0.0763"},
    {"doc": "docs/09_em_anchor_provenance.md", "key": "em_eps_milc", "label": "em(ε=0.78) = 2.2646"},
    {"doc": "docs/09_em_anchor_provenance.md", "key": "z_flag45_only", "label": "z(只 ε) = −25.47"},
    {"doc": "docs/09_em_anchor_provenance.md", "key": "z_flag45_total", "label": "z(ε+X) = −11.16"},
    # docs/10（G0）
    {"doc": "docs/10_g0_x_operator_verdict.md", "key": "G0_dM_D_hf", "label": "δM(D) = 0.01811"},
    {"doc": "docs/10_g0_x_operator_verdict.md", "key": "G0_deficit_D", "label": "deficit = 141.3×"},
    # docs/11（路線 C）
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "routeC_a", "label": "a = 5.6581"},
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "routeC_b", "label": "b = 3.1793"},
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "routeC_chi2_dof", "label": "χ²/dof = 3154"},
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "routeC_O_eff", "label": "⟨O⟩_eff = 1.2419"},
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "degen_np_S0Splus", "label": "退化差 1.9787 MeV"},
    {"doc": "docs/11_route_c_total_dm_verdict.md", "key": "degen_S_Xi_D", "label": "退化差 2.0430 MeV"},
    # 01_problem_record（時序／開放問題表）
    {"doc": "01_problem_record.md", "key": "k0_from_a0", "label": "R6 行 k₀ 3.815262"},
    {"doc": "01_problem_record.md", "key": "routeC_chi2_dof", "label": "Q16 行 χ²/dof 3154"},
    {"doc": "01_problem_record.md", "key": "G0_deficit_D", "label": "Q15 行 141"},
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


def main():
    hdr("§0 讀 L2（期望值唯一來源）")
    if not L2PATH.exists():
        print(f"✗ 缺 {L2PATH}")
        return 2
    l2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    der = l2.get("_derived", {})
    chk("E1 L2 有 _derived 且非空", bool(der), f"{len(der)} keys")

    hdr("§1 逐 claim 交叉核（唯讀）")
    cache, by_doc = {}, {}
    for c in CLAIMS:
        doc, key, label = c["doc"], c["key"], c.get("label", c["key"])
        p = ROOT / doc
        if not p.exists():
            chk(f"C.{doc} 存在", False, "檔案不存在")
            continue
        if doc not in cache:
            cache[doc] = p.read_text(encoding="utf-8")
        if key not in der:
            chk(f"C.{doc}.{label}（L2 key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1
        exp = der[key].get("float")
        if exp is None:
            chk(f"C.{doc}.{label}", False, "L2 值超 float 範圍")
            continue
        ok, got, rel = match(tokens(cache[doc]), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{Path(doc).name[:18]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok
            else f"✗ L2 {exp:.6g}，文檔最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率")
    for d, (h, t) in by_doc.items():
        print(f"    {d:44s} {h}/{t}")
        chk(f"F.{Path(d).name} 全命中（{h}/{t}）", h == t)

    hdr("§3 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "by_doc": by_doc}
    (DATA / "l3_doc_crosscheck.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2),
                                                 encoding="utf-8")
    print(f"\n寫出: {DATA/'l3_doc_crosscheck.json'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())