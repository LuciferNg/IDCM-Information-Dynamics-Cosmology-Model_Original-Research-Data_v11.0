#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — M13 文檔交叉核（idcm-audit-verify §2.4；唯讀）

鐵律
  R1 唯讀：對 docs 只 open 讀；只寫 data/ 下兩個輸出。
  R2 期望值**只來自 L2**（`_derived[key]`）；不抄 docs 內數字。
  R3 容差依文檔有效位：≥6 位 → 1e−4；≥4 位 → 5e−4；其餘 → 5e−3。
  R4 負值以絕對值比（文檔用 U+2212）。
  R5 每個 headline 數字都要有 claim；逐檔命中率必須 100%。
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
    """數值 token：處理 U+2212、千分位逗號（僅 3 位組）、及 a/b 簡單分數。"""
    t = text.replace("\u2212", "-")
    t = re.sub(r"(?<=\d),(?=\d{3}(?!\d))", "", t)      # 105,952 → 105952
    out = []
    for m in re.finditer(r"(?<![\d.])(\d+)\s*/\s*(\d+)(?![\d.])", t):
        a, b = int(m.group(1)), int(m.group(2))
        if b:
            sig = max(len(str(a)), len(str(b)))
            out.append((a/b, min(sig, 6)))
    for m in re.finditer(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?", t):
        mant, exp = m.group(1), m.group(2)
        try:
            v = float(mant)*(10.0**int(exp)) if exp else float(mant)
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
        rel = abs(abs(v) - t)/t
        tol = tol_override if tol_override else (1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3))
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


# ============================================================================
# CLAIMS —— (doc, L2 key)；每個 headline 數字一條
# ============================================================================
CLAIMS = [
    # docs/01 重現閘
    {"doc": "docs/01_reproduction_gate.md", "key": "rho_d", "label": "ρ_d"},
    {"doc": "docs/01_reproduction_gate.md", "key": "rho_u", "label": "ρ_u（負）"},
    {"doc": "docs/01_reproduction_gate.md", "key": "sigma_rho_d", "label": "σ_ρd"},
    {"doc": "docs/01_reproduction_gate.md", "key": "delta_s", "label": "δ_s"},
    {"doc": "docs/01_reproduction_gate.md", "key": "dX", "label": "Δ_Ξ"},
    # docs/02 空容量（枚舉類容量）
    {"doc": "docs/02_null_capacity.md", "key": "nclass_K2", "label": "K=2 類大小"},
    {"doc": "docs/02_null_capacity.md", "key": "nclass_K3", "label": "K=3 類大小"},
    {"doc": "docs/02_null_capacity.md", "key": "nclass_K4", "label": "K=4 類大小"},
    {"doc": "docs/02_null_capacity.md", "key": "nclass_K5", "label": "K=5 類大小"},
    {"doc": "docs/02_null_capacity.md", "key": "nclass_K6", "label": "K=6 類大小"},
    # docs/03 閉包掃描（A₂ Cartan）
    {"doc": "docs/03_closure_scan.md", "key": "sigma_rho_u", "label": "σ_ρu"},
    {"doc": "docs/03_closure_scan.md", "key": "rd_int3", "label": "ρ_d(α_s=3)=0.5"},
    # docs/04 候選
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rd_int1", "label": "α_s=1 → 1.5"},
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rd_int2", "label": "α_s=2 → 0.75"},
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rd_int4", "label": "α_s=4 → 0.375"},
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rd_int5", "label": "α_s=5 → 0.3"},
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rd_int8", "label": "α_s=8 → 0.1875"},
    {"doc": "docs/04_carry_by_2_candidate.md", "key": "rho_d", "label": "ρ_d 觀測"},
    # docs/05 普遍性（α_s/α_u 窗；線性傳播約定）
    {"doc": "docs/05_universality.md", "key": "alpha_ratio", "label": "α_s/α_u 中心"},
    {"doc": "docs/05_universality.md", "key": "alpha_ratio_lin_lo", "label": "窗下界（線性）"},
    {"doc": "docs/05_universality.md", "key": "alpha_ratio_lin_hi", "label": "窗上界（線性）"},
    # docs/06 樹層核
    {"doc": "docs/06_ktree_alpha_probe.md", "key": "rho_u_a1416", "label": "α=(1,4,16) 之 ρ_u"},
    {"doc": "docs/06_ktree_alpha_probe.md", "key": "rho_d_a1416", "label": "α=(1,4,16) 之 ρ_d"},
    # docs/07 甲路線（216 試驗）
    {"doc": "docs/07_alpha_source.md", "key": "n_trials_alpha", "label": "N_trials = 216"},
    # docs/08 精度預算
    {"doc": "docs/08_precision_budget.md", "key": "delta_s", "label": "δ_s"},
    {"doc": "docs/08_precision_budget.md", "key": "sigma_floor_sigma", "label": "Σ 硬底 σ_ρd ≥"},
    {"doc": "docs/08_precision_budget.md", "key": "dS", "label": "Δ_Σ"},
    {"doc": "docs/08_precision_budget.md", "key": "dX", "label": "Δ_Ξ"},
    {"doc": "docs/08_precision_budget.md", "key": "sigma_rho_d", "label": "σ_ρd"},
    # docs/09 完備性審計（家族掃描）
    {"doc": "docs/09_completeness_audit.md", "key": "chi2_x2", "label": "χ²(x=2)"},
    {"doc": "docs/09_completeness_audit.md", "key": "chi2_x3", "label": "χ²(x=3)"},
    {"doc": "docs/09_completeness_audit.md", "key": "chi2_x4", "label": "χ²(x=4)"},
    {"doc": "docs/09_completeness_audit.md", "key": "z_r15_4", "label": "z(r=15/4)"},
    {"doc": "docs/09_completeness_audit.md", "key": "z_r11_3", "label": "z(r=11/3)"},
    {"doc": "docs/09_completeness_audit.md", "key": "z_r7_2", "label": "z(r=7/2)"},
    # docs/10 扇區
    {"doc": "docs/10_sector_extension.md", "key": "delta_s", "label": "δ_s"},
    {"doc": "docs/10_sector_extension.md", "key": "delta_d", "label": "δ_d"},
    # docs/11 T1 定理
    {"doc": "docs/11_t1_provenance_theorem.md", "key": "z_konly_rho_d", "label": "α 依 k → z=7.74"},
    {"doc": "docs/11_t1_provenance_theorem.md", "key": "z_donly_rho_d", "label": "α 依 d → z=−9.01"},
    # docs/12 競爭者格
    {"doc": "docs/12_rival_lattice_audit.md", "key": "rd_r15_4", "label": "ρ_d(15/4)=0.4"},
    {"doc": "docs/12_rival_lattice_audit.md", "key": "rd_r11_3", "label": "ρ_d(11/3)=9/22"},
    {"doc": "docs/12_rival_lattice_audit.md", "key": "rd_r7_2", "label": "ρ_d(7/2)=3/7"},
    {"doc": "docs/12_rival_lattice_audit.md", "key": "z_r15_4", "label": "z(15/4)"},
    # docs/13 迹閉合
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "trace", "label": "迹 ρ_u+ρ_d"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "sigma_trace", "label": "σ_迹"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "corr", "label": "corr（負）"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "ratio_dN_dX", "label": "Δ_N/Δ_Ξ"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "target_3_16", "label": "3/16"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "diff", "label": "ρ_d−ρ_u"},
    # docs/14 四閘驗證
    {"doc": "docs/14_hypothesis_verification.md", "key": "trace_sigma_chan", "label": "迹（Σ 通道）"},
    {"doc": "docs/14_hypothesis_verification.md", "key": "z_trace_sigma_chan", "label": "z（Σ 通道）"},
    {"doc": "docs/14_hypothesis_verification.md", "key": "z_trace", "label": "z（Ξ 通道）"},
    {"doc": "docs/14_hypothesis_verification.md", "key": "need_sigmaX0_3sigma_id", "label": "識別所需 σ(Ξ⁰)"},
    {"doc": "docs/14_hypothesis_verification.md", "key": "dev_over_base_q16", "label": "dev/基準"},
    {"doc": "docs/14_hypothesis_verification.md", "key": "dev_vs_3_8", "label": "絕對 dev"},
    # docs/15 審計報告（headline = 審計四數 ＋ 兩項真錯之值）
    {"doc": "docs/15_audit_report.md", "key": "rho_d", "label": "ρ_d（L2↔doc）"},
    {"doc": "docs/15_audit_report.md", "key": "sigma_rho_d", "label": "σ_ρd"},
    {"doc": "docs/15_audit_report.md", "key": "corr", "label": "corr（負）"},
    {"doc": "docs/15_audit_report.md", "key": "dev_over_base_q16", "label": "更正後 dev/基準 3.571"},
    {"doc": "docs/15_audit_report.md", "key": "need_sigmaX0_3sigma_id", "label": "識別所需 σ(Ξ⁰)"},
    {"doc": "docs/15_audit_report.md", "key": "target_3_8", "label": "3/8"},
    # docs/16 傳播審計
    {"doc": "docs/16_propagation_audit.md", "key": "rho_d", "label": "ρ_d（canon）"},
    {"doc": "docs/16_propagation_audit.md", "key": "trace", "label": "迹 0.377615281752"},
    {"doc": "docs/16_propagation_audit.md", "key": "dev_over_base_q16", "label": "3.571"},
    {"doc": "docs/16_propagation_audit.md", "key": "sigma_rho_d", "label": "σ_ρd"},
    {"doc": "docs/16_propagation_audit.md", "key": "need_sigmaX0_3sigma_id", "label": "0.0158"},
    # docs/17 handoff 回覆
    {"doc": "docs/17_handoff_replies.md", "key": "ratio_dN_dX", "label": "Δ_N/Δ_Ξ"},
    {"doc": "docs/17_handoff_replies.md", "key": "sigma_ratio_dN_dX", "label": "σ(Δ_N/Δ_Ξ)"},
    {"doc": "docs/17_handoff_replies.md", "key": "trace", "label": "迹 = 2q"},
    {"doc": "docs/17_handoff_replies.md", "key": "z_konly_rho_d", "label": "α 依 k 7.74σ"},
    {"doc": "docs/17_handoff_replies.md", "key": "z_donly_rho_d", "label": "α 依 d −9.01σ"},
    {"doc": "docs/17_handoff_replies.md", "key": "nclass_K5", "label": "K=5 類大小"},
    {"doc": "docs/17_handoff_replies.md", "key": "alpha_ratio", "label": "r = α_s/α_u"},
    {"doc": "docs/17_handoff_replies.md", "key": "alpha_sum_ratio", "label": "α_s/(α_u+α_d)"},
    {"doc": "docs/17_handoff_replies.md", "key": "sigma_alpha_sum_ratio", "label": "σ(α 組合)"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "alpha_sum_ratio", "label": "α 組合（附註5）"},
    {"doc": "docs/13_symmetric_trace_closure.md", "key": "sigma_alpha_sum_ratio", "label": "σ（附註5）"},
    # final_result（階段性收尾，2026-11-03）
    {"doc": "final_result.md", "key": "alpha_ratio", "label": "α_s/α_u 中心"},
    {"doc": "final_result.md", "key": "rho_d", "label": "ρ_d"},
    {"doc": "final_result.md", "key": "trace", "label": "迹"},
    {"doc": "final_result.md", "key": "ratio_dN_dX", "label": "Δ_N/Δ_Ξ"},
    {"doc": "final_result.md", "key": "alpha_sum_ratio", "label": "α 組合"},
    {"doc": "final_result.md", "key": "sigma_alpha_sum_ratio", "label": "σ(α 組合)"},
    {"doc": "final_result.md", "key": "sigma_rho_d", "label": "σ_ρd"},
    {"doc": "final_result.md", "key": "sigma_trace", "label": "σ_迹"},
    # docs/18 全推（R1／R2／R3）
    {"doc": "docs/18_full_push.md", "key": "trace", "label": "迹（Ξ 通道）"},
    {"doc": "docs/18_full_push.md", "key": "Xi_c_split", "label": "Ξc⁰−Ξc⁺"},
    {"doc": "docs/18_full_push.md", "key": "Xi_c_sigma", "label": "σ(Ξc 分裂)"},
    {"doc": "docs/18_full_push.md", "key": "sigma_Xi0_needed_B5", "label": "σ(Ξ⁰) 需 ≤（vs B=5）"},
]


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
            chk(f"C.{Path(doc).stem[:16]}.{label}（key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1
        exp = der[key].get("float")
        if exp is None:
            chk(f"C.{Path(doc).stem[:16]}.{label}", False, "L2 值超 float")
            continue
        ok, got, rel = match(tokens(cache[doc]), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{Path(doc).stem[:16]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok else
            f"✗ L2 {exp:.6g}；文檔最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率（必須 100%）")
    for d, (h, t) in sorted(by_doc.items()):
        print(f"    {d:44s} {h}/{t}")
        chk(f"F.{Path(d).stem} 全命中（{h}/{t}）", h == t)

    hdr("§3 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": n, "fail": NFAIL,
                       "by_doc": {d: f"{h}/{t}" for d, (h, t) in by_doc.items()}}
    (DATA / "l3_doc_crosscheck.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / 'l3_doc_crosscheck.json'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())