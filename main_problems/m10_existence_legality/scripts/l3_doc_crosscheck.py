#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — 文檔交叉核（M10；idcm-audit-verify §2.4）
============================================================
鐵律：R1 唯讀（只寫 data/ 兩個輸出）；R2 期望值**只來自 L2** `_derived`；R3 容差依有效位；
      R4 負值取絕對值；R5 容差要分辨真錯；R6 細整數必須帶 anchor。
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG, RES = [], {}
NFAIL = 0


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


CLAIMS = [
    # ---- 判準與合法集
    {"doc": "docs/09_channel_completion.md", "key": "P_size", "label": "合法集 90", "anchor": "𝒞_v2合法集", "span": 40},
    {"doc": "docs/09_channel_completion.md", "key": "A_max", "label": "A_max 92", "anchor": "合法集=90條；總收縮33.8%；`A_max`", "span": 40},
    # ---- 檢定
    {"doc": "docs/10_falsification_nubase.md", "key": "SP_size", "label": "S\\P = 182", "anchor": "S\\P"},
    {"doc": "docs/10_falsification_nubase.md", "key": "tier3", "label": "Tier III = 146", "anchor": "**TierIII**（rel<1e−4）|", "span": 30},
    {"doc": "docs/10_falsification_nubase.md", "key": "radius_rel", "label": "失效半徑 2.84e−4", "anchor": "max rel"},
    {"doc": "docs/10_falsification_nubase.md", "key": "span", "label": "時間盲 8.3e53", "anchor": "跨度", "span": 30, "tol": 6e-3},
    {"doc": "docs/12_primordial_label.md", "key": "radius_rel", "label": "丰度標籤同半徑", "anchor": "max rel"},
    {"doc": "docs/17_surviving_results.md", "key": "radius_rel", "label": "S2.4 半徑", "anchor": "S2.4"},
    {"doc": "docs/17_surviving_results.md", "key": "tier3", "label": "X3 之 Tier III 條數", "anchor": "146 條之", "span": 30},
    # ---- 機制分類
    {"doc": "docs/15_mechanism_axis.md", "key": "split_C3", "label": "𝒞3 位壘 126", "anchor": "**𝒞3**", "span": 60},
    {"doc": "docs/15_mechanism_axis.md", "key": "split_C2", "label": "𝒞2 二階弱 51", "anchor": "**𝒞2**", "span": 40},
    {"doc": "docs/15_mechanism_axis.md", "key": "split_C1sel", "label": "𝒞1-禁戒 5", "anchor": "𝒞1-禁戒", "span": 60},
    {"doc": "docs/15_mechanism_axis.md", "key": "split_unexpl", "label": "未解釋 0", "anchor": "無法解釋", "span": 30},
    # ---- W1 約定閘
    {"doc": "docs/19_em_sign_downstream.md", "key": "em_ident", "label": "em(p)−em(n) = 1.2667", "anchor": "1.2667"},
    {"doc": "docs/19_em_sign_downstream.md", "key": "x_n", "label": "x(n) = 0.7823", "anchor": "0.7823"},
    {"doc": "docs/20_status_audit.md", "key": "A_max_v1", "label": "A_max(v1) 1.27175", "anchor": "1.27175"},
    {"doc": "docs/20_status_audit.md", "key": "A_max_v2", "label": "A_max(v2) 1.24937", "anchor": "1.24937"},
    {"doc": "docs/19_em_sign_downstream.md", "key": "gap_v2", "label": "缺口 0.01733", "anchor": "缺口"},
    {"doc": "docs/19_em_sign_downstream.md", "key": "sigma_gap_v2", "label": "0.184σ", "anchor": "0.184"},
    {"doc": "docs/21_r21_own_residuals.md", "key": "suppression_target", "label": "抑制目標 0.276548", "anchor": "0.276548"},
    {"doc": "docs/21_r21_own_residuals.md", "key": "fam_hits_1pct", "label": "±1% 命中 5", "anchor": "±1%"},
    {"doc": "docs/21_r21_own_residuals.md", "key": "fam_size", "label": "式族 4538", "anchor": "4538"},
    {"doc": "docs/20_status_audit.md", "key": "gamma_over_v2", "label": "γ 超 v2 11.76%", "anchor": "11.76", "scale": 100, "span": 20},
]

TOK = re.compile(r"(\d+(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")


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
    if t == 0:
        hit = next((v for v, _ in toks if v == 0), None)
        return (hit is not None), hit, 0.0
    for v, sig in toks:
        if v == 0:
            continue
        rel = abs(abs(v) - t) / t
        tol = tol_override if tol_override else (1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3))
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


def window_around(text, anchor, span=90):
    flat = re.sub(r"[\s\u00a0]+", "", text)
    a = re.sub(r"[\s\u00a0]+", "", str(anchor))
    i = flat.find(a)
    if i < 0:
        return ""
    return flat[max(0, i - 20): i + len(a) + span]


def main():
    hdr("§0 讀 L2")
    if not L2PATH.exists():
        print("✗ 缺 L2 輸出，先跑 scripts/l2_independent_verify.py")
        return 2
    l2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    der = l2.get("_derived", {})
    chk("E1 L2 _derived 非空", bool(der), f"{len(der)} keys")

    hdr("§1 逐 claim（唯讀）")
    cache, by_doc = {}, {}
    for c in CLAIMS:
        doc, key, label = c["doc"], c["key"], c.get("label", c["key"])
        p = ROOT / doc
        if not p.exists():
            chk(f"C.{doc} 存在", False, "檔案不存在")
            continue
        cache.setdefault(doc, p.read_text(encoding="utf-8"))
        if key not in der:
            chk(f"C.{label}（L2 key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1
        exp = der[key]["float"] * float(c.get("scale", 1))
        txt = cache[doc]
        if c.get("anchor"):
            txt = window_around(cache[doc], c["anchor"], c.get("span", 90))
            if not txt:
                chk(f"C.{label}", False, f"錨點「{c['anchor']}」不存在")
                continue
        ok, got, rel = match(tokens(txt), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{doc.split('/')[-1][:18]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok
            else f"✗ L2 {exp:.6g}，最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率")
    for d, (h, t) in sorted(by_doc.items()):
        print(f"    {d:46s} {h}/{t}")
        chk(f"F.{d} 全命中", h == t, f"{h}/{t}")

    hdr("§3 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    print(f"{n - NFAIL}/{n} PASS（L3）")
    LOG.append(f"{n - NFAIL}/{n} PASS（L3）")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "by_doc": by_doc}
    RES["pass"] = n - NFAIL
    RES["fail"] = NFAIL
    (DATA / "l3_doc_crosscheck.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l3.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())