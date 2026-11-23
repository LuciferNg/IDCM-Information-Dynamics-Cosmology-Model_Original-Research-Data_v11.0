#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/13 — 盲測後診斷（**非判據**）：族內最近目標之值
================================================================================
· 族由 `scripts/12_count_family_blindtest.py` 之**同一 code path** 重建（import 其
  `ATOMS`/`build_D`/`maps`/`composites`），並自核重建值數 = scripts/12 記錄之 14442 ⟹ 同族。
· 目的：報「最近幾個」與跨度，供人判斷族有幾貼近；**不參與任何判詞**（docs/10 §2.4 判詞已於跑前寫死）。
· 唯讀族定義；零寫入族外檔案。

執行: /usr/bin/python3 scripts/13_nearest_diag.py
輸出: data/13_nearest_diag.json、data/run_13.log
"""
import importlib.util
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
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


def load_m12():
    spec = importlib.util.spec_from_file_location(
        "m12_frozen", HERE / "12_count_family_blindtest.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # 只有 def/常數，__main__ 守衛 ⟹ 不會執行
    return mod


def build(module):
    D = module.build_D()
    bases = ([(f"{v}[{lab}]", mp.mpf(v)) for v, lab in module.ATOMS]
             + [(lab, mp.mpf(v)) for lab, v, _ in D])
    entries = []
    for lab, x in bases:
        L = mp.log10(x)
        for ml, Lv in module.maps(L, x):
            if abs(Lv) <= 300:
                entries.append((f"{lab}·{ml}", Lv))
    for lab, x in [(lab, mp.mpf(v)) for lab, v, _ in D]:
        L = mp.log10(x)
        for cl, Lv in module.composites(L, x):
            if abs(Lv) <= 300:
                entries.append((f"{lab}·{cl}", Lv))
    entries.sort(key=lambda t: t[1])
    uniq = []
    for lab, Lv in entries:
        if uniq and abs(Lv - uniq[-1][1]) < mp.mpf("1e-30"):
            continue
        uniq.append((lab, Lv))
    return D, uniq


def main():
    m = load_m12()
    ref = json.loads((DATA / "12_count_family.json").read_text(encoding="utf-8"))
    MU = mp.mpf(ref["_target"]["mu"])
    L_mu = mp.log10(MU)

    hdr("§1 由 scripts/12 同一 code path 重建族（自核同族）")
    D, uniq = build(m)
    print(f"  |D| = {len(D)}｜唯一值 = {len(uniq)}")
    chk("G1 重建值數 = scripts/12 記錄（14442）", len(uniq) == ref["_family"]["n_unique"],
        f"{len(uniq)} vs {ref['_family']['n_unique']}")
    chk("G2 |D| = 1368", len(D) == ref["_family"]["n_D"], f"{len(D)}")

    hdr("§2 最近 10 個（診斷，非判據）")
    sc = sorted(((abs(mp.power(10, Lv - L_mu) - 1), lab, Lv) for lab, Lv in uniq),
                key=lambda t: t[0])
    top = []
    for dev, lab, Lv in sc[:10]:
        top.append({"expr": lab, "dev": mp.nstr(dev, 8), "dev_pct": mp.nstr(dev * 100, 6),
                    "log10_delta": mp.nstr(Lv - L_mu, 10), "log10_value": mp.nstr(Lv, 20)})
        print(f"  dev = {mp.nstr(dev * 100, 5):>9}%  log₁₀Δ = {mp.nstr(Lv - L_mu, 8):>12}  {lab}")
    chk("G3 最近者 dev > 1e−3（與 scripts/12 之 0 命中一致）",
        mp.mpf(sc[0][0]) > mp.mpf("1e-3"), f"最近 dev = {mp.nstr(sc[0][0]*100, 5)}%")

    hdr("§3 族跨度")
    span = [mp.nstr(uniq[0][1], 8), mp.nstr(uniq[-1][1], 8)]
    print(f"  log₁₀ ∈ [{span[0]}, {span[1]}]")
    chk("G4 目標落於族跨度內（唔係外插）", mp.mpf(span[0]) < L_mu < mp.mpf(span[1]),
        f"目標 log₁₀ = {mp.nstr(L_mu, 10)}")

    hdr("§4 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_nearest"] = top
    RES["_span_log10"] = span
    RES["_target"] = {"mu": ref["_target"]["mu"], "log10_mu": ref["_target"]["log10_mu"]}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL}
    (DATA / "13_nearest_diag.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_13.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / '13_nearest_diag.json'}\n寫出: {DATA / 'run_13.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
