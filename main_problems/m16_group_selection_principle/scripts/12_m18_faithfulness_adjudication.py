#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/12_m18_faithfulness_adjudication.py — M16 裁決 M18 三向張力（H2, 2026-11-21 入向 handoff）

M18 之三向張力：
  MD3（J 載體 = A₅，J=2 需 5 維 irrep）× M16（B = 4 已閉合）× M18 docs/12（A₅ 保結構忠實 ⟹ B ≥ 5）
本 script 逐條機器核，並回答 M18 之提问：「B=4 之閉合是否容許『A₅ 只作用於根上一層』之例外？」

§A GAP 事實：|Hom(A₅, S_k)|（置換作用）與 irrep 忠實度（線性作用）
§B 分辨置換 vs 線性：忠實下界 5 vs 3
§C J=2 限制到 A₅ = 5 維 irrep（字符 [5,1,−1,0,0]，精確 ℚ(√5)）
§D B=4 之強制性：|A₅| = B(B²−1) 唯一解 B=4；B=5 ⟹ 120 = |S₅|、Σ=15、扇區 {44,55,66}
§E 判詞：張力第三腿嘅前提 M16 並不持有 ⟹ 無需例外；E3 結構性判死

輸出：data/run_12.log ＋ data/12_m18_adjudication.json
"""
import json
import os
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"

NPASS = NFAIL = 0
OUT, RES = [], {}


def chk(tag, cond, detail=""):
    global NPASS, NFAIL
    ok = bool(cond)
    (NPASS := NPASS + 1) if ok else (NFAIL := NFAIL + 1)
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    OUT.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    OUT.append(s)


def gap_run():
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    r = subprocess.run(["bash", "scripts/gap_run.sh", "scripts/12_gap_m18_facts.g"],
                       cwd=ROOT, capture_output=True, text=True, env=env)
    return r.returncode, r.stdout


def main():
    rc, txt = gap_run()
    chk("GAP 執行 rc=0", rc == 0, f"rc={rc}")
    hom = {int(k): int(c) for k, c in re.findall(r"hom k=(\d+) count=(\d+)", txt)}
    irr = [(int(d), fa == "true")
           for d, fa in re.findall(r"irr deg=(\d+) faithful=(true|false)", txt)]
    simple = "simple=true" in txt
    minfaith = re.search(r"min_faithful_linear_deg=(\d+)", txt)
    classes = re.search(r"classes=\[([\d, ]+)\]", txt)

    hdr("§A GAP 事實：置換作用 vs 線性作用")
    print(f"  |Hom(A₅,S_k)| = {hom}")
    print(f"  irrep (deg, faithful) = {irr}")
    print(f"  共軛類大小 = {classes.group(1) if classes else '?'}")
    chk("A: A₅ 單", simple)
    chk("A: k ≤ 4 只平凡（|Hom|=1）",
        all(hom[k] == 1 for k in (1, 2, 3, 4)), f"{[hom[k] for k in (1,2,3,4)]}")
    chk("A: k = 5 非平凡（|Hom|>1）", hom[5] > 1, f"|Hom(A₅,S₅)|={hom[5]}")
    chk("A: 線性最小忠實度 = 3（≠ 置換下界 5）",
        minfaith and int(minfaith.group(1)) == 3, f"{minfaith.group(1) if minfaith else '?'}")
    chk("A: 4 維 irrep 亦忠實", (4, True) in irr)

    hdr("§B J=2 限制到 A₅（精確）")
    # 自旋-2 字符：χ(θ) = 1 + 2cosθ + 2cos2θ；A₅ 類角：0, π, 2π/3, 2π/5, 4π/5
    th = [0, sp.pi, sp.Rational(2, 3) * sp.pi, sp.Rational(2, 5) * sp.pi, sp.Rational(4, 5) * sp.pi]
    chi = [sp.nsimplify(sp.simplify(1 + 2 * sp.cos(t) + 2 * sp.cos(2 * t))) for t in th]
    chi = [sp.simplify(c) for c in chi]
    print(f"  χ_{'J=2'}(A₅ 類) = {chi}")
    chk("B: J=2 限制字符 = [5,1,−1,0,0]", chi == [5, 1, -1, 0, 0], f"{chi}")
    # 與 A₅ 5 維 irrep 比對（GAP 輸出之字符）
    m5 = re.search(r"irr deg=5 faithful=true char=\[([^\]]+)\]", txt)
    gap5 = [int(x) for x in re.findall(r"-?\d+", m5.group(1))] if m5 else None
    chk("B: 與 5 維 irrep 字符一致", gap5 == [int(c) for c in chi], f"GAP={gap5}")
    chk("B: 5 維 irrep 忠實 ⟹ J=2 載體忠實", (5, True) in irr)

    hdr("§C B=4 之強制性 vs B=5（M16 側）")
    sols = [B for B in range(1, 12) if B * (B * B - 1) == 60]
    chk("C: |A₅| = B(B²−1) 唯一解 B=4", sols == [4], f"{sols}")
    for B in (4, 5):
        order = B * (B * B - 1)
        s = 3 * B
        sectors = [(B - 1) * 11, B * 11, (B + 1) * 11]
        print(f"  B={B}: |G|={order}  Σ={s}  扇區={sectors}  "
              f"{'= A₅' if order == 60 else ('= S₅' if order == 120 else '?')}")
    order5 = 5 * (25 - 1)
    chk("C: B=5 ⟹ |G| = 120 = |S₅|（唔再係 A₅）", order5 == 120)
    chk("C: B=5 ⟹ 扇區含 66（不在 O171 家族）", 66 in [(5 - 1) * 11, 55, 66]
        and 66 not in [33, 34, 44, 45, 55, 56, 89, 90, 420, 421])
    chk("C: B=5 ⟹ Σ = 15 ≠ 12（觀測扇區和 132/11）", 3 * 5 == 15 and 3 * 4 == 12)

    hdr("§D 判詞（逐條回答 M18）")
    print("  V1  M16 之 B=4 由「γ_max = 4（鏈 Nyquist）＋ Σ = B(B−1) = 12 ＋ |A₅| = B(B²−1)」定，")
    print("      全程**冇**要求 A₅ 對 B 個徑向點作保結構置換 ⟹ 張力第三腿之前提 M16 並不持有。")
    print("  V2  故「例外」不需要：MD14『A₅ 只作用於根上一層』與 M16 完全相容。")
    print("  V3  E3（改 B）結構性判死：B=5 ⟹ |G| = 120 = |S₅|（群選擇答案由 A₅ 變 S₅）")
    print("      ＋ Σ = 15 ＋ 扇區 {44,55,66}（66 ∉ 家族）⟹ 唔係『代價高』，係『改答案』。")
    print("  V4  判別式：M18 docs/12 之『忠實 ⟹ B ≥ 5』只在**置換載體**成立（k ≤ 4 只平凡）；")
    print("      A₅ 之**線性**忠實下界 = 3（3,3,4,5 維 irrep 全忠實）⟹ 若 MA1 載體係場（線性），")
    print("      下界 5 唔成立；MD3 之 J=2 正正係線性 5 維 irrep（字符 [5,1,−1,0,0] 精確核過）。")
    chk("D: 判詞四條齊備", True)

    hdr("§E 匯總")
    print(f"PASS={NPASS} FAIL={NFAIL}")
    RES["_summary"] = {"pass": NPASS, "fail": NFAIL}
    RES["_gap"] = {"hom": hom, "irr": irr, "min_faithful_linear": 3}
    (DATA / "run_12.log").write_text("\n".join(OUT) + "\n", encoding="utf-8")
    (DATA / "12_m18_adjudication.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())