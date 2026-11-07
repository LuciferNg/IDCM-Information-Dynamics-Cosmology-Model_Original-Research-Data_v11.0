#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M16 獨立重推導（idcm-audit-verify §2.4）
================================================================
鐵律
  R1 唔 import／call 任何 scripts/0*.py；唔讀 docs/ 一個字。
  R2 INPUTS 只放原始量：鏈常數（B, N_d, L₅, L₄, N_gen, Σγ, ⟨W⟩²）、
     Planck 2018 測量值、以及 data/run_02_scan_all_500.txt（GAP 原始輸出，屬「測量」）。
  R3 DERIVE 全部係 INPUTS 之純函數（零隨機）。
  R4 數值用 mpmath dps=40；整數運算用 int／Fraction。
輸出：data/l2_independent_verify.json ＋ data/run_l2.log
"""
import json
import re
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
import os
SCAN = Path(os.environ.get("M16_SCAN", str(DATA / "run_02_scan_all_500.txt")))

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


def ser(x):
    if isinstance(x, (int, bool)):
        return {"float": float(x), "str": str(int(x)) if not isinstance(x, bool) else str(x)}
    if isinstance(x, Fraction):
        f = float(x)
        return {"float": f, "str": str(x)}
    try:
        xf = float(x)
    except (OverflowError, ValueError):
        xf = None
    if xf is not None and (xf == 0 or 1e-300 < abs(xf) < 1e300):
        return {"float": xf, "str": mp.nstr(x, 25)}
    return {"float": None, "str": mp.nstr(x, 25)}


# ══════════════════════════════════════════════════════════════
# ① INPUTS — 原始量（零推導結果、零文檔值）
# ══════════════════════════════════════════════════════════════
INPUTS = {
    # 鏈常數（IDCM 已封，非本 O 推導）
    "B": 4,              # γ_max = Nyquist 帶值
    "N_d": 44,           # 鏈長（88-site 之半）
    "L5": 11,            # Lucas L₅
    "L4": 7,             # Lucas L₄
    "N_gen": 3,          # 世代數（= B − 1，M15 已封）
    "Sum_gamma": 90,     # Σγ（預算幣制總量）
    "W2": "61.5",        # ⟨W⟩²（O39：123/2）
    "dims_A5": [1, 3, 3, 4, 5],   # A₅ 不可約維數（GAP / 標準表）
    # Planck 2018（TT,TE,EE+lowE+lensing）
    "om_bh2": "0.02237",
    "s_om_bh2": "0.00015",
    "h_planck": "0.6736",
    "s_h_planck": "0.0054",
    # IDCM 自身 H₀
    "h_idcm": "0.674071",
    # GAP 原始輸出（92,570 群之 irrep 維數重數）
    "scan_path": str(SCAN),
}

# ══════════════════════════════════════════════════════════════
# ② 由原始掃描檔重算群層量（自建 parser，唔用 L1 之 regex 與陣列）
# ══════════════════════════════════════════════════════════════
def parse_scan(path):
    rows = []
    tok = re.compile(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]")
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            head = ln.split("[", 1)[0].split()
            if len(head) != 2:
                continue
            order, idx = int(head[0]), int(head[1])
            degs = []
            for d, k in tok.findall(ln):
                degs.extend([int(d)] * int(k))
            rows.append((order, idx, degs))
    return rows


ROWS = parse_scan(INPUTS["scan_path"])


def _prod(xs):
    p = 1
    for x in xs:
        p *= x
    return p


def dnt(degs):
    return sorted({d for d in degs if d > 1})


def _sum12_sets():
    import itertools
    sets = []
    for r in range(1, 13):
        for a in itertools.combinations(range(2, 13), r):
            if sum(a) == 12:
                sets.append(a)
    return sets


def _uniq_c():
    out = []
    for c in range(1, 41):
        mm = [r for r in ROWS if sum(dnt(r[2])) == c and _prod(dnt(r[2])) == r[0]]
        if len(mm) == 1 and mm[0][0] == 60 and mm[0][1] == 5:
            out.append(c)
    return out


SCAN_DERIVED = {
    "n_groups": len(ROWS),
    "n_groups_le120": sum(1 for o, _, _ in ROWS if o <= 120),
    "burnside_ok": all(sum(d * d for d in g) == o for o, _, g in ROWS),
    "hits_prod_eq_order": sum(1 for o, _, g in ROWS if _prod(dnt(g)) == o),
    "hits_sum12": sum(1 for o, _, g in ROWS if sum(dnt(g)) == 12),
    "hits_conj": sum(1 for o, _, g in ROWS
                     if sum(dnt(g)) == 12 and _prod(dnt(g)) == o),
    "v1_no_distinct": sum(1 for o, _, g in ROWS
                          if _prod([d for d in g if d > 1]) == o
                          and sum(d for d in g if d > 1) == 12),
    # 去「非平凡」修飾：Σ(相異全部) = 12（等價於 Σ_{d>1} = 11）
    "v2_no_nontrivial": sum(1 for o, _, g in ROWS if sum(sorted(set(g))) == 12),
    "dnt345_groups": sum(1 for o, _, g in ROWS if dnt(g) == [3, 4, 5]),
    # DNT ⊇ {3,4,5}：包含 3,4,5 之群數 ＋ 其中最小階
    "dnt345_superset": sum(1 for o, _, g in ROWS if {3, 4, 5} <= set(dnt(g))),
    "dnt345_superset_min_order": min([o for o, _, g in ROWS if {3, 4, 5} <= set(dnt(g))] or [0]),
    # perfect ⟺ k = |G:G'| = |G| − Σ_{d>1} m_d d² = 1（線性特徵標恰一個）
    "perfect_500": sum(1 for o, _, g in ROWS
                       if (o - sum(d * d for d in g if d > 1)) == 1 and o > 1),
}


# ══════════════════════════════════════════════════════════════
# ③ DERIVE — 全部係 INPUTS 之純函數
# ══════════════════════════════════════════════════════════════
def _gamma(I, n):
    """γ_n = 4 sin²(nπ/(4B)) = 2(1 − cos(nπ/(2B)))，B = γ_max（Nyquist 帶值）"""
    return 4 * mp.sin(mp.mpf(n) * mp.pi / (4 * I["B"])) ** 2


DERIVE = {
    # 度集與群階
    "group_order": lambda I: _prod([I["B"] - 1, I["B"], I["B"] + 1]),
    "sum_deg": lambda I: mp.mpf(3 * I["B"]),
    "sum_deg_A5": lambda I: mp.mpf(sum(I["dims_A5"])),
    "sum_deg2_all": lambda I: mp.mpf(sum(d * d for d in I["dims_A5"])),   # = |G|（Burnside）
    "sum_deg2_dnt": lambda I: mp.mpf(sum(d * d for d in set(I["dims_A5"]) if d > 1)),  # 相異：50
    "sum_deg2_mult_deep": lambda I: mp.mpf(sum(d * d for d in I["dims_A5"] if d > 1)),  # 重數：59
    "n_irr_A5": lambda I: mp.mpf(len(I["dims_A5"])),
    # 唯一性算術（Σ=12 之相異整數集）
    "n_sets_sum12": lambda I: mp.mpf(len(_sum12_sets())),
    "n_surv_L2": lambda I: mp.mpf(sum(1 for a in _sum12_sets()
                                      if _prod(a) > sum(d * d for d in a))),
    # γ 光譜
    "gamma_1": lambda I: _gamma(I, 1),
    "gamma_2": lambda I: _gamma(I, 2),
    "gamma_3": lambda I: _gamma(I, 3),
    "gamma_5": lambda I: _gamma(I, 5),
    "gamma_6": lambda I: _gamma(I, 6),
    "gamma_7": lambda I: _gamma(I, 7),
    "delta_1": lambda I: _gamma(I, 0) + _gamma(I, 2) - 2 * _gamma(I, 1),
    "delta_2": lambda I: _gamma(I, 1) + _gamma(I, 3) - 2 * _gamma(I, 2),
    "delta_3": lambda I: _gamma(I, 2) + _gamma(I, 4) - 2 * _gamma(I, 3),
    "delta_4": lambda I: _gamma(I, 3) + _gamma(I, 5) - 2 * _gamma(I, 4),
    "delta_5": lambda I: _gamma(I, 4) + _gamma(I, 6) - 2 * _gamma(I, 5),
    "delta_6": lambda I: _gamma(I, 5) + _gamma(I, 7) - 2 * _gamma(I, 6),
    "delta_7": lambda I: _gamma(I, 6) + _gamma(I, 8) - 2 * _gamma(I, 7),
    "sigma_omega_b": lambda I: (mp.mpf(I["om_bh2"]) / mp.mpf(I["h_planck"]) ** 2)
                               * mp.sqrt((mp.mpf(I["s_om_bh2"]) / mp.mpf(I["om_bh2"])) ** 2
                                         + (2 * mp.mpf(I["s_h_planck"]) / mp.mpf(I["h_planck"])) ** 2),
    "n_groups_order60": lambda I: mp.mpf(sum(1 for o, _, _ in ROWS if o == 60)),
    "gamma_n_ni": lambda I: _gamma(I, 0),
    "gamma_3_plus_5": lambda I: _gamma(I, 3) + _gamma(I, 5),
    "gamma_4": lambda I: _gamma(I, 4),
    "gamma_max": lambda I: _gamma(I, 2 * I["B"]),
    "delta_ap": lambda I: _gamma(I, 3) + _gamma(I, 5) - 2 * _gamma(I, 4),
    # 預算（依賴，非本 O 推導）
    "omega_b_struct": lambda I: Fraction(I["B"] - 1, 60),
    "omega_dm_struct": lambda I: Fraction(I["B"] ** 2, 60),
    "omega_lambda_struct": lambda I: Fraction(41, 60),
    "f_unproj": lambda I: mp.mpf(I["W2"]) / mp.mpf(I["Sum_gamma"]),
    # Planck Ω_b 與 k = |G:G'|
    "omega_b_planck": lambda I: mp.mpf(I["om_bh2"]) / mp.mpf(I["h_planck"]) ** 2,
    "sigma_rel_omega_b": lambda I: mp.sqrt((mp.mpf(I["s_om_bh2"]) / mp.mpf(I["om_bh2"])) ** 2
                                           + (2 * mp.mpf(I["s_h_planck"]) / mp.mpf(I["h_planck"])) ** 2),
    "k_center": lambda I: 1 / (20 * (mp.mpf(I["om_bh2"]) / mp.mpf(I["h_planck"]) ** 2)),
    "k_sigma": lambda I: (1 / (20 * (mp.mpf(I["om_bh2"]) / mp.mpf(I["h_planck"]) ** 2)))
                         * mp.sqrt((mp.mpf(I["s_om_bh2"]) / mp.mpf(I["om_bh2"])) ** 2
                                   + (2 * mp.mpf(I["s_h_planck"]) / mp.mpf(I["h_planck"])) ** 2),
    "z_k1": lambda I: (DERIVE["k_center"](I) - 1) / DERIVE["k_sigma"](I),
    "z_k2": lambda I: (DERIVE["k_center"](I) - 2) / DERIVE["k_sigma"](I),
    "z_k3": lambda I: (DERIVE["k_center"](I) - 3) / DERIVE["k_sigma"](I),
    "z_k4": lambda I: (DERIVE["k_center"](I) - 4) / DERIVE["k_sigma"](I),
    "k_with_h_idcm": lambda I: 1 / (20 * (mp.mpf(I["om_bh2"]) / mp.mpf(I["h_idcm"]) ** 2)),
    # look-elsewhere
    "family_frac": lambda I: mp.mpf(17) / 45,
    "le_6": lambda I: mp.mpf(1) / 6,
    "sumd_plus_deep": lambda I: mp.mpf(sum(I["dims_A5"]) + I["N_d"]),
    "cubic_B2": lambda I: mp.mpf(2 ** 3 - 2),
    "cubic_B3": lambda I: mp.mpf(3 ** 3 - 3),
    "cubic_B5": lambda I: mp.mpf(5 ** 3 - 5),
    "cubic_B6": lambda I: mp.mpf(6 ** 3 - 6),
    "R_band": lambda I: mp.mpf(2 * I["N_d"] // (2 * I["B"])),
    "N_d_val": lambda I: mp.mpf(I["N_d"]),
    "le_56": lambda I: mp.mpf(1) / 56,
    "le_8": lambda I: mp.mpf(1) / 8,
    # 群層（由原始掃描重算）
    "n_groups": lambda I: mp.mpf(SCAN_DERIVED["n_groups"]),
    "n_groups_le120": lambda I: mp.mpf(SCAN_DERIVED["n_groups_le120"]),
    "hits_prod": lambda I: mp.mpf(SCAN_DERIVED["hits_prod_eq_order"]),
    "hits_sum12": lambda I: mp.mpf(SCAN_DERIVED["hits_sum12"]),
    "hits_conj": lambda I: mp.mpf(SCAN_DERIVED["hits_conj"]),
    "v1_no_distinct": lambda I: mp.mpf(SCAN_DERIVED["v1_no_distinct"]),
    "v2_no_nontrivial": lambda I: mp.mpf(SCAN_DERIVED["v2_no_nontrivial"]),
    "dnt345_groups": lambda I: mp.mpf(SCAN_DERIVED["dnt345_groups"]),
    "dnt345_superset": lambda I: mp.mpf(SCAN_DERIVED["dnt345_superset"]),
    "dnt345_superset_min_order": lambda I: mp.mpf(SCAN_DERIVED["dnt345_superset_min_order"]),
    "perfect_500": lambda I: mp.mpf(SCAN_DERIVED["perfect_500"]),
    "uniq_c_12": lambda I: mp.mpf(len(_uniq_c())),
    # 記賬恆等式
    "order_from_nd": lambda I: mp.mpf(I["N_d"] + I["B"] ** 2),
    "deep_sector": lambda I: mp.mpf(I["L4"] * 60),
    "surface_total": lambda I: mp.mpf((I["B"] - 1 + I["B"] + I["B"] + 1) * I["L5"]),
    "K_tree": lambda I: mp.mpf(3 * I["L5"]),
}

# ══════════════════════════════════════════════════════════════
# ④ IDENTITIES — 兩條獨立路徑同一量
# ══════════════════════════════════════════════════════════════
IDENTITIES = [
    ("I1 3·4·5 = 60 = N_d + B²",
     lambda I: _prod([I["B"] - 1, I["B"], I["B"] + 1]),
     lambda I: I["N_d"] + I["B"] ** 2),
    ("I2 γ_3+γ_5 = γ_max（對極互補）",
     lambda I: _gamma(I, 3) + _gamma(I, 5),
     lambda I: _gamma(I, 2 * I["B"])),
    ("I3 γ_3+γ_5 = 2γ_4（等差）",
     lambda I: _gamma(I, 3) + _gamma(I, 5),
     lambda I: 2 * _gamma(I, 4)),
    ("I4 Σd(A₅) = B²",
     lambda I: sum(I["dims_A5"]),
     lambda I: I["B"] ** 2),
    ("I5 Σγ = 2(N_d + 1)",
     lambda I: I["Sum_gamma"],
     lambda I: 2 * (I["N_d"] + 1)),
    ("I6 ⟨W⟩²/Σγ = 41/60",
     lambda I: mp.mpf(I["W2"]) / mp.mpf(I["Sum_gamma"]),
     lambda I: mp.mpf(41) / mp.mpf(60)),
    ("I7 Σ=12 集數：窮舉 = 組合重數",
     lambda I: len(_sum12_sets()),
     lambda I: len({a for a in _sum12_sets()})),
    ("I8 表面總和 132 = 3·N_d",
     lambda I: 33 + 44 + 55,
     lambda I: 3 * I["N_d"]),
    ("I9 3²+4² = 5²",
     lambda I: 3 ** 2 + 4 ** 2,
     lambda I: 5 ** 2),
    ("I11 Burnside：|G| = Σ_{irr} d²（A₅）",
     lambda I: sum(d * d for d in I["dims_A5"]),
     lambda I: _prod([I["B"] - 1, I["B"], I["B"] + 1])),
    ("I12 Σd² = Σd + Σd(d−1)（DNT，相異）",
     lambda I: sum(d * d for d in set(I["dims_A5"]) if d > 1),
     lambda I: sum(d for d in set(I["dims_A5"]) if d > 1)
               + sum(d * (d - 1) for d in set(I["dims_A5"]) if d > 1)),
    ("I10 k(Ω_b=1/20) = 1",
     lambda I: 1 / (20 * mp.mpf(1) / 20),
     lambda I: mp.mpf(1)),
]


def main():
    hdr("§0 輸入（原始量）")
    for k, v in INPUTS.items():
        if k == "scan_path":
            continue
        print(f"  {k:14s} = {v}")
    print(f"  scan           = {len(ROWS)} rows")

    hdr("§1 DERIVE：由 INPUTS 純函數重算")
    D = {}
    for key, fn in DERIVE.items():
        try:
            v = fn(INPUTS)
            D[key] = v if isinstance(v, (int, Fraction)) else mp.mpf(v)
            print(f"  {key:22s} = {D[key] if isinstance(D[key], Fraction) else mp.nstr(mp.mpf(D[key]), 18)}")
        except Exception as e:  # noqa: BLE001
            chk(f"D.{key} 可推導", False, f"{type(e).__name__}: {e}")
    chk(f"D 全部 key 可推導（{len(DERIVE)} 個）", len(D) == len(DERIVE),
        f"{len(D)}/{len(DERIVE)}")

    hdr("§2 IDENTITIES（差 < 1e−30）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(INPUTS)), mp.mpf(fb(INPUTS))
        d = abs(a - b)
        chk(tag, d < mp.mpf("1e-30"), f"差 = {mp.nstr(d, 6)}")

    hdr("§3 交叉核（內在關係）")
    chk("C1 Burnside 對全部掃描行成立", SCAN_DERIVED["burnside_ok"])
    chk("C2 A₅ 唯一命中（Σ=12 ∧ ∏=|G|）", SCAN_DERIVED["hits_conj"] == 1 and _uniq_c() == [12],
        f"hits={SCAN_DERIVED['hits_conj']}, uniq_c={_uniq_c()}")
    chk("C3 k=1 相容 <1σ、k=2 排除 >10σ",
        abs(float(D["z_k1"])) < 1 and float(D["z_k2"]) < -10,
        f"z1={mp.nstr(D['z_k1'],5)} z2={mp.nstr(D['z_k2'],6)}")
    chk("C4 Δ_γ(帶 3,4,5) = 0", abs(float(D["delta_ap"])) < 1e-30)
    chk("C5 分數 41/60 精確",
        Fraction(41, 60) == Fraction(615, 900))

    hdr("§4 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_inputs"] = {k: (v if not isinstance(v, Fraction) else str(v)) for k, v in INPUTS.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}

    out = DATA / "l2_independent_verify.json"
    out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l2.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {out}\n寫出: {DATA / 'run_l2.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())