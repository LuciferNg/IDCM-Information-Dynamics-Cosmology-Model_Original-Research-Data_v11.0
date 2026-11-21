#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L2 — M6 獨立重推導（idcm-audit-verify §2.4）
================================================================================
由**原始輸入**重算本記錄全部 headline 數字；唔 import L1、唔讀 docs。
產出 data/l2_independent_verify.json ＋ data/run_l2.log（供 L3）。
"""
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent if HERE.name == "scripts" else HERE
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

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
    try:
        xf = float(x)
    except (OverflowError, ValueError):
        xf = None
    if xf is not None and (xf == 0 or 1e-300 < abs(xf) < 1e300):
        return {"float": xf, "str": mp.nstr(x, 25)}
    return {"float": None, "str": mp.nstr(x, 25)}


# ① INPUTS —— 原始量（計數／測量／上游已封結論；零推導結果）
INPUTS = {
    "N": 88, "N_d": 44, "B": 4,
    "M0_MeV": "27.900915",
    "hbar_c_MeV_fm": "197.3269804",
    "M_Pl_MeV": "1.220890e22",
    "m_p_MeV": "938.27208816",
    "v_SM_MeV": "246219.6",
    "a0_v1_fm": "0.377422", "a0_v2_fm": "0.384184",
    "t_v1_MeV": "0.7884", "t_v2_MeV": "1.7711",
    "kappa_lo_MeV": "68.5696", "kappa_hi_MeV": "70.5744",
    "H0_pred": "67.4071", "H0_planck": "67.36", "H0_sigma": "0.54",
}


def gam(n, I):
    return 4 * mp.sin(mp.mpf(n) * mp.pi / I["N"]) ** 2


def g_kernel(d, I):
    N = I["N"]
    return (mp.mpf(N * N - 1) - 6 * d * (N - d)) / (12 * N)


# ② DERIVE —— 純函數
DERIVE = {
    "sum_gamma":  lambda I: mp.fsum(gam(n, I) for n in range(I["N_d"] + 1)),
    "pair_dev":   lambda I: mp.fsum(abs(gam(n, I) + gam(I["N_d"] - n, I) - 4) for n in range(I["N_d"] + 1)),
    "lam_max_1":  lambda I: 4 * 1 - gam(0, I),
    "n_strict":   lambda I: mp.mpf(sum(1 for k in range(2, I["N_d"] + 1)
                                      if (4 * k - mp.fsum(gam(j, I) for j in range(k))) < 4 * k)),
    "zero_root":  lambda I: mp.mpf(I["N_d"]) - mp.sqrt(mp.mpf(2582)) / 2,
    "g18":        lambda I: g_kernel(18, I),
    "g19":        lambda I: g_kernel(19, I),
    "capacity_MeV": lambda I: mp.fsum(gam(n, I) for n in range(I["N_d"] + 1)) * mp.mpf(I["M0_MeV"]),
    "nq_mp":      lambda I: mp.ceil(mp.mpf(I["m_p_MeV"]) / mp.mpf(I["M0_MeV"])
                                    / mp.fsum(gam(n, I) for n in range(I["N_d"] + 1))),
    "nq_vsm":     lambda I: mp.ceil(mp.mpf(I["v_SM_MeV"]) / mp.mpf(I["M0_MeV"])
                                    / mp.fsum(gam(n, I) for n in range(I["N_d"] + 1))),
    "mpl_over_m0": lambda I: mp.mpf(I["M_Pl_MeV"]) / mp.mpf(I["M0_MeV"]),
    "nq_mpl":     lambda I: mp.mpf(I["M_Pl_MeV"]) / mp.mpf(I["M0_MeV"])
                             / mp.fsum(gam(n, I) for n in range(I["N_d"] + 1)),
    "TA_v1":      lambda I: mp.mpf(I["hbar_c_MeV_fm"]) / (mp.mpf(I["a0_v1_fm"]) * mp.mpf(I["M0_MeV"])),
    "TA_v2":      lambda I: mp.mpf(I["hbar_c_MeV_fm"]) / (mp.mpf(I["a0_v2_fm"]) * mp.mpf(I["M0_MeV"])),
    "Tt_v1":      lambda I: mp.mpf(I["t_v1_MeV"]) / mp.mpf(I["M0_MeV"]),
    "Tt_v2":      lambda I: mp.mpf(I["t_v2_MeV"]) / mp.mpf(I["M0_MeV"]),
    "kappa_lo":   lambda I: mp.mpf(I["kappa_lo_MeV"]) / mp.mpf(I["M0_MeV"]),
    "kappa_hi":   lambda I: mp.mpf(I["kappa_hi_MeV"]) / mp.mpf(I["M0_MeV"]),
    "H0_z":       lambda I: abs(mp.mpf(I["H0_pred"]) - mp.mpf(I["H0_planck"])) / mp.mpf(I["H0_sigma"]),
}

# ③ IDENTITIES —— 兩條獨立路徑
IDENTITIES = [
    ("I1 Σ_{0..44}γ = N+2（數值 vs 閉式）",
     lambda I: mp.fsum(gam(n, I) for n in range(I["N_d"] + 1)),
     lambda I: mp.mpf(I["N"] + 2)),
    ("I2 γ_n + γ_(44−n) = 4 全 n（最大偏差 = 0）",
     lambda I: mp.fsum(abs(gam(n, I) + gam(I["N_d"] - n, I) - 4) for n in range(I["N_d"] + 1)),
     lambda I: mp.mpf(0)),
    ("I3 零交叉：g(d*) = 0（閉式根代入）",
     lambda I: g_kernel(mp.mpf(I["N_d"]) - mp.sqrt(mp.mpf(2582)) / 2, I),
     lambda I: mp.mpf(0)),
    ("I4 λ_max(1) = γ₄₄ = 4k（兩路徑）",
     lambda I: 4 * 1 - gam(0, I),
     lambda I: gam(I["N_d"], I)),
    ("I5 每量子容量 = Σγ·M₀（兩路徑）",
     lambda I: mp.fsum(gam(n, I) for n in range(I["N_d"] + 1)) * mp.mpf(I["M0_MeV"]),
     lambda I: mp.mpf(I["N"] + 2) * mp.mpf(I["M0_MeV"])),
]


def main():
    hdr("§0 輸入（原始量）")
    I = dict(INPUTS)
    for k, v in I.items():
        print(f"  {k:16s} = {v}")

    hdr("§1 DERIVE")
    D = {}
    for key, fn in DERIVE.items():
        try:
            D[key] = mp.mpf(fn(I))
            print(f"  {key:14s} = {mp.nstr(D[key], 20)}")
        except Exception as e:  # noqa: BLE001
            chk(f"D.{key} 可推導", False, f"{type(e).__name__}: {e}")
    chk(f"D 全部 key 可推導（{len(DERIVE)} 個）", len(D) == len(DERIVE), f"{len(D)}/{len(DERIVE)}")

    hdr("§2 IDENTITIES（差 < 1e−30）")
    for tag, fa, fb in IDENTITIES:
        a, b = mp.mpf(fa(I)), mp.mpf(fb(I))
        chk(tag, abs(a - b) < mp.mpf("1e-30"), f"差 = {mp.nstr(abs(a-b), 6)}")

    hdr("§3 內在關係")
    chk("C1 Σγ 為整數 90", D.get("sum_gamma") == 90, f"= {mp.nstr(D.get('sum_gamma', mp.mpf(0)), 6)}")
    chk("C2 每量子容量 = 90·M₀ = 2511.08 MeV",
        abs(D.get("capacity_MeV", mp.mpf(0)) - mp.mpf("2511.082350")) < mp.mpf("1e-5"),
        mp.nstr(D.get("capacity_MeV", mp.mpf(0)), 12))
    chk("C3 n_q(m_p) = 1 ∧ n_q(v_SM) = 99",
        D.get("nq_mp") == 1 and D.get("nq_vsm") == 99,
        f"{mp.nstr(D.get('nq_mp', 0), 4)} / {mp.nstr(D.get('nq_vsm', 0), 4)}")
    chk("C4 零交叉落 (18,19)，g(18)>0>g(19)",
        mp.mpf(18) < D.get("zero_root", 0) < mp.mpf(19) and D.get("g18", 0) > 0 > D.get("g19", 0))
    chk("C5 κ 窗 > 1 ⟹ κ 屬 M₀ 標度", D.get("kappa_lo", 0) > 1 and D.get("kappa_hi", 0) > 1,
        f"[{mp.nstr(D.get('kappa_lo',0),8)}, {mp.nstr(D.get('kappa_hi',0),8)}] M₀")

    hdr("§4 匯總")
    nchk = len([k for k in RES])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    _pas = nchk - NFAIL
    print(f"{_pas}/{nchk} PASS（L2）")
    LOG.append(f"{_pas}/{nchk} PASS（L2）")
    RES["_inputs"] = {k: (v if not isinstance(v, mp.mpf) else mp.nstr(v, 30)) for k, v in I.items()}
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}
    RES["pass"] = nchk - NFAIL
    RES["fail"] = NFAIL
    (DATA / "l2_independent_verify.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l2.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: data/l2_independent_verify.json\ndata/run_l2.log")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())