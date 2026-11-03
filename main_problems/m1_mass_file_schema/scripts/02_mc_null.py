#!/usr/bin/env python3
"""M1 scripts/02 — R2 補充：MC null（10⁷）與判別力位元數

問題：觀測到「一個編碼」或「ρ = 61」算唔算證據？用隨機目標（均勻／對數均勻）作 null。

方法：
  - 以 k ≤ 4 全部子集和（149,986 個）為池，對 10⁷ 個隨機目標計 ρ（ε=1e-3）
  - null A：m ~ U(0,90)（遍歷整個帳本）
  - null B：m ~ 對數均勻於 [m_e/M₀, m_t/M₀]（模擬「觀測質量都係對數分佈」的選擇效應）
  - 對每個費米子報 p = P(ρ_null ≥ ρ_obs) 及 p0 = P(ρ_null = 0)
  - 證據位元 = −log₂ p（上限 1 bit 之說法以 p0 檢查）
輸出：data/02_mc_null.json
"""
import itertools
import json
import os

import numpy as np
from mpmath import mp, mpf, cos, pi

mp.dps = 80
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
os.makedirs(DATA, exist_ok=True)

N = 44
M0 = mpf("27.9009145112")
EPS = 1e-3
KMAX = 4
DRAWS = 10_000_000

GAMMA = np.array([float(2 - 2 * cos(mpf(n) * pi / 44)) for n in range(1, N + 1)])
FERMIONS = {"e": "0.51099895000", "mu": "105.6583755", "tau": "1776.86",
            "u": "2.16", "d": "4.67", "s": "93.4", "c": "1273.0", "b": "4180.0", "t": "172570.0"}
OBS = {"e": 0, "mu": 61, "tau": 0, "u": 0, "d": 0, "s": 40, "c": 0, "b": 0, "t": 0}  # 由 scripts/01 得（k≤4, 1e-3）


def main():
    sums = np.sort(np.concatenate([
        np.array([float(sum(GAMMA[i] for i in c))
                  for c in itertools.combinations(range(N), kk)], dtype=float) if kk else np.zeros(1)
        for kk in range(0, KMAX + 1)]))

    def rho_of(m):
        lo = np.searchsorted(sums, m * (1 - EPS), side="left")
        hi = np.searchsorted(sums, m * (1 + EPS), side="right")
        return hi - lo

    rng = np.random.default_rng(20261028)
    draws_A = rng.uniform(0.0, 90.0, DRAWS)
    lo_m, hi_m = float(mpf("0.51099895000") / M0), float(mpf("172570.0") / M0)
    draws_B = np.exp(rng.uniform(np.log(lo_m), np.log(hi_m), DRAWS))

    out = {"meta": {"script": "02_mc_null", "draws": DRAWS, "eps": EPS, "kmax": KMAX,
                    "pool_size": int(sums.size), "seeds": 20261028,
                    "nullA": "U(0,90)", "nullB": f"log-uniform[{lo_m:.6g},{hi_m:.6g}]"},
           "obs": OBS}
    for tag, draws in (("A", draws_A), ("B", draws_B)):
        r = rho_of(draws)
        hist = {str(int(v)): int((r == v).sum()) for v in (0, 1, 2, 3, 5, 10, 20, 40, 60)}
        zeros = float((r == 0).mean())
        rows = {}
        for name, obs in OBS.items():
            p_ge = float((r >= obs).mean())
            rows[name] = {"rho_obs": obs, "p_ge": p_ge,
                          "bits": (float(np.log2(1 / p_ge)) if p_ge > 0 else None)}
        out[f"null_{tag}"] = {"mean_rho": float(r.mean()), "median_rho": float(np.median(r)),
                              "p_rho_zero": zeros, "hist": hist, "per_fermion": rows}

    checks = [("nullA_p_zero_between_0_1", 0.0 <= out["null_A"]["p_rho_zero"] <= 1.0),
              ("mean_rho_gt_0", out["null_A"]["mean_rho"] > 0),
              ("pool_size", sums.size == 149986)]
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["check_summary"] = {"total": len(checks), "pass": sum(1 for _, p in checks if p)}

    with open(os.path.join(DATA, "02_mc_null.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print(json.dumps({"checks": out["check_summary"], "pool": int(sums.size)}, ensure_ascii=False))
    for tag in ("A", "B"):
        d = out[f"null_{tag}"]
        print(f"[null {tag}] mean_rho={d['mean_rho']:.2f} median={d['median_rho']:.0f} "
              f"P(rho=0)={d['p_rho_zero']:.4f} hist={d['hist']}")
        print("   " + "  ".join(f"{k}:p={v['p_ge']:.3g}" for k, v in d["per_fermion"].items()))


if __name__ == "__main__":
    main()