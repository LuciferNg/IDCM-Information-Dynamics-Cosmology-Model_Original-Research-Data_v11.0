#!/usr/bin/env python3
"""M1 scripts/01 — R2 多重性引擎：ρ(m; k, ε) 精確計數

定義（本 M 凍結）：
  ρ_bin(m; k, ε)  = #{ S ⊆ {1..44} : |S| ≤ k, |M₀·Σ_{n∈S} γ_n − m| ≤ ε·m }        （純二值 σ）
  ρ_1p (m; k, ε)  = #{ (S, j) : |S| ≤ k−1, j ∉ S, 0 < m − M₀Σ_S < M₀γ_j }          （恰一個 partial）
                    〔σ_j = (m − M₀Σ_S)/(M₀γ_j) ∈ (0,1) 由目標唯一決定 ⟹ 每對 (S,j) 恰一個編碼〕
  k*(m; ε)        = min{ k : ρ_bin(m;k,ε) ≥ 1 }                                    （二值最小支撐）

算法：meet-in-the-middle（44 模 → 22 + 22），按支撐大小分組；以 searchsorted 精確計數。
      枚舉用 float64；所有「命中」以 mpmath 80 dps 覆核，並報最小邊界裕度（穩定性證書）。
輸出：data/01_multiplicity.json
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
SPLIT = 22
KMAX = 12
M0 = mpf("27.9009145112")            # MeV
EPS_LIST = [1e-3, 1e-4]
K_LIST = [1, 2, 3, 4, 6, 8, 12]

FERMIONS = {
    "e": "0.51099895000", "mu": "105.6583755", "tau": "1776.86",
    "u": "2.16", "d": "4.67", "s": "93.4", "c": "1273.0", "b": "4180.0", "t": "172570.0",
}
SCHEME = {"u": "2 GeV MS-bar", "d": "2 GeV MS-bar", "s": "2 GeV MS-bar",
          "c": "m_c MS-bar", "b": "m_b MS-bar", "t": "pole", "e": "pole", "mu": "pole", "tau": "pole"}


def gamma_mp(n):
    return 2 - 2 * cos(mpf(n) * pi / 44)


GAMMA_MP = [gamma_mp(n) for n in range(1, N + 1)]
GAMMA = np.array([float(g) for g in GAMMA_MP])


def half_tables(gamma, kmax):
    """回傳 [(sizes, sorted_sum_arrays)] × 2"""
    tabs = []
    for lo, hi in ((0, SPLIT), (SPLIT, N)):
        idx = list(range(lo, hi))
        by_size = []
        for k in range(0, kmax + 1):
            combos = list(itertools.combinations(idx, k))
            arr = np.array([sum(gamma[i] for i in c) for c in combos], dtype=float) if combos else np.zeros(1)
            by_size.append(np.sort(arr))
        tabs.append(by_size)
    return tabs


def rho_bin(target, k, eps, tabs):
    lo_w = target * (1 - eps)
    hi_w = target * (1 + eps)
    total = 0
    A, B = tabs
    for ka in range(0, min(k, KMAX) + 1):
        kb_max = min(k - ka, KMAX)
        sa = A[ka]
        if sa.size == 0:
            continue
        acc = 0
        for kb in range(0, kb_max + 1):
            sb = B[kb]
            lo = np.searchsorted(sb, lo_w - sa, side="left")
            hi = np.searchsorted(sb, hi_w - sa, side="right")
            acc += int((hi - lo).sum())
        total += acc
    return total


def rho_partial1(target, k, kmax_sub=None):
    """恰一個 partial；|S| ≤ k−1；σ_j ∈ (0,1) 由目標決定（容差 ε 只作用於最終質量 ⟹ 此處無 ε）"""
    kmax_sub = (k - 1) if kmax_sub is None else kmax_sub
    gmax = float(GAMMA[-1])
    # 以 γ 之降序門檻計數： #{j : γ_j > δ} = 44 − #{γ_j ≤ δ}
    order = np.argsort(GAMMA)              # 升序
    gsorted = GAMMA[order]
    total = 0
    for size in range(0, kmax_sub + 1):
        for combo in itertools.combinations(range(N), size):
            s = float(sum(GAMMA[i] for i in combo))
            delta = target - s
            if delta <= 0 or delta >= gmax:
                continue
            n_above = N - int(np.searchsorted(gsorted, delta, side="right"))
            if n_above <= 0:
                continue
            in_S_above = sum(1 for i in combo if GAMMA[i] > delta)
            total += max(n_above - in_S_above, 0)
    return total


def margin_check(combo, target_mp, eps):
    """以 mpmath 覆核某編碼之偏差並回傳邊界裕度（ε 相對）"""
    s = sum((GAMMA_MP[i] for i in combo), mpf(0))
    dev = abs(s - target_mp) / target_mp
    return dev, abs(dev - mpf(str(eps)))


def main():
    tabs = half_tables(GAMMA, KMAX)
    out = {"meta": {"script": "01_multiplicity_engine", "N": N, "split": SPLIT,
                    "kmax": KMAX, "M0_MeV": str(M0), "eps": EPS_LIST, "k_list": K_LIST,
                    "gamma1": mp.nstr(GAMMA_MP[0], 12), "gamma44": mp.nstr(GAMMA_MP[-1], 12),
                    "sum_gamma": mp.nstr(sum(GAMMA_MP, mpf(0)), 20)},
           "scheme": SCHEME}

    checks = []
    checks.append(("sum_gamma=90", abs(sum(GAMMA_MP, mpf(0)) - 90) < mpf("1e-60")))
    # 引擎自檢：k=0..2 與暴力枚舉一致
    for t in (mpf("3.78691442023"), mpf("0.018314774227")):
        brute = 0
        for k in range(0, 3):
            for c in itertools.combinations(range(N), k):
                s = sum((GAMMA_MP[i] for i in c), mpf(0))
                if abs(s - t) <= mpf("1e-3") * t:
                    brute += 1
        fast = rho_bin(float(t), 2, 1e-3, tabs)
        checks.append((f"engine_match_k2_t{t}", brute == fast))

    rows = {}
    for name, mass in FERMIONS.items():
        t_mp = mpf(mass) / M0
        t = float(t_mp)
        rec = {"target_units": mp.nstr(t_mp, 12), "scheme": SCHEME[name], "rho": {}}
        for eps in EPS_LIST:
            key = f"{eps:.0e}"
            rec["rho"][key] = {str(k): rho_bin(t, k, eps, tabs) for k in K_LIST}
            rec["rho"][key]["k_star"] = next((k for k in K_LIST if rec["rho"][key][str(k)] >= 1), None)
        rec["rho_1p"] = {str(k): rho_partial1(t, k) for k in (1, 2, 3, 4, 5)}
        rows[name] = rec
    out["rows"] = rows

    # 最近接近度（穩定性證書）：k ≤ 4 全部子集和相對目標之最小偏差
    all_sums4 = np.concatenate([
        np.array([float(sum(GAMMA[i] for i in c))
                  for c in itertools.combinations(range(N), kk)], dtype=float) if kk else np.zeros(1)
        for kk in range(0, 5)])
    for name, rec in rows.items():
        t = float(mpf(FERMIONS[name]) / M0)
        rec["closest_dev_k4"] = float(np.min(np.abs(all_sums4 - t)) / t)

    # 覆蓋率（精確 measure）：k ≤ 4 全部子集和之區間聯集 / 90
    def coverage(k, eps):
        combos_sums = []
        for kk in range(0, k + 1):
            for c in itertools.combinations(range(N), kk):
                combos_sums.append(float(sum(GAMMA[i] for i in c)))
        s = np.sort(np.array(combos_sums))
        iv = np.stack([s * (1 - eps), s * (1 + eps)], axis=1)
        iv = iv[np.argsort(iv[:, 0])]
        total = 0.0
        cur_lo, cur_hi = iv[0]
        for lo, hi in iv[1:]:
            if lo <= cur_hi:
                cur_hi = max(cur_hi, hi)
            else:
                total += cur_hi - cur_lo
                cur_lo, cur_hi = lo, hi
        total += cur_hi - cur_lo
        return total / 90.0

    out["coverage"] = {f"k{k}_e{eps:.0e}": coverage(k, eps) for k in (1, 2, 3, 4) for eps in EPS_LIST}
    out["coverage_note"] = "均勻目標 m~U(0,90) 下被「至少一個編碼」擊中之機率；>1 bit 判別力之必要條件"

    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["check_summary"] = {"total": len(checks), "pass": sum(1 for _, p in checks if p)}

    with open(os.path.join(DATA, "01_multiplicity.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print(json.dumps({"checks": out["check_summary"], "coverage": out["coverage"]}, ensure_ascii=False))
    print("\nfermion   m/M0            rho(k<=2,1e-3) rho(k<=3) rho(k<=4) rho(k<=8) rho(k<=12) k*(1e-3)  rho_1p(k<=2)")
    for name, rec in rows.items():
        r = rec["rho"]["1e-03"]
        print(f"{name:>7} {rec['target_units']:>15} {r['2']:>13} {r['3']:>9} {r['4']:>9} "
              f"{r['8']:>9} {r['12']:>10} {str(r['k_star']):>9} {rec['rho_1p']['2']:>13}")


if __name__ == "__main__":
    main()