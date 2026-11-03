#!/usr/bin/env python3
"""M1 scripts/06 — R6 跨量子組合律：可識別性判定

背景：b、t 之 m/M₀ 超過單量子帳本 Σγ = 90 ⟹ n_q ≥ ⌈t/90⌉（M6 §3）。
檢定：
  A. 兩條零參數候選律（等分 E／滿載+餘額 F）能否精確重現 b、t（可識別性）
  B. 由候選律導出之量子含量（h/r/c）在譜層有無結構身份（密度閘 p ≤ 0.1）
  C. 預註冊 R6-Q1（每量子 ≥ 86）、R6-Q2（兩餘額同源）之對帳

輸出：data/06_multiquantum.json
"""
import itertools
import json
import os

import numpy as np
from mpmath import mp, mpf, cos, pi, ceil

mp.dps = 80
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "data"))
os.makedirs(DATA, exist_ok=True)

N = 44
EPS = mpf("1e-3")
M0 = mpf("27.9009145112")
SUM_G = mpf(90)

GAMMA_MP = [2 - 2 * cos(mpf(n) * pi / 44) for n in range(1, N + 1)]
GAMMA = np.array([float(g) for g in GAMMA_MP])
G_DESC = sorted(GAMMA_MP, reverse=True)
LAM = {k: sum(G_DESC[:k], mpf(0)) for k in range(1, 45)}


def k_min(t):
    for k in range(1, 45):
        if LAM[k] >= t:
            return k
    return None


TABS = []
for lo, hi in ((0, 22), (22, 44)):
    sums = np.zeros(1)
    sizes = np.zeros(1, dtype=np.int8)
    for mid in range(lo, hi):
        sums = np.concatenate([sums, sums + GAMMA[mid]])
        sizes = np.concatenate([sizes, sizes + 1])
    TABS.append([np.sort(sums[sizes == k]) for k in range(0, hi - lo + 1)])


def count_window(lo_t, hi_t, kmax):
    A, B = TABS
    total = 0
    for ka in range(0, min(kmax, len(A) - 1) + 1):
        kb_max = min(kmax - ka, len(B) - 1)
        sa = A[ka]
        if sa.size == 0:
            continue
        for kb in range(0, kb_max + 1):
            sb = B[kb]
            total += int((np.searchsorted(sb, hi_t - sa, side="right") -
                          np.searchsorted(sb, lo_t - sa, side="left")).sum())
    return total


def coverage(k, eps):
    """精確 measure：k 支撐全部子集和之 ε-區間聯集 / 90（暴力枚舉，k ≤ 3）"""
    vals = []
    for kk in range(0, k + 1):
        for c in itertools.combinations(range(N), kk):
            vals.append(float(sum(GAMMA[i] for i in c)))
    s = np.sort(np.array(vals))
    iv = np.stack([s * (1 - eps), s * (1 + eps)], axis=1)
    iv = iv[np.argsort(iv[:, 0])]
    total, cur_lo, cur_hi = 0.0, iv[0][0], iv[0][1]
    for lo, hi in iv[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            total += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    return (total + cur_hi - cur_lo) / 90.0


def main():
    out = {"meta": {"script": "06_multiquantum", "M0_MeV": str(M0), "sum_gamma": 90,
                    "eps": str(EPS)}, "checks": []}
    checks = []

    m_b = mpf("4180.0") / M0
    m_t = mpf("172570.0") / M0
    nq_b = int(ceil(m_b / SUM_G))
    nq_t = int(ceil(m_t / SUM_G))
    checks.append(("nq_b=2", nq_b == 2))
    checks.append(("nq_t=69", nq_t == 69))

    targets = {
        "h_b": nq_b * SUM_G - m_b,          # 滿載讀法之「洞」
        "r_b": m_b - (nq_b - 1) * SUM_G,    # 滿載+餘額讀法之餘額
        "h_t": nq_t * SUM_G - m_t,
        "r_t": m_t - (nq_t - 1) * SUM_G,
        "c_b": m_b / nq_b,                  # 等分讀法之每量子含量
        "c_t": m_t / nq_t,
    }
    for k, v in targets.items():
        checks.append((f"{k}_in_range", 0 < v <= 90))

    eps_f = float(EPS)
    rows = {}
    for name, t in targets.items():
        rec = {"value_units": mp.nstr(t, 14), "k_min": k_min(t)}
        devs = [(abs(g - t) / t, n) for n, g in enumerate(GAMMA_MP, 1)]
        bd, bn = min(devs)
        rec["single_mode_best"] = {"n": bn, "dev": mp.nstr(bd, 6)}
        tf = float(t)
        rec["rho_kmin"] = count_window(tf * (1 - eps_f), tf * (1 + eps_f), rec["k_min"])
        rec["rho_k2"] = count_window(tf * (1 - eps_f), tf * (1 + eps_f), 2)
        rec["t_over_gamma1"] = mp.nstr(t / GAMMA_MP[0], 10)
        rows[name] = rec
    out["targets"] = rows

    cov = {f"k{k}": coverage(k, eps_f) for k in (1, 2, 3)}
    out["coverage"] = cov

    out["laws"] = {
        "E_equal_split": {"b": "2 × 74.90793892 = 149.815877839", "t": "69 × 89.63915831 = 6185.10192312",
                          "dof": "σ 共用（44 維向量），分布律零新常數"},
        "F_fill_plus_remainder": {"b": "90 + 59.815877839", "t": "68 × 90 + 65.10192312",
                                  "dof": "σ 共用；僅末量子非滿"},
    }

    ident = {}
    for name, rec in rows.items():
        p_single = cov["k1"] if float(rec["single_mode_best"]["dev"]) <= eps_f else None
        p_k2 = cov["k2"] if rec["rho_k2"] >= 1 else None
        if p_single is not None and p_single <= 0.1:
            judged = "命中（單模，p=%.4f）" % p_single
        elif rec["rho_k2"] == 1 and p_k2 is not None and p_k2 <= 0.1:
            judged = "命中（k≤2 且唯一）"
        else:
            judged = "無鑑別力"
        ident[name] = {"single_mode_dev": rec["single_mode_best"]["dev"], "p_single": p_single,
                       "k2_hit": rec["rho_k2"] >= 1, "p_k2": p_k2, "judged": judged}
    out["structural_identity"] = ident

    out["R6"] = {
        "R6_Q1": {"claim": "每量子 ≥ 86 = 90 − γ₄₄",
                  "E_contents": {"b": mp.nstr(targets["c_b"], 10), "t": mp.nstr(targets["c_t"], 10)},
                  "F_last_quantum": {"b": mp.nstr(targets["r_b"], 10), "t": mp.nstr(targets["r_t"], 10)},
                  "verdict": "🔴 FAIL",
                  "reason": "E 下 b 每量子 = 74.9079 < 86；F 下末量子 59.8159 < 86 ⟹ 「≥86」係由 t 均值外推，無結構依據"},
        "R6_Q2": {"claim": "兩餘額同源（同一算子作用 n_q 次）",
                  "per_quantum_deficit_b": mp.nstr(targets["h_b"] / nq_b, 8),
                  "per_quantum_deficit_t": mp.nstr(targets["h_t"] / nq_t, 8),
                  "ratio": mp.nstr((targets["h_b"] / nq_b) / (targets["h_t"] / nq_t), 8),
                  "verdict": "🔴 FAIL",
                  "reason": "常數虧損律要求 h_b/2 = h_t/69；實測 15.0920611 vs 0.3608243（比 41.826）⟹ 非同源"},
    }

    out["identifiability"] = {
        "verdict": "🔴 兩律（E／F）皆以零新常數精確重現 b、t ⟹ 質量值層不可識別",
        "reason": "σ ∈ [0,1]^44 連續 ⟹ 任何分配皆可表；n_q 本身不可觀測",
    }
    out["conclusion"] = {
        "verdict": "🔴 跨量子組合律：no-go（不可識別）",
        "status": "本 M 維持 partial；觸發器 = 未來獨立結構提供量子分配機制（如 M11 通道算子）",
        "禁區": ["D1 不得由 b/t 質量值反解量子分配",
                "D2 不得以 γ₁ 格點命中當結構（間距 γ₁ = 0.0051 ⟹ 1e−3 容差下必然多重命中）"],
    }

    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["check_summary"] = {"total": len(checks), "pass": sum(1 for _, p in checks if p)}

    with open(os.path.join(DATA, "06_multiquantum.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print(json.dumps({"checks": out["check_summary"], "coverage": cov}, ensure_ascii=False))
    print(f"\n{'量':>4} {'值(M₀)':>14} {'k_min':>5} {'單模 dev':>10} {'ρ(k_min)':>9} {'ρ(k≤2)':>7}  判")
    for name, rec in rows.items():
        print(f"{name:>4} {rec['value_units']:>14} {rec['k_min']:>5} {rec['single_mode_best']['dev']:>10} "
              f"{rec['rho_kmin']:>9} {rec['rho_k2']:>7}  {ident[name]['judged']}")
    print("\nR6-Q1", out["R6"]["R6_Q1"]["verdict"], "｜R6-Q2", out["R6"]["R6_Q2"]["verdict"])
    print("per-quantum deficit b / t =", out["R6"]["R6_Q2"]["per_quantum_deficit_b"], "/",
          out["R6"]["R6_Q2"]["per_quantum_deficit_t"], "｜ratio", out["R6"]["R6_Q2"]["ratio"])


if __name__ == "__main__":
    main()