#!/usr/bin/env python3
"""M1 scripts/03 — R3 正規形與等價類：顯式編碼 + 簡併度

正規形全序（本 M 凍結）：
    (n_q, p, k, lex) ↑   n_q 量子數；p partial 數（0/1）；k 二值支撐；lex 索引序列字典序

編碼族（可精確計數）：
  F0（p=0）：S ⊆ {1..44}, |S| ≤ k, |Σ_S − t| ≤ 1e−3·t
  F1（p=1，top-mode 規範）：S ⊆ {1..43}, |S| ≤ k−1, 0 < t − Σ_S < γ₄₄ = 4, σ₄₄ = (t−Σ_S)/4 ∈ (0,1)
       ⟹ 每 S 恰一個編碼；σ₄₄ 由目標唯一決定（零額外自由度）

實現：子集和以 bitmask 倍增法一次生成（numpy 向量化），按支撐大小分組排序；計數用 searchsorted。
輸出：data/03_canonical.json
"""
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
KCAP = 20
FERMIONS = {"e": "0.51099895000", "mu": "105.6583755", "tau": "1776.86",
            "u": "2.16", "d": "4.67", "s": "93.4", "c": "1273.0", "b": "4180.0", "t": "172570.0"}

GAMMA_MP = [2 - 2 * cos(mpf(n) * pi / 44) for n in range(1, N + 1)]
GAMMA = np.array([float(g) for g in GAMMA_MP])
GMAX = float(GAMMA_MP[-1])


def build(mode_ids):
    """回傳 by_size[k] = (sorted sums, masks, mode_ids)；倍增法 O(2^m) 向量化"""
    sums = np.zeros(1)
    masks = np.zeros(1, dtype=np.uint64)
    sizes = np.zeros(1, dtype=np.int8)
    for bit, mid in enumerate(mode_ids):
        g = GAMMA[mid]
        sums = np.concatenate([sums, sums + g])
        masks = np.concatenate([masks, masks | (np.uint64(1) << np.uint64(bit))])
        sizes = np.concatenate([sizes, sizes + 1])
    out = []
    for k in range(0, len(mode_ids) + 1):
        sel = sizes == k
        s = sums[sel]
        order = np.argsort(s, kind="mergesort")
        out.append((s[order], masks[sel][order], np.array(mode_ids, dtype=np.int64)))
    return out


def count_window(tabs, lo_t, hi_t, kmax, strict_lo=False, strict_hi=False):
    A, B = tabs["A"], tabs["B"]
    total = 0
    for ka in range(0, min(kmax, len(A) - 1) + 1):
        kb_max = min(kmax - ka, len(B) - 1)
        sa = A[ka][0]
        if sa.size == 0:
            continue
        for kb in range(0, kb_max + 1):
            sb = B[kb][0]
            lo_i = np.searchsorted(sb, lo_t - sa, side="right" if strict_lo else "left")
            hi_i = np.searchsorted(sb, hi_t - sa, side="left" if strict_hi else "right")
            total += int(np.maximum(hi_i - lo_i, 0).sum())
    return total


def first_in_window(tabs, lo_t, hi_t, kmax, strict_lo=False, strict_hi=False):
    A, B = tabs["A"], tabs["B"]
    for ka in range(0, min(kmax, len(A) - 1) + 1):
        kb_max = min(kmax - ka, len(B) - 1)
        sa, ma, ida = A[ka]
        if sa.size == 0:
            continue
        for kb in range(0, kb_max + 1):
            sb, mb, idb = B[kb]
            lo_i = np.searchsorted(sb, lo_t - sa, side="right" if strict_lo else "left")
            hi_i = np.searchsorted(sb, hi_t - sa, side="left" if strict_hi else "right")
            hit = np.nonzero(hi_i > lo_i)[0]
            if hit.size:
                ia = int(hit[0])
                ib = int(lo_i[ia])
                return decode(int(ma[ia]), ida), decode(int(mb[ib]), idb), ka + kb
    return None, None, None


def decode(mask, ids):
    return [int(ids[b]) + 1 for b in range(len(ids)) if (mask >> np.uint64(b)) & np.uint64(1)]


def main():
    ids_all = list(range(N))
    tabs_all = {"A": build(ids_all[:22]), "B": build(ids_all[22:])}
    ids_n44 = list(range(N - 1))
    tabs_n44 = {"A": build(ids_n44[:22]), "B": build(ids_n44[22:])}

    out = {"meta": {"script": "03_canonical_form", "eps": "1e-3", "M0_MeV": str(M0),
                    "order": "(n_q, p, k, lex) up", "kcap": KCAP,
                    "F1": "top-mode 規範 σ₄₄=(t−Σ_S)/4 ∈ (0,1)，S ⊆ {1..43}"},
           "rows": {}}
    checks = [("tables_all_size", int(sum(a[0].size for a in tabs_all["A"]) + sum(a[0].size for a in tabs_all["B"])) == 2 ** 23),
              ("tables_n44_size", int(sum(a[0].size for a in tabs_n44["A"]) + sum(a[0].size for a in tabs_n44["B"])) == 2 ** 22 + 2 ** 21),
              ("gamma44=4", abs(GAMMA_MP[-1] - 4) < mpf("1e-60"))]

    for name, mass in FERMIONS.items():
        t_mp = mpf(mass) / M0
        t = float(t_mp)
        n_q = int(ceil(t_mp / 90))
        rec = {"target_units": mp.nstr(t_mp, 14), "n_q": n_q, "single_quantum": n_q == 1}
        if n_q != 1:
            rec["status"] = "🔴 跨量子：單量子 λ ≤ 90 不足，組合律未定（開放）"
            rec["per_quantum_mean_units"] = mp.nstr(t_mp / n_q, 10)
            out["rows"][name] = rec
            continue

        lo_w, hi_w = t * (1 - float(EPS)), t * (1 + float(EPS))
        k_star = next((k for k in range(1, KCAP + 1)
                       if count_window(tabs_all, lo_w, hi_w, k) >= 1), None)
        if k_star is not None:
            ca, cb, _ = first_in_window(tabs_all, lo_w, hi_w, k_star)
            enc = sorted(ca + cb)
            s = sum((GAMMA_MP[i - 1] for i in enc), mpf(0))
            rec.update({"p": 0, "k": k_star, "encoding": enc,
                        "sum_units": mp.nstr(s, 14), "dev": mp.nstr(abs(s - t_mp) / t_mp, 8),
                        "degeneracy": count_window(tabs_all, lo_w, hi_w, k_star),
                        "status": "✅ 二值編碼存在（F0）"})
        else:
            k_p = next((k for k in range(2, KCAP + 1)
                        if count_window(tabs_n44, t - GMAX, t, k - 1, strict_lo=True, strict_hi=True) >= 1), None)
            if k_p is None:
                rec.update({"p": None, "k": None, "status": "🔴 k ≤ 20 兩族皆無編碼"})
            else:
                ca, cb, _ = first_in_window(tabs_n44, t - GMAX, t, k_p - 1, strict_lo=True, strict_hi=True)
                enc = sorted(ca + cb)
                s = sum((GAMMA_MP[i - 1] for i in enc), mpf(0))
                sigma = (t_mp - s) / mpf(4)
                rec.update({"p": 1, "k": k_p, "encoding": enc, "sigma44": mp.nstr(sigma, 10),
                            "sum_units": mp.nstr(s, 14),
                            "dev": mp.nstr(abs(s + 4 * sigma - t_mp) / t_mp, 8),
                            "degeneracy": count_window(tabs_n44, t - GMAX, t, k_p - 1, strict_lo=True, strict_hi=True),
                            "status": "🟡 二值無解（k≤20）⟹ F1：1 partial（σ₄₄）"})
        out["rows"][name] = rec

    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["check_summary"] = {"total": len(checks), "pass": sum(1 for _, p in checks if p)}
    out["verdict"] = {
        "F0_exists": sum(1 for r in out["rows"].values() if r.get("status", "").startswith("✅")),
        "F1_forced": sum(1 for r in out["rows"].values() if r.get("status", "").startswith("🟡")),
        "multi_quantum": sum(1 for r in out["rows"].values() if not r.get("single_quantum", True)),
        "none": sum(1 for r in out["rows"].values() if r.get("status", "").startswith("🔴")),
    }

    with open(os.path.join(DATA, "03_canonical.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print(json.dumps({"checks": out["check_summary"], "verdict": out["verdict"]}, ensure_ascii=False))
    for n, r in out["rows"].items():
        if r.get("single_quantum"):
            print(f"{n:>3} n_q={r['n_q']} p={r['p']} k={r['k']} deg={r.get('degeneracy')} "
                  f"dev={r.get('dev','-')} σ44={r.get('sigma44','-')} enc={r.get('encoding')}")
        else:
            print(f"{n:>3} n_q={r['n_q']} {r['status']} (每量子均值 {r['per_quantum_mean_units']} M₀)")


if __name__ == "__main__":
    main()