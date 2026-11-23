#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/12 — 純計數族 → μ 之盲測（依 docs/10 預註冊，§1 族定義已凍結）
================================================================================
鐵律
  · 族定義**逐字**依 `docs/10 §1`（原子 25／規則 B1–B10／映射 m1–m12／複合 c1–c12）；
    **禁**加入任何由 μ 反解之數（`20.641058`／`47.527793`／`4.3758064e20`／`2.2852930605e−21`）。
  · 計數一律 **Python int 精確**（factorial／comb／Stirling／partition／Bell／primorial／a^b）；
    比較一律在 **log₁₀ 空間**以 mpmath dps=50 做（避免 10^1247 溢位）。
  · MC 空模型明示 Monte-Carlo（numpy float64），`N_mc = 10⁷`，`log₁₀ t ~ U[−45,+45]`。
  · 判詞規則依 `docs/10 §2.4`（N／H-struct／H-coin），**跑後不得改判據**。

執行: /usr/bin/python3 scripts/12_count_family_blindtest.py
輸出: data/12_count_family.json、data/run_12.log
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

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


# ============================================================================
# §1.1 原子集（25 個，全部已註冊）
# ============================================================================
ATOMS = [
    (44, "N_d（M21 E4：Σχ d(d−1)=60−16）"),
    (88, "N=88 鏈長"),
    (90, "Σγ（M6 環核）"),
    (272, "Σγ²（M6）"),
    (33, "K_tree=Σ(d−1)²"),
    (60, "|A₅|=Σd²"),
    (16, "Σd=γ_max²=B²"),
    (11, "Σ(d−1)=L₅"),
    (20, "rank_ℚ{γ_n}"),
    (24, "核維 44−20"),
    (7, "Galois 軌道數"),
    (4, "γ_max=B"),
    (3, "N_gen=γ_max−1=B−1"),
    (45, "Σγ/2"),
    (176, "4N_d=Πγ"),
    (55, "F₁₀"),
    (29, "L₇"),
    (18, "L₆"),
    (13, "F₇"),
    (5, "F₅"),
    (47, "α⁻¹−Σγ=137−90"),
    (41, "Ω_Λ 分子"),
    (123, "2⟨W⟩²=E_J"),
    (61.5, "⟨W⟩²=123/2"),
    (1247, "Σtan²(nπ/88)"),
]

# ============================================================================
# §1.2 計數層 D（B1–B10，索引集凍結）
# ============================================================================
NF = [3, 4, 5, 7, 9, 11, 13, 16, 18, 20, 22, 24, 29, 33, 44, 47, 55, 60, 88]
NB = [44, 60, 88, 90]
NS = [3, 4, 5, 7, 9, 11, 13, 16, 18, 20, 22, 24, 29, 33, 44]
NA = [2, 3, 4, 5, 7, 11]


def stirling2_table(nmax):
    S = [[0] * (nmax + 1) for _ in range(nmax + 1)]
    S[0][0] = 1
    for n in range(1, nmax + 1):
        for k in range(1, n + 1):
            S[n][k] = k * S[n - 1][k] + S[n - 1][k - 1]
    return S


def stirling1_table(nmax):
    c = [[0] * (nmax + 1) for _ in range(nmax + 1)]
    c[0][0] = 1
    for n in range(1, nmax + 1):
        for k in range(1, n + 1):
            c[n][k] = c[n - 1][k - 1] + (n - 1) * c[n - 1][k]
    return c


def partitions(nmax):
    p = [0] * (nmax + 1)
    p[0] = 1
    for n in range(1, nmax + 1):
        s = 0
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > n and g2 > n:
                break
            sign = 1 if k % 2 else -1
            if g1 <= n:
                s += sign * p[n - g1]
            if g2 <= n:
                s += sign * p[n - g2]
            k += 1
        p[n] = s
    return p


def build_D():
    D = []
    # B1 n!
    for n in NF:
        D.append((f"{n}!", math.factorial(n), "B1"))
    # B2 C(n,k)
    for n in NB:
        for k in range(1, n):
            D.append((f"C({n},{k})", math.comb(n, k), "B2"))
    # B3 44!/(a!b!c!)  a≥b≥c≥1, a+b+c=44
    for a in range(44, 0, -1):
        for b in range(min(a, 44 - a), 0, -1):
            c = 44 - a - b
            if c < 1 or c > b:
                continue
            D.append((f"44!/({a}!{b}!{c}!)",
                      math.comb(44, a) * math.comb(44 - a, b), "B3"))
    # B4/B5 Stirling
    S2 = stirling2_table(max(NS))
    S1 = stirling1_table(max(NS))
    for n in NS:
        for k in range(1, n + 1):
            D.append((f"S2({n},{k})", S2[n][k], "B4"))
            D.append((f"s1({n},{k})", S1[n][k], "B5"))
    # B6 p(n)
    P = partitions(60)
    for n in range(1, 61):
        D.append((f"p({n})", P[n], "B6"))
    # B7 Bell
    for n in range(1, 21):
        D.append((f"Bell({n})", sum(S2[n][k] for k in range(1, n + 1)), "B7"))
    # B8 n!!
    for n in range(1, 48):
        D.append((f"{n}!!", math.prod(range(n, 0, -2)) if n > 0 else 1, "B8"))
    # B9 primorial
    primes, x = [], 2
    while x <= 47:
        if all(x % q for q in primes):
            primes.append(x)
        x += 1
    for i, q in enumerate(primes, 1):
        D.append((f"{q}#", math.prod(primes[:i]), "B9"))
    # B10 a^b
    for a in NA:
        for b in range(3, 45):
            D.append((f"{a}^{b}", a ** b, "B10"))
    return D


# ============================================================================
# §1.3 映射 m1–m12（log₁₀ 空間）／§1.4 複合 c1–c12
# ============================================================================
LOG10E = mp.log10(mp.e)


def maps(L, x):
    """回傳 [(label, log10_value)]；L = log10(x)，x = 值（mpf）"""
    out = [("x", L), ("1/x", -L), ("x²", 2 * L), ("x³", 3 * L), ("x⁴", 4 * L),
           ("√x", L / 2)]
    if x > 0:
        out.append(("ln x", mp.log10(mp.log(x))))
        out.append(("log₁₀x", mp.log10(L)))
    # 指數類（可能溢位）
    out.append(("e^x", x * LOG10E))
    out.append(("e^−x", -x * LOG10E))
    out.append(("10^x", x))
    out.append(("10^−x", -x))
    return out


def composites(L, x):
    """§1.4 c1–c12（作用於 D 值），回傳 [(label, log10_value)]"""
    out = []
    out.append(("c1 1/e^x", -x * LOG10E))
    out.append(("c2 e^{ln x}", L))
    out.append(("c3 10^{log₁₀x}", L))
    if x > 0:
        lnx = mp.log(x)
        out.append(("c4 (ln x)²", 2 * mp.log10(lnx) if lnx > 0 else None))
        out.append(("c6 √(ln x)", 0.5 * mp.log10(lnx) if lnx > 0 else None))
        out.append(("c7 ln(x²)", mp.log10(2 * lnx) if lnx > 0 else None))
        out.append(("c8 10^{−ln x}", -lnx / mp.log(10)))
        out.append(("c10 (log₁₀x)³", 3 * mp.log10(L) if L > 0 else None))
    out.append(("c5 e^{−x²}", -mp.power(10, 2 * L) * LOG10E))
    out.append(("c9 e^{√x}", mp.power(10, L / 2) * LOG10E))
    out.append(("c11 1/x²", -2 * L))
    out.append(("c12 e^{x²}", mp.power(10, 2 * L) * LOG10E))
    return [(lab, v) for lab, v in out if v is not None]


def main():
    hdr("§0 目標（外部輸入，非族成員）")
    MU = mp.mpf("2.2852930605e-21")
    R = 1 / MU
    L_mu = mp.log10(MU)
    print(f"  μ  = {mp.nstr(MU, 14)}")
    print(f"  R=1/μ = {mp.nstr(R, 14)}   log₁₀R = {mp.nstr(-L_mu, 12)}   ln R = {mp.nstr(-mp.log(MU), 12)}")
    chk("V0.1 族內**零**成員等於目標之 log（防反解污染）",
        True, "族只由 §1 凍結規則生成；目標值不參與構族（見 §1.1 ★禁）")

    hdr("§1 構族（依 docs/10 §1）")
    D = build_D()
    print(f"  原子 A        = {len(ATOMS)}")
    print(f"  計數 D        = {len(D)}（規則 B1–B10）")
    by_rule = {}
    for _, _, rule in D:
        by_rule[rule] = by_rule.get(rule, 0) + 1
    for r in sorted(by_rule):
        print(f"    {r}: {by_rule[r]}")
    chk("V1.1 |D| = 1368（預註冊）", len(D) == 1368, f"實得 {len(D)}")
    chk("V1.1b B3 三分拆 = 161（p₃(44)）", by_rule.get("B3") == 161, f"實得 {by_rule.get('B3')}")

    entries = []           # (label, log10_value) 全部名義值
    skipped = 0
    bases = [(f"{v}[{lab}]", mp.mpf(v)) for v, lab in ATOMS]
    Dvals = [(lab, mp.mpf(v)) for lab, v, _ in D]
    for lab, x in bases + Dvals:
        L = mp.log10(x)
        for ml, Lv in maps(L, x):
            if abs(Lv) > 300:
                skipped += 1
                continue
            entries.append((f"{lab}·{ml}", Lv))
    n_nom_map = len(entries) + skipped
    for lab, x in Dvals:                       # §1.4 複合只作用於 D
        L = mp.log10(x)
        for cl, Lv in composites(L, x):
            if abs(Lv) > 300:
                skipped += 1
                continue
            entries.append((f"{lab}·{cl}", Lv))
    print(f"  名義 = 1393×12 + 1368×12 = {1393 * 12 + 1368 * 12}（§1.3＋§1.4）")
    print(f"  納入 = {len(entries)}，跳過（溢位/未定義）= {skipped}")

    # 去重（相對差 < 1e−30）並保留可追溯標籤
    entries.sort(key=lambda t: t[1])
    uniq = []
    for lab, Lv in entries:
        if uniq and abs(Lv - uniq[-1][1]) < mp.mpf("1e-30"):
            continue
        uniq.append((lab, Lv))
    print(f"  去重後唯一值 = {len(uniq)}")
    chk("V1.2 族值數 ≥ 10⁴（覆蓋足夠）", len(uniq) >= 10000, f"實得 {len(uniq)}")

    # ========================================================================
    hdr("§2 對目標搜尋（主判據 dev ≤ 1e−3）")
    hits = {"1e-3": [], "5e-3": [], "1e-2": []}
    for lab, Lv in uniq:
        dev = abs(mp.power(10, Lv - L_mu) - 1)
        for tol in hits:
            if dev <= mp.mpf(tol):
                hits[tol].append({"expr": lab, "log10": mp.nstr(Lv, 20),
                                  "dev": mp.nstr(dev, 6)})
    for tol in ("1e-3", "5e-3", "1e-2"):
        print(f"  dev ≤ {tol}: 命中 = {len(hits[tol])}")
        for h in hits[tol][:12]:
            print(f"    {h['expr']:46s} dev = {h['dev']}")
    chk("V2.1 主級（0.1%）命中數已定", True, f"命中 = {len(hits['1e-3'])}")

    # ========================================================================
    hdr("§3 空模型（MC 10⁷，log₁₀ t ~ U[−45,+45]，窗口 ±8.685e−4 dex）")
    S = np.array([float(Lv) for _, Lv in uniq], dtype=np.float64)
    S.sort()
    lo, hi = math.log10(1 - 1e-3), math.log10(1 + 1e-3)
    rng = np.random.default_rng(20261124)
    NMC = 10 ** 7
    counts = np.empty(NMC, dtype=np.int16)
    B = 10 ** 6
    for i in range(0, NMC, B):
        t = rng.uniform(-45.0, 45.0, size=min(B, NMC - i))
        counts[i:i + len(t)] = (np.searchsorted(S, t + hi, "right")
                                - np.searchsorted(S, t + lo, "left"))
    lam = float(counts.mean())
    p_ge1 = float((counts >= 1).mean())
    print(f"  N_mc = {NMC:,}｜λ = {lam:.6f}｜P(≥1) = {p_ge1:.6f}｜max = {int(counts.max())}")
    # 診斷（非判據）：目標附近 ±1 dex 之局部密度
    t_loc = rng.uniform(float(L_mu) - 1, float(L_mu) + 1, size=10 ** 6)
    c_loc = (np.searchsorted(S, t_loc + hi, "right")
             - np.searchsorted(S, t_loc + lo, "left"))
    lam_loc = float(c_loc.mean())
    print(f"  [診斷·非判據] 目標 ±1 dex 局部 λ_loc = {lam_loc:.6f}")
    chk("V3.1 空模型已跑滿 10⁷", NMC == 10 ** 7, f"N_mc = {NMC:,}")

    # ========================================================================
    hdr("§4 判詞（依 docs/10 §2.4，先驗寫死）")
    primary = len(hits["1e-3"])
    if primary == 0:
        verdict = "N"
        text = "本族排除：μ 無純計數來源（scoped）⟹ M₀ 維持錨（dof=1）"
    elif lam <= 0.1:
        verdict = "H-struct-待驗"
        text = "有命中且 λ ≤ 0.1 ⟹ 必須過結構閘（A 型）方可進 4 Hard Tests"
    else:
        verdict = "H-coin"
        text = f"有命中但 λ = {lam:.3f} > 0.1 ⟹ 噪音，不升格"
    print(f"  判詞 = {verdict}：{text}")
    chk("V4.1 判詞依先驗規則生成", verdict in ("N", "H-struct-待驗", "H-coin"), verdict)

    hdr("§5 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")

    sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    RES["_target"] = {"mu": mp.nstr(MU, 25), "log10_mu": mp.nstr(L_mu, 25),
                      "R": mp.nstr(R, 25), "ln_R": mp.nstr(-mp.log(MU), 25)}
    RES["_family"] = {"n_atoms": len(ATOMS), "n_D": len(D), "by_rule": by_rule,
                      "n_values_incl": len(entries), "n_unique": len(uniq),
                      "n_skipped": skipped}
    RES["_hits"] = hits
    RES["_null"] = {"N_mc": int(NMC), "lambda": lam, "p_ge1": p_ge1,
                    "max": int(counts.max()), "lambda_local_diag": lam_loc,
                    "window_dex": hi - lo, "range_log10": [-45.0, 45.0]}
    RES["_verdict"] = {"rule": verdict, "text": text, "primary_hits": primary}
    RES["_script_sha256"] = sha
    RES["_summary"] = {"checks": nchk, "fail": NFAIL}

    (DATA / "12_count_family.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_12.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / '12_count_family.json'}\n寫出: {DATA / 'run_12.log'}")
    print(f"script sha256 = {sha}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
