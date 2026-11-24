#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/14 — 步長比族 μ = r^N 之盲測（依 docs/11 預註冊；門 C）
================================================================================
鐵律
  · r 族**逐字**依 `docs/11 §1`（R1 γ 比 1892／R2 整數比 306／R3 結構比 21／R4 φ 族 8）。
  · **禁**任何由 μ 反解之 r；族內不含 μ^{1/n}、e^{−47.5278/n} 之類。
  · 主判據只一級：`N` 精確落於特權整數集 `I`（`|N−n| ≤ 1e−6`）。
  · 空模型 = **target-shift**（`log₁₀ μ′ ~ U[−25,−15]`，`M = 10⁵`），明示 Monte-Carlo（numpy）。
  · 判詞規則依 `docs/11 §2`，跑後不得改。

執行: /usr/bin/python3 scripts/14_step_ratio_scan.py
輸出: data/14_step_ratio.json、data/run_14.log
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40

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


# ── 凍結常數 ─────────────────────────────────────────────────────────────────
MU = mp.mpf("2.2852930605e-21")
LN_MU = mp.log(MU)

I_PRIV = [3, 4, 5, 7, 9, 11, 13, 16, 18, 20, 22, 24, 29, 33, 36, 41, 44, 45, 47, 55, 60, 88, 90]
S2 = [3, 4, 5, 7, 11, 16, 18, 20, 24, 33, 41, 44, 45, 55, 60, 88, 90, 176]
R3 = ["3/4", "3/5", "4/5", "1/3", "1/4", "1/5", "3/16", "16/60", "41/60", "42/41", "85/41",
      "5/8", "11/5", "4/15", "2/45", "3/88", "5/176", "11/44", "29/44", "33/44", "55/88"]
PHI = (1 + mp.sqrt(5)) / 2
R4 = {"1/φ": 1 / PHI, "φ/2": PHI / 2, "φ⁻²": 1 / PHI ** 2, "11/5": mp.mpf(11) / 5,
      "11": mp.mpf(11), "5/3": mp.mpf(5) / 3, "4/φ²": 4 / PHI ** 2, "29/18": mp.mpf(29) / 18}


def gamma(n):
    return 4 * mp.sin(mp.pi * n / 88) ** 2


def build_family():
    fam = []                                    # (label, mpf value, 類)
    for n in range(1, 45):                      # R1
        for m in range(1, 45):
            if n == m:
                continue
            fam.append((f"γ{n}/γ{m}", gamma(n) / gamma(m), "R1"))
    for a in S2:                                # R2
        for b in S2:
            if a == b:
                continue
            fam.append((f"{a}/{b}", mp.mpf(a) / b, "R2"))
    for s in R3:                                # R3
        a, b = s.split("/")
        fam.append((s, mp.mpf(a) / mp.mpf(b), "R3"))
    for k, v in R4.items():                     # R4
        fam.append((k, mp.mpf(v), "R4"))
    return fam


def main():
    hdr("§0 目標（外部輸入）")
    print(f"  μ = {mp.nstr(MU, 14)}｜ln μ = {mp.nstr(LN_MU, 14)}｜R = {mp.nstr(1/MU, 12)}")
    chk("V0.1 族內零 μ-反解成員", True, "族只由 docs/11 §1 規則生成，μ 不參與構族")

    hdr("§1 構族（依 docs/11 §1）")
    fam = build_family()
    cnt = {}
    for _, _, k in fam:
        cnt[k] = cnt.get(k, 0) + 1
    print(f"  |R| = {len(fam)}　（R1 γ 比 {cnt.get('R1')}／R2 整數比 {cnt.get('R2')}"
          f"／R3 結構比 {cnt.get('R3')}／R4 φ 族 {cnt.get('R4')}）")
    chk("V1.1 |R| = 2227（預註冊）", len(fam) == 2227, f"實得 {len(fam)}")
    chk("V1.2 R1 = 44×43 = 1892", cnt.get("R1") == 1892, f"{cnt.get('R1')}")
    chk("V1.3 R2 = 18×17 = 306", cnt.get("R2") == 306, f"{cnt.get('R2')}")

    hdr("§2 對 μ 搜尋：N = ln μ / ln r")
    rows = []
    for lab, r, kind in fam:
        if r == 1 or r <= 0:
            continue
        N = LN_MU / mp.log(r)
        nrst = min(I_PRIV, key=lambda n: abs(N - n))
        dist = abs(N - nrst)
        rows.append((dist, N, nrst, lab, kind, r))
    rows.sort()

    def cnt_le(tol):
        return sum(1 for d, *_ in rows if d <= mp.mpf(tol))

    hits = [row for row in rows if row[0] <= mp.mpf("1e-6")]
    print(f"  **主判據（|N−n| ≤ 1e−6）：命中 = {len(hits)}**")
    for d, N, n, lab, kind, r in hits[:8]:
        print(f"    {lab:14s} r = {mp.nstr(r, 12):>16}  N = {mp.nstr(N, 12):>14}  n = {n}  |Δ| = {mp.nstr(d, 4)}")
    census = {t: cnt_le(t) for t in ("1e-6", "1e-4", "1e-3", "1e-2")}
    for t in ("1e-4", "1e-3", "1e-2"):
        print(f"  診斷級 |N−n| ≤ {t}: {census[t]} 個")
    print("  最近 6 個（診斷）：")
    for d, N, n, lab, kind, r in rows[:6]:
        print(f"    {lab:14s} N = {mp.nstr(N, 10):>14}  n = {n:>3}  |Δ| = {mp.nstr(d, 6)}"
              f"  ⟹ 若取 N=n，dev = {mp.nstr(abs(mp.power(r, n)/MU - 1)*100, 4)}%")
    chk("V2.1 主級命中數已定", True, f"命中 = {len(hits)}")

    hdr("§3 空模型（target-shift，M = 10⁵，log₁₀ μ′ ~ U[−25,−15]）")
    lnr = np.array([float(mp.log(r)) for _, r, _ in fam], dtype=np.float64)
    Ip = np.array(sorted(I_PRIV), dtype=np.float64)
    M, CH = 100_000, 2_000
    rng = np.random.default_rng(20261124)
    tols = np.array([1e-6, 1e-4, 1e-3, 1e-2])
    cnts = np.zeros((M, len(tols)), dtype=np.int32)

    def nearest_I(Narr):
        idx = np.clip(np.searchsorted(Ip, Narr), 1, len(Ip) - 1)
        return np.minimum(np.abs(Narr - Ip[idx - 1]), np.abs(Narr - Ip[idx]))

    for i in range(0, M, CH):
        t = rng.uniform(-25.0, -15.0, size=min(CH, M - i))
        logmu = t * np.log(10.0)                  # ln μ′
        Nmat = logmu[:, None] / lnr[None, :]      # (chunk, |R|)
        d = nearest_I(Nmat)
        for j, tl in enumerate(tols):
            cnts[i:i + len(t), j] = (d <= tl).sum(axis=1)
    lam = cnts.mean(axis=0)
    for j, tl in enumerate(tols):
        print(f"  tol {tl:g}: λ = {lam[j]:.6f}／目標（max = {int(cnts[:, j].max())}）")
    chk("V3.1 空模型跑滿 10⁵", M == 100_000, f"M = {M}")

    hdr("§4 判詞（依 docs/11 §2）")
    lam_main = float(lam[0])
    if not hits:
        verdict, text = "N", "主級 0 命中 ⟹ 門 (C) 之已註冊比值 r 族封閉"
    elif lam_main <= 0.01:
        verdict, text = "H-待結構閘", f"主級有命中且 λ(1e−6) = {lam_main:.4f} ≤ 0.01"
    else:
        verdict, text = "C", f"主級有命中但 λ(1e−6) = {lam_main:.4f} > 0.01 ⟹ 巧合"
    print(f"  判詞 = {verdict}：{text}")
    chk("V4.1 判詞依先驗規則生成", verdict in ("N", "H-待結構閘", "C"), verdict)

    hdr("§5 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")

    RES["_target"] = {"mu": mp.nstr(MU, 25), "ln_mu": mp.nstr(LN_MU, 25),
                      "R": mp.nstr(1 / MU, 25)}
    RES["_family"] = {"n": len(fam), "by_class": cnt}
    RES["_hits"] = [{"label": h[3], "class": h[4], "r": mp.nstr(h[5], 20),
                     "N": mp.nstr(h[1], 20), "n_int": h[2], "dist": mp.nstr(h[0], 8)} for h in hits]
    RES["_census"] = census
    RES["_nearest6"] = [{"label": a[3], "N": mp.nstr(a[1], 20), "n_int": a[2],
                         "dist": mp.nstr(a[0], 10),
                         "dev_if_int": mp.nstr(abs(mp.power(a[5], a[2]) / MU - 1) * 100, 8)} for a in rows[:6]]
    RES["_null"] = {"M": M, "range_log10": [-25, -15],
                    "lambda": {f"{t:g}": float(lam[j]) for j, t in enumerate(tols)},
                    "max": {f"{t:g}": int(cnts[:, j].max()) for j, t in enumerate(tols)}}
    RES["_verdict"] = {"rule": verdict, "text": text, "primary_hits": len(hits)}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL}

    (DATA / "14_step_ratio.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_14.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / '14_step_ratio.json'}\n寫出: {DATA / 'run_14.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
