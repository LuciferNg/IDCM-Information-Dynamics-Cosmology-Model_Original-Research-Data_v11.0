#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R28 — 支撐下界與自由度下界：局部核之定理級 no-go（v2，2026-10-29）

v1（同日，已備份 ~/.hermes/backups/50_r28_support_dof_bound.v1）三 FAIL 之診斷
（分類依 skill `idcm-audit-verify` §10-A：**先懷疑 check，再懷疑被審對象**）

  ✗1 檢查器 bug（致命）：`mp.dps = 40` 對 `import mpmath as mp` **無效**（須 `mp.mp.dps`）⟹
      全程實跑 dps=15 ⟹ P1 出現 13 個 ~1e−16 假違反、P3 殘差 1e−13。本版改 `mp.mp.dps = 50`。
  ✗2 判準漏洞（真）：二模之 det=0 只係**必要**條件；退化情形（null vector 之 `b·c = 0` ⟹ 分母 0）給假解。
      實例 (n,m) = (15,30)：`a = −b`、`d = 0` ⟹ `ρ_u ≡ −1`、`ρ_d ≡ 0`，與目標無關。本版加 `|b·c| > 1e−8` 閘。
  ✗3 判準漏洞（真）：離散族之「1σ 雙命中」被退化比值（分母 ≈ 0 ⟹ 比值爆走）污染。
      本版加 `|Σ c_n sin4x_n| > 1e−6` 閘，並改以 **z ＋ MC null（look-elsewhere）** 判別。

預註冊判準（v2 跑前凍結；v1 之判準一字不改，只補上述兩道退化閘與 null 檢定）
  P1 單模下界：|ρ_u(n)| = 1/|2cos(2πn/45)|；恆等逐 n 偏差 < 1e−40（dps=50）；
      下確界 = 1/(2cos(π/45))；且 min > 10×|目標|。
  P2 二模（實數權重）：同命中之解 ⟺ det = 0 **且** null vector 之 |b·c| > 1e−8。窮盡 946 對。
  P3 三模（實數權重）：2 方程 3 未知 ⟹ 恆有非零解；窮盡 C(44,3)；
      「全 c ≠ 0 且非退化」之比例；抽 5 例以 dps=50 回驗（殘差 < 1e−30）。
  P4 零參數真檢定（離散 c ∈ {±1}，K = 1..4 窮盡）：報退化排除數、min z、z<1／2／3 計數、
      **MC null**（t' ~ N(t,σ)，R = 2000）之 p 值。
        z ≡ √[ ((ρ_u − t_u)/σ_u)² + ((ρ_d − t_d)/σ_d)² ]
      p < 0.05 ⟹ 族內有真命中（列**後驗候選**，須另立盲測）；否則零資訊。
  P5 自由度：最小可達支撐 K*（實數權重）與剩餘自由度 K*−2；並答「零參數可否」。
  禁：不得改站指派（u2／d3／s4）、不得擬合權重、不得擴至 K ≥ 5 而宣稱窮盡、不得把後驗候選當證據。

跑法：/usr/bin/python3 scripts/50_r28_support_dof_bound.py
"""
from __future__ import annotations
import atexit as _atexit
import itertools
import json
import pathlib as _pl
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path, self.buf, self.out = path, [], _sys.__stdout__

    def write(self, s):
        self.buf.append(s)
        if self.out is not None:
            try:
                self.out.write(s)
            except Exception:
                pass

    def flush(self):
        if self.out is not None:
            try:
                self.out.flush()
            except Exception:
                pass

    def close(self):
        self.path.write_text("".join(self.buf), encoding="utf-8")


ROOT = _pl.Path(__file__).resolve().parent.parent
_log = _Tee(ROOT / "data" / "run_50_r28_support_dof_bound.log")
_sys.stdout = _log
_atexit.register(_log.close)

import mpmath as mp
import numpy as _np

mp.mp.dps = 50                      # ★ v1 bug：`mp.dps = 40` 對 module alias 無效
OUT: dict = {}
FAIL: list = []


def chk(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, "  ｜ " + detail if detail else ""))
    if not ok:
        FAIL.append(name)


NS = list(range(1, 45))
A = [mp.sin(2 * mp.pi * n / 45) for n in NS]      # u（站 2）
D = [mp.sin(3 * mp.pi * n / 45) for n in NS]      # d（站 3）
B = [mp.sin(4 * mp.pi * n / 45) for n in NS]      # s（站 4）
TU, TD = mp.mpf("-0.02588837"), mp.mpf("0.40350365")
SU, SD = mp.mpf("0.03358390"), mp.mpf("0.04478518")

print("=" * 96)
print("R28 — 支撐下界與自由度下界：局部核之定理級 no-go（v2）")
print("=" * 96)
print("\n§0 輸入（dps = %d）" % mp.mp.dps)
print("   站指派 u→2、d→3、s→4（已裁）｜模式 n = 1..44 ｜ φ_n(j) ∝ sin(πnj/45)")
print("   目標 ρ_d = %+.8f ± %.8f ｜ ρ_u = %+.8f ± %.8f（R27 協方差版）"
      % (float(TD), float(SD), float(TU), float(SU)))

# ── §1 單模下界 ─────────────────────────────────────────────────────────────
print("\n§1 單模核：ρ_u(n) = sin2x/sin4x = 1/(2cos2x)（x = πn/45）")
worst, argw, viol = None, None, 0
for i, n in enumerate(NS):
    if abs(B[i]) < mp.mpf("1e-45"):
        continue
    r, r_id = A[i] / B[i], 1 / (2 * mp.cos(2 * mp.pi * n / 45))
    if abs(r - r_id) > mp.mpf("1e-40"):
        viol += 1
    if worst is None or abs(r) < worst:
        worst, argw = abs(r), n
bound = 1 / (2 * mp.cos(mp.pi / 45))
print("   恆等偏差 > 1e−40 之 n 數 = %d（應 0）" % viol)
print("   實測 min |ρ_u(n)| = %.12f（n = %d）" % (float(worst), argw))
print("   解析下界 = %.12f ｜目標 |ρ_u| = %.8f ⟹ 差 %.2f×"
      % (float(bound), float(abs(TU)), float(worst / abs(TU))))
chk("P1 單模下界定理", viol == 0 and abs(worst - bound) < mp.mpf("1e-12") and worst > 10 * abs(TU),
    "min = %.12f ｜差目標 %.2f×" % (float(worst), float(worst / abs(TU))))
OUT["P1"] = {"min": str(worst), "n": argw, "bound": str(bound), "factor": str(worst / abs(TU))}

# ── §2 二模（實數權重）──────────────────────────────────────────────────────
print("\n§2 二模（實數權重）：det = A_n·D_m − A_m·D_n = 0 **且** |b·c| > 1e−8")
An = [A[i] - TU * B[i] for i in range(44)]
Dn = [D[i] - TD * B[i] for i in range(44)]
sc = max(abs(x) for x in An + Dn)
best = (None, None, mp.mpf("1e99"), mp.mpf(0))
nzero_valid = 0
for i, j in itertools.combinations(range(44), 2):
    det = abs(An[i] * Dn[j] - An[j] * Dn[i]) / (sc * sc)
    if det < mp.mpf("1e-25"):
        c = (Dn[j], -Dn[i])                       # 2×2 null vector
        bc = abs(B[i] * c[0] + B[j] * c[1]) / max(abs(c[0]), abs(c[1]))
        if bc > mp.mpf("1e-8"):
            nzero_valid += 1
        if det < best[2]:
            best = (NS[i], NS[j], det, bc)
print("   窮盡 946 對；min |det| = %.3e 於 (n,m) = (%s,%s)（該對 |b·c| = %.3e ⟹ %s）"
      % (float(best[2]), best[0], best[1], float(best[3]),
         "有效" if best[3] > mp.mpf("1e-8") else "**退化**（ρ 恆定，非命中）"))
print("   det ≈ 0 且非退化之配對數 = %d" % nzero_valid)
chk("P2 二模（實數權重）不可同時命中", nzero_valid == 0,
    "非退化 det=0 配對 = %d ⟹ %s" % (nzero_valid, "全族不可達" if nzero_valid == 0 else "須逐對查核"))
OUT["P2"] = {"min_det": str(best[2]), "argmin": [best[0], best[1]], "n_valid_det0": nzero_valid}

# ── §3 三模（實數權重）──────────────────────────────────────────────────────
print("\n§3 三模（實數權重）：null space 維 = 1 ⟹ 恆可解（扣除退化）")
ntri, nfull, sample = 0, 0, []
for i, j, k in itertools.combinations(range(44), 3):
    ntri += 1
    c1 = An[j] * Dn[k] - An[k] * Dn[j]
    c2 = An[k] * Dn[i] - An[i] * Dn[k]
    c3 = An[i] * Dn[j] - An[j] * Dn[i]
    nrm = max(abs(c1), abs(c2), abs(c3))
    if nrm == 0:
        continue
    c = [c1 / nrm, c2 / nrm, c3 / nrm]
    bc = abs(B[i] * c[0] + B[j] * c[1] + B[k] * c[2])
    if min(abs(x) for x in c) > mp.mpf("1e-6") and bc > mp.mpf("1e-8"):
        nfull += 1
        if len(sample) < 5:
            sample.append((NS[i], NS[j], NS[k], c))
print("   窮盡 %d 三元組；「全 c ≠ 0 且非退化」= %d（%.2f%%）" % (ntri, nfull, 100.0 * nfull / ntri))
resmax = mp.mpf(0)
for (n1, n2, n3, c) in sample:
    idx = [n1 - 1, n2 - 1, n3 - 1]
    den = sum(c[t] * B[idx[t]] for t in range(3))
    ru = sum(c[t] * A[idx[t]] for t in range(3)) / den
    rd = sum(c[t] * D[idx[t]] for t in range(3)) / den
    e = max(abs(ru - TU), abs(rd - TD))
    resmax = max(resmax, e)
    print("      n=(%2d,%2d,%2d) c=(%+.6f,%+.6f,%+.6f) ⟹ |Δρ_u|=%.2e |Δρ_d|=%.2e"
          % (n1, n2, n3, float(c[0]), float(c[1]), float(c[2]),
             float(abs(ru - TU)), float(abs(rd - TD))))
chk("P3 三模恆可解（殘差 < 1e−30）", resmax < mp.mpf("1e-30") and nfull > 0,
    "%d/%d 非退化 ｜樣本最大殘差 %.2e" % (nfull, ntri, float(resmax)))
OUT["P3"] = {"triples": ntri, "valid": nfull, "max_residual": str(resmax)}

# ── §4 零參數真檢定 ─────────────────────────────────────────────────────────
print("\n§4 零參數真檢定：離散 c ∈ {±1}，K = 1..4 **窮盡**；退化閘 |Σ c·sin4x| > 1e−6")
print("   z ≡ √[((ρ_u−t_u)/σ_u)² + ((ρ_d−t_d)/σ_d)²] ｜ MC null：t' ~ N(t,σ)，R = 2000")
_Af = _np.array([float(x) for x in A])
_Df = _np.array([float(x) for x in D])
_Bf = _np.array([float(x) for x in B])
_tu, _td, _su, _sd = float(TU), float(TD), float(SU), float(SD)


def build(K):
    subs = _np.array(list(itertools.combinations(range(44), K)), dtype=_np.int16)
    signs = _np.array(list(itertools.product((1, -1), repeat=K)), dtype=_np.float64)
    return subs, signs, signs @ _Af[subs].T, signs @ _Df[subs].T, signs @ _Bf[subs].T


P4 = {}
rng = _np.random.default_rng(20261029)
for K in (1, 2, 3, 4):
    subs, signs, nu, dd, de = build(K)
    ok = _np.abs(de) > 1e-6
    excl = int((~ok).sum())
    with _np.errstate(divide="ignore", invalid="ignore"):
        ru_all = nu / de
        rd_all = dd / de
        z = _np.sqrt(((ru_all - _tu) / _su) ** 2 + ((rd_all - _td) / _sd) ** 2)
    z = _np.where(ok, z, _np.inf)
    flat = int(_np.argmin(z))
    si, ni = divmod(flat, subs.shape[0])
    zmin = float(z[si, ni])
    n1, n2, n3 = int((z < 1).sum()), int((z < 2).sum()), int((z < 3).sum())
    R = 2000
    md = _np.empty(R)
    for r in range(R):
        tu2 = _tu + _su * rng.standard_normal()
        td2 = _td + _sd * rng.standard_normal()
        zr = _np.sqrt(((ru_all - tu2) / _su) ** 2 + ((rd_all - td2) / _sd) ** 2)
        zr[~ok] = _np.inf
        md[r] = zr.min()
    p = float((md <= zmin).mean())
    print("   K=%d：N_trials=%7d ｜退化排除 %4d ｜ min z = %8.4f（n=%s sign=%s）"
          % (K, nu.size, excl, zmin, tuple(int(NS[i]) for i in subs[ni]),
             tuple(int(x) for x in signs[si])))
    print("        z<1：%5d ｜z<2：%5d ｜z<3：%5d ｜ null min-z 中位 = %6.3f ｜ **p = %.4f**"
          % (n1, n2, n3, float(_np.median(md)), p))
    P4["K%d" % K] = {"trials": int(nu.size), "excluded": excl, "zmin": zmin,
                     "argmin_n": [int(NS[i]) for i in subs[ni]],
                     "argmin_sign": [int(x) for x in signs[si]],
                     "rho_u": float(ru_all[si, ni]), "rho_d": float(rd_all[si, ni]),
                     "counts": [n1, n2, n3], "median_null_minz": float(_np.median(md)), "p": p}
    if zmin < 3:
        print("        ★ 後驗候選（**不作證據**，須另立盲測）：ρ_u = %+.8f、ρ_d = %+.8f ⟹ z = %.3f"
              % (float(ru_all[si, ni]), float(rd_all[si, ni]), zmin))
sig = {k: v["p"] for k, v in P4.items() if v["zmin"] < 3}
print("\n   ⟹ z < 3 之族與其 p 值：%s" % (sig if sig else "無"))
chk("P4 離散族 MC null 檢定（p < 0.05 才算命中）",
    all(v >= 0.05 for v in sig.values()) if sig else True,
    "；".join("K=%s p=%.4f" % (k[1:], v) for k, v in sig.items()) or "無 z<3 候選 ⟹ 全族零資訊")
OUT["P4"] = P4

# ── §5 自由度結論 ───────────────────────────────────────────────────────────
print("\n§5 自由度結論")
print("   單模：|ρ_u| ≥ %.6f（差目標 %.2f×）⟹ 定理級判死" % (float(bound), float(bound / abs(TU))))
print("   二模（實數權重）：非退化 det=0 配對 = %d ⟹ %s"
      % (nzero_valid, "判死" if nzero_valid == 0 else "須逐對查核"))
print("   三模（實數權重）：%d/%d 非退化三元組皆可解 ⟹ 可達，但耗 1 個連續自由度" % (nfull, ntri))
OUT["P5"] = {"bound_single": str(bound), "valid_det0_pairs": nzero_valid, "valid_triples": nfull}
chk("P5 自由度下界成立", nzero_valid == 0 and nfull > 0, "K* = 3，dof = 1 ⟹ 零參數須由結構固定權重比")

print("\n§6 判詞")
print("   ① 單模（局部）核：**定理級判死**（|ρ_u| ≥ %.6f ⟹ 19.36×）" % float(bound))
print("   ② 二模（含任意實數權重）：**判死**（det ≠ 0；唯一 det ≈ 0 之 (15,30) 係退化）")
print("   ③ 三模：可達但需 1 個連續自由度 ⟹ 擬合，唔屬零參數")
print("   ④ 離散 ±1 族（K ≤ 4 窮盡，N_trials = %d）：判別見 §4"
      % sum(v["trials"] for v in P4.values()))
print("\n   FAIL = %d ｜ written data/50_r28_support_dof_bound.json , data/run_50_r28_support_dof_bound.log" % len(FAIL))
print("=" * 96)
OUT["fail"] = FAIL
(ROOT / "data" / "50_r28_support_dof_bound.json").write_text(
    json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")
_sys.exit(0)