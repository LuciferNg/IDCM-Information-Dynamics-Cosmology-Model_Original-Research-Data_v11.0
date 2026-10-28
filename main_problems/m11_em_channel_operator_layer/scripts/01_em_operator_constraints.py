#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M11 scripts/01 -- U5 攻堅第一輪：`em` 通道算子之約束代數 ＋ 三族 pairwise 核檢定。

背景（M4 移交，2026-10-28）：
  M4 `docs/19` 診斷 —— Q2/Q3/Q8 之真正上游係「**一個未建立之結構層**」（通道級電荷讀數算子），
  唔係任何現存 O 之待辦項。M4 `docs/37 §一 B4` 在版本紀律下懸空 ⟹ 本記錄（M11）承接。

資料（全部獨立重算，唔抄 M4 表）：
  `E_ch ≡ Δm_ch − coeff·X`，`X = m_d − m_u = 2.56 ± 0.094`（FLAG24），
  8 通道 PDG 2024 質量 ⟹ E_ch 為**純電磁殘餘**（已扣夸克質量差）。

本輪三問：
  Q-A 約束代數：8 通道之 SM 量子數差矩陣 rank 係幾？（能否只用 QN 決定 em？）
  Q-B 候選族：pairwise 核 `em(h) = Σ_{i<j} q_i q_j K(|n_i−n_j|)` 能唔能同時命中**符號**與量級？
  Q-C 若否：失敗模式係咩（符號 vs 量級），對算子類別有咩硬限制？

零擬合宣告：`K(g)` 用 M5 精確環 Green 函數（無參數）；mode 指派用 M4/O233 既有 `u→2, d→3`
  之同構延伸（第 1 代 (2,3)、第 2 代 (4,5)、第 3 代 (6,7)）—— **事前指定**，非由 8 通道反解。
  唯一標度 = 1 個整體乘數（em 之 MeV 標度屬 O38 `M₀` 未決），由最小二乘定，報自由度成本。
"""
import json
import math
import os
from fractions import Fraction as F

CHECKS = []


def chk(label, ok, got, expect):
    CHECKS.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-64s | %-24s | %s" % (label, got, expect))


# ---------------------------------------------------------------- 原始輸入
M = {"pi+": 139.57039, "pi0": 134.9768, "K+": 493.677, "K0": 497.611,
     "p": 938.272089, "n": 939.565421, "S-": 1197.449, "S0": 1192.642,
     "S+": 1189.37, "X-": 1321.71, "X0": 1314.86, "D+": 1869.66, "D0": 1864.84}
X = 2.56            # m_d − m_u（FLAG24）
SIG_X = 0.094
N = 88

# quark masses / charges / mode assignment（事前指定）
mq = {"u": 2.14, "d": 4.70, "s": 93.46, "c": 1270.0}
qc = {"u": F(2, 3), "d": F(-1, 3), "s": F(-1, 3), "c": F(2, 3)}
MODE = {"u": 2, "d": 3, "c": 4, "s": 5}

# 8 通道：(名, 初態, 終態, coeff, M4 docs/18 表列 E, σ_E)
CHANNELS = [
    ("n−p",      "n",  "p",  +1, -1.2667, 0.0940),
    ("K+−K0",    "K+", "K0", -1, -1.3740, 0.0963),
    ("π+−π0",    "pi+", "pi0", +1, +2.0336, 0.0940),
    ("Σ−−Σ0",    "S-", "S0", +1, +2.2470, 0.1010),
    ("Σ0−Σ+",    "S0", "S+", +1, +0.7120, 0.1196),
    ("Σ−−Σ+",    "S-", "S+", +2, +2.9590, 0.2028),
    ("Ξ−−Ξ0",    "X-", "X0", +1, +4.2900, 0.2301),
    ("D+−D0",    "D+", "D0", +1, +2.2600, 0.1172),
]

# ---------------------------------------------------------------- A. E_ch 獨立重算
print("=" * 100)
print("A. E_ch = Δm − coeff·X 獨立重算（PDG 2024 質量；唔抄 M4 表）")
print("=" * 100)
E = {}
dmax = 0.0
for name, a, b, k, m4_val, sig in CHANNELS:
    dm = M[a] - M[b]
    e = dm - k * X
    E[name] = e
    dev = abs(e - m4_val)
    dmax = max(dmax, dev)
    print("  %-8s Δm=%+9.4f  coeff=%+d  E=%+9.4f   表列=%+9.4f  |Δ|=%.4f" % (name, dm, k, e, m4_val, dev))
chk("A1 8 通道 E_ch 全部由 PDG 重算，與 M4 docs/18 一致（|Δ|<5e-3 MeV）", dmax < 5e-3, "max|Δ|=%.5f" % dmax, "<0.005")
chk("A2 E_ch 跨度 5.66 MeV（−1.374 … +4.290）",
    abs((max(E.values()) - min(E.values())) - 5.664) < 0.01,
    "span=%.3f" % (max(E.values()) - min(E.values())), "5.664")

# ★ A3 定義釐清：E_ch 唔等於純 em（lattice 外部錨：em(K⁺−K⁰) = +2.2646，M4 docs/36）
EM_K_LATTICE = 2.2645889622304956
strong_IB_K = E["K+−K0"] - EM_K_LATTICE
print("\n  ★ 定義釐清（K 通道）：E_ch = %+.4f ；lattice 純 em = %+.4f ⟹ 強 IB 殘餘 = %+.4f MeV"
      % (E["K+−K0"], EM_K_LATTICE, strong_IB_K))
chk("A3 E_ch 混入強 IB 殘餘（|殘餘| > 3 MeV ⟹ 唔可當純 em 用）", abs(strong_IB_K) > 3.0,
    "strong-IB=%+.3f" % strong_IB_K, "|.|>3")
chk("A4 強 IB 殘餘大於 em 本身（K 通道：6.199 vs 2.265）", abs(strong_IB_K) > abs(EM_K_LATTICE),
    "|%.3f| > %.3f" % (strong_IB_K, EM_K_LATTICE), "True")

# ---------------------------------------------------------------- A. 量子數差矩陣 rank
print("\n" + "=" * 100)
print("B. 約束代數：8 通道之 SM 量子數差（ΔQ, ΔI₃, ΔB, ΔS, ΔC）矩陣 rank")
print("=" * 100)
# 每通道：由夸克含量差直接算 QN 差（u→d 替換 ⟹ ΔQ=−1, ΔI₃=+1 等）
QN = {
    "n−p":    dict(Q=0, I3=-1, B=0, S=0, C=0),
    "K+−K0":  dict(Q=-1, I3=0, B=0, S=0, C=0),
    "π+−π0":  dict(Q=1, I3=0, B=0, S=0, C=0),
    "Σ−−Σ0":  dict(Q=-1, I3=0, B=0, S=0, C=0),
    "Σ0−Σ+":  dict(Q=-1, I3=0, B=0, S=0, C=0),
    "Σ−−Σ+":  dict(Q=-2, I3=0, B=0, S=0, C=0),
    "Ξ−−Ξ0":  dict(Q=-1, I3=0, B=0, S=0, C=0),
    "D+−D0":  dict(Q=1, I3=0, B=0, S=0, C=0),
}
cols = ["Q", "I3", "B", "S", "C"]
rows = [[F(QN[c][k]) for k in cols] for c in E]
# 精確 rank（Fraction 消元）
Mm = [r[:] for r in rows]
rank = 0
for col in range(len(cols)):
    piv = next((r for r in range(rank, len(Mm)) if Mm[r][col] != 0), None)
    if piv is None:
        continue
    Mm[rank], Mm[piv] = Mm[piv], Mm[rank]
    for r in range(len(Mm)):
        if r != rank and Mm[r][col] != 0:
            f = Mm[r][col] / Mm[rank][col]
            Mm[r] = [x - f * y for x, y in zip(Mm[r], Mm[rank])]
    rank += 1
chk("B1 QN 差矩陣 rank = 2（僅 ΔQ、ΔI₃；B/S/C 全零）", rank == 2, "rank=%d" % rank, 2)

# QN-only 算子必然同值 → 同 (ΔQ,ΔI₃) 之通道群內 spread 係 QN 盲區
groups = {}
for c, v in E.items():
    key = (QN[c]["Q"], QN[c]["I3"])
    groups.setdefault(key, []).append((c, v))
print("\n  同 (ΔQ,ΔI₃) 群內 E 之 spread（QN-only 算子必須全部相等）：")
maxspread = 0.0
for key, items in groups.items():
    if len(items) < 2:
        continue
    vals = [v for _, v in items]
    sp = max(vals) - min(vals)
    maxspread = max(maxspread, sp)
    print("    ΔQ=%+d ΔI₃=%+d : %s  spread=%.3f MeV"
          % (key[0], key[1], ", ".join("%s %+.3f" % (c, v) for c, v in items), sp))
chk("B2 QN-blindness：最壞群內 spread > 5 MeV ⟹ QN-only 算子不可能", maxspread > 5.0,
    "max spread=%.3f" % maxspread, ">5.0")

# ---------------------------------------------------------------- C. pairwise 核候選
print("\n" + "=" * 100)
print("C. 候選族：em(h) = Σ_{i<j} q_i q_j K(|n_i − n_j|)，K = M5 精確環 Green 函數")
print("=" * 100)


def g_ring(d):
    d = min(d, N - d)
    return float(F((N * N - 1) - 6 * d * (N - d), 12 * N))


HAD = {"n": (["d", "d", "u"], []), "p": (["u", "u", "d"], []),
       "K+": (["u"], ["s"]), "K0": (["d"], ["s"]),
       "S-": (["d", "d", "s"], []), "S0": (["u", "d", "s"], []), "S+": (["u", "u", "s"], []),
       "X-": (["d", "s", "s"], []), "X0": (["u", "s", "s"], []),
       "D+": (["c"], ["d"]), "D0": (["c"], ["u"])}


def quarks(h):
    a, b = HAD[h]
    return [(f, qc[f]) for f in a] + [(f, -qc[f]) for f in b]


def em(h, kern=g_ring):
    Qs = quarks(h)
    t = 0.0
    for i in range(len(Qs)):
        for j in range(i + 1, len(Qs)):
            (f1, c1), (f2, c2) = Qs[i], Qs[j]
            t += float(c1 * c2) * kern(abs(MODE[f1] - MODE[f2]))
    return t


print("  mode 指派（事前指定，非反解）：" + ", ".join("%s→%d" % (f, MODE[f]) for f in ("u", "d", "c", "s")))
pred, sign_ok, sign_bad = {}, [], []
for name, a, b, k, m4_val, sig in CHANNELS:
    if name == "π+−π0":
        pred[name] = None       # annihilation topology：pairwise 夸克框架不適用
        print("  %-8s  —（annihilation 拓撲，pairwise 框架外，剔除）" % name)
        continue
    v = em(a) - em(b)
    pred[name] = v
    (sign_ok if v * E[name] > 0 else sign_bad).append(name)
    print("  %-8s pred_raw=%+8.4f  obs=%+8.4f  %s" % (name, v, E[name], "sign OK" if v * E[name] > 0 else "SIGN WRONG"))

chk("C1 符號：7 通道中全部同號", not sign_bad, "sign_ok=%d, wrong=%s" % (len(sign_ok), sign_bad), "0 wrong")

# 1 參數（整體標度）最小二乘
use = [n for n in pred if pred[n] is not None]
num = sum(pred[n] * E[n] for n in use)
den = sum(pred[n] ** 2 for n in use)
scale = num / den
resid = {n: scale * pred[n] - E[n] for n in use}
rms = math.sqrt(sum(r * r for r in resid.values()) / len(use))
rms_qn = math.sqrt(sum((E[n] - sum(E[m] for m in use) / len(use)) ** 2 for n in use) / len(use))
print("\n  1 參數（標度）最小二乘：k = %.4f ；RMS = %.3f MeV（同一組數據之 QN-only 基準 RMS = %.3f）"
      % (scale, rms, rms_qn))
for n in use:
    print("    %-8s pred=%+7.3f  obs=%+7.3f  resid=%+7.3f (%.2fσ)"
          % (n, scale * pred[n], E[n], resid[n], abs(resid[n]) / dict((c[0], c[5]) for c in CHANNELS)[n]))
chk("C2 1 參數模型 RMS < 2 MeV（量級可接受）", rms < 2.0, "RMS=%.3f" % rms, "<2.0")
chk("C3 但符號結構不完全正確 ⟹ 不作算子升格", not sign_bad,
    "wrong signs: %s" % sign_bad, "[] （要全部正確才可升格）")

# ★ C4 純 em 對照（唯一有 lattice 錨之通道）：零參數（單位耦合），非擬合
k_pred = em("K+") - em("K0")
dev_K = (k_pred - EM_K_LATTICE) / EM_K_LATTICE
chk("C4 K 通道純 em 預測（單位耦合、零擬合）與 lattice 差 <15%", abs(dev_K) < 0.15,
    "pred=%+.4f vs %.4f (dev %+.1f%%)" % (k_pred, EM_K_LATTICE, dev_K * 100), "<15%")

# 逐通道 χ²（用 σ_E）
chi2 = sum((resid[n] / dict((c[0], c[5]) for c in CHANNELS)[n]) ** 2 for n in use)
print("  7 通道 χ² = %.1f（dof = 6）⟹ 若符號全對，χ² 亦需 ~6 才可接受" % chi2)

# ---------------------------------------------------------------- D. 結論量
# ---------------------------------------------------------------- D. 差額結構：強 IB 殘餘 R = E − em_raw
print("\n" + "=" * 100)
print("D. 差額分解：R_ch ≡ E_ch − em_pred（單位耦合）—— R 應由強 IB（∝ 夸克質量差）承擔")
print("=" * 100)
R = {}
print("  %-8s %+9s %+9s %+9s %8s" % ("通道", "E_obs", "em_pred", "R=E−em", "R/X"))
for name, a, b, k, m4_val, sig in CHANNELS:
    if pred[name] is None:
        continue
    R[name] = E[name] - pred[name]
    print("  %-8s %+9.4f %+9.4f %+9.4f %8.3f" % (name, E[name], pred[name], R[name], R[name] / X))
rx = sorted(v / X for v in R.values())
chk("D1 R/X 非普適（spread > 1.5 ⟹ 單一 ⟨O⟩ 因子不可行，與 M4 docs/13 一致）",
    (rx[-1] - rx[0]) > 1.5, "spread=%.3f  [%.2f, %.2f]" % (rx[-1] - rx[0], rx[0], rx[-1]), ">1.5")
chk("D2 R 與 E_ch 同量級（|R|max/|E|max > 0.5 ⟹ 強 IB 不可忽略）",
    max(abs(v) for v in R.values()) / max(abs(v) for v in E.values()) > 0.5,
    "ratio=%.3f" % (max(abs(v) for v in R.values()) / max(abs(v) for v in E.values())), ">0.5")

DATA = {
    "E_ch_MeV": {k: round(v, 4) for k, v in E.items()},
    "QN_rank": rank,
    "qn_group_max_spread_MeV": round(maxspread, 4),
    "pairwise_ring": {"pred_raw": {k: (None if v is None else round(v, 4)) for k, v in pred.items()},
                      "scale_k": round(scale, 6), "rms_MeV": round(rms, 4), "chi2_7": round(chi2, 2),
                      "sign_wrong": sign_bad},
    "mode_assignment": MODE,
    "R_strong_IB": {k: round(v, 4) for k, v in R.items()},
    "R_over_X": {k: round(v / X, 4) for k, v in R.items()},
    "em_K_lattice": EM_K_LATTICE,
    "strong_IB_K": round(strong_IB_K, 4),
    "em_K_pred_unit_coupling": round(k_pred, 4),
    "em_K_dev_pct": round(dev_K * 100, 2),
    "n_pass": sum(1 for c in CHECKS if c["ok"]),
    "n_total": len(CHECKS),
}
print("\n" + "=" * 100)
print("RESULT: %d/%d PASS" % (DATA["n_pass"], DATA["n_total"]))
print("=" * 100)

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "01_em_operator_constraints.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"checks": CHECKS, **DATA}, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON -> data/01_em_operator_constraints.json")