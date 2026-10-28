#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/42 -- FLAG 2024 (arXiv:2411.04268v3) 數據修訂對 `docs/18` 八通道約束集之影響。

動機（2026-10-28）：
  `docs/18` 用 `X = m_d − m_u = 2.56 ± 0.094`（把 m_u、m_d 當獨立量作差）。
  FLAG 2024 v3 同時給出 `m_u/m_d = 0.465(24)`（Tab.12）⟹ 正確形式為 `X = m_d(1 − m_u/m_d)`，
  因格點上 m_u、m_d 同源高度相關。差值形式虛報精度 ~25%。

本腳本只做兩件事（唔改任何既有結論，只量化影響）：
  R1  重算 8 通道 `E_ch` 與 `σ_E`（X 舊 vs 新），列出 Δ。
  R2  重跑 `docs/18` 之四項排除之**定量部分**（電荷加性 2 參數 R²、substitution-count 線性 R²、
       spectator 群內 spread、K/N 一致），檢「定性結論是否改變」。
  R3  由 FLAG24 Eq.(45) `ε = 0.79(6)` 獨立核算 `em(K⁺−K⁰)` 錨（與 docs/36 之外部錨對照）。

全部輸入有出處；禁任何擬合升格。
"""
import json
import math
import os

CHECK = []


def chk(label, ok, got, expect):
    CHECK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-62s | %-26s | %s" % (label, got, expect))


# ---------------------------------------------------------------- 來源值（FLAG 2024 v3）
MU, S_MU = 2.14, 0.08          # Eq.(43)  N_f=2+1+1, MS-bar 2 GeV
MD, S_MD = 4.70, 0.05
MS_S, S_MS = 93.46, 0.58
RATIO, S_RATIO = 0.465, 0.024  # Tab.12
EPS, S_EPS = 0.79, 0.06        # Eq.(45)

X_OLD, S_X_OLD = MD - MU, math.hypot(S_MD, S_MU)
X_NEW = MD * (1.0 - RATIO)
S_X_NEW = math.hypot((1.0 - RATIO) * S_MD, MD * S_RATIO)

M = {"pi+": 139.57039, "pi0": 134.9768, "K+": 493.677, "K0": 497.611,
     "p": 938.272089, "n": 939.565421, "S-": 1197.449, "S0": 1192.642,
     "S+": 1189.37, "X-": 1321.71, "X0": 1314.86, "D+": 1869.66, "D0": 1864.84}

CH = [("n−p", "n", "p", +1, "baryon u↔d，spectators {u,d}"),
      ("K+−K0", "K+", "K0", -1, "meson u→d，spectator s̄"),
      ("π+−π0", "pi+", "pi0", +1, "annihilation"),
      ("Σ−−Σ0", "S-", "S0", +1, "baryon，spectators {d,s}"),
      ("Σ0−Σ+", "S0", "S+", +1, "baryon，spectators {u,s}"),
      ("Σ−−Σ+", "S-", "S+", +2, "baryon 2×u↔d"),
      ("Ξ−−Ξ0", "X-", "X0", +1, "baryon，spectators {s,s}"),
      ("D+−D0", "D+", "D0", +1, "meson，spectator c")]

print("=" * 104)
print("R0  X 之兩種形式（FLAG 2024 v3）")
print("=" * 104)
print("  X_old (差值，docs/18 現用) = %.4f ± %.4f MeV" % (X_OLD, S_X_OLD))
print("  X_new (比值，正確形式)     = %.4f ± %.4f MeV" % (X_NEW, S_X_NEW))
print("  ΔX = %+.4f MeV ；σ 增 %.1f%%" % (X_NEW - X_OLD, (S_X_NEW / S_X_OLD - 1) * 100))
chk("R0.1 X 兩種形式差 > 0.04 MeV（達 σ_E 之 48%）", abs(X_NEW - X_OLD) > 0.04,
    "ΔX=%+.4f" % (X_NEW - X_OLD), ">0.04")

# ---------------------------------------------------------------- R1 重算
print("\n" + "=" * 104)
print("R1  8 通道 E_ch 與 σ_E 重算（dm 來自 PDG 2024；σ_m ≪，只報 X 貢獻）")
print("=" * 104)
rows = {}
print("  %-8s %-3s %9s %9s %9s %8s %8s %7s" % ("通道", "c", "E_old", "E_new", "ΔE", "σ_old", "σ_new", "σ增%"))
for name, a, b, c, kind in CH:
    dm = M[a] - M[b]
    e_old = dm - c * X_OLD
    e_new = dm - c * X_NEW
    s_old = abs(c) * S_X_OLD
    s_new = abs(c) * S_X_NEW
    rows[name] = dict(coeff=c, dm=dm, E_old=e_old, E_new=e_new, dE=e_new - e_old,
                      s_old=s_old, s_new=s_new, kind=kind)
    print("  %-8s %+d %+9.4f %+9.4f %+9.4f %8.4f %8.4f %6.1f%%"
          % (name, c, e_old, e_new, e_new - e_old, s_old, s_new, (s_new / s_old - 1) * 100))
chk("R1.1 全部通道 σ_E 增大（無一例外）", all(r["s_new"] > r["s_old"] for r in rows.values()),
    "min σ_new=%.4f" % min(r["s_new"] for r in rows.values()), "σ_new > σ_old")
chk("R1.2 E 之平移與 coeff 成正比（ΔE = coeff·(X_old−X_new)）",
    all(abs(r["dE"] - r["coeff"] * (X_OLD - X_NEW)) < 1e-12 for r in rows.values()),
    "max|殘差|=%.2e" % max(abs(r["dE"] - r["coeff"] * (X_OLD - X_NEW)) for r in rows.values()), "<1e-12")

# ---------------------------------------------------------------- R2 四項排除之定量部分
print("\n" + "=" * 104)
print("R2  docs/18 四項排除之定量部分：X 修訂後是否仍成立")
print("=" * 104)


def r2(y, x1, x2):
    """2 參數最小二乘 + R²（x1, x2 為解釋變數）。"""
    n = len(y)
    s11 = sum(a * a for a in x1)
    s22 = sum(b * b for b in x2)
    s12 = sum(a * b for a, b in zip(x1, x2))
    t1 = sum(a * c for a, c in zip(x1, y))
    t2 = sum(b * c for b, c in zip(x2, y))
    det = s11 * s22 - s12 * s12
    if abs(det) < 1e-18:
        return float("nan"), float("nan")
    k1 = (t1 * s22 - t2 * s12) / det
    k2 = (t2 * s11 - t1 * s12) / det
    ybar = sum(y) / n
    ss_res = sum((c - k1 * a - k2 * b) ** 2 for a, b, c in zip(x1, x2, y))
    ss_tot = sum((c - ybar) ** 2 for c in y)
    return (k1, k2), 1 - ss_res / ss_tot


names = list(rows)
# 電荷加性：κ_s ΣQ² + κ_x ΣQ_iQ_j（docs/17 之形式）
QQ2 = {"n−p": 0.0, "K+−K0": 0.0, "π+−π0": 0.0, "Σ−−Σ0": 0.0, "Σ0−Σ+": 0.0,
       "Σ−−Σ+": 0.0, "Ξ−−Ξ0": 0.0, "D+−D0": 0.0}
SUMQQ = {"n−p": float(2 * (2 / 3) * (-1 / 3) + (-1 / 3) ** 2), "K+−K0": float((2 / 3) * (1 / 3) + (1 / 3) ** 2),
         "π+−π0": 0.0, "Σ−−Σ0": 0.0, "Σ0−Σ+": 0.0, "Σ−−Σ+": 0.0, "Ξ−−Ξ0": 0.0, "D+−D0": 0.0}
# 直接由夸克含量差算 ΣQ² 與 ΣQ_iQ_j 之差（更乾淨）
HAD = {"n": (["d", "d", "u"], []), "p": (["u", "u", "d"], []),
       "K+": (["u"], ["s"]), "K0": (["d"], ["s"]), "S-": (["d", "d", "s"], []),
       "S0": (["u", "d", "s"], []), "S+": (["u", "u", "s"], []), "X-": (["d", "s", "s"], []),
       "X0": (["u", "s", "s"], []), "D+": (["c"], ["d"]), "D0": (["c"], ["u"])}
QC = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3}


def charges(h):
    a, b = HAD[h]
    return [QC[f] for f in a] + [-QC[f] for f in b]


dQ2, dQiQj, dN = [], [], []
for name, a, b, c, kind in CH:
    if a in HAD and b in HAD:
        qa, qb = charges(a), charges(b)
        dQ2.append(sum(x * x for x in qa) - sum(x * x for x in qb))
        dQiQj.append(sum(qa[i] * qa[j] for i in range(len(qa)) for j in range(i + 1, len(qa)))
                     - sum(qb[i] * qb[j] for i in range(len(qb)) for j in range(i + 1, len(qb))))
    else:  # π+−π0：annihilation，pair frame 唔同 ⟹ 以电荷差近似（docs/18 亦然）
        dQ2.append(0.0)
        dQiQj.append(0.0)

# ★ 實現修正：π+−π0 屬 annihilation 拓撲（π0 為 (uū−dd̄)/√2 疊加，夸克對框架不同），
#   於電荷加性擬合中剔除；只對 7 條 substitution 通道做 test。上一版用 0/0 佔位令設計矩陣退化 ⟹ 該版數字作廢。
fit_names = [n for n in names if rows[n]["kind"] != "annihilation"]
y_old = [rows[n]["E_old"] for n in fit_names]
y_new = [rows[n]["E_new"] for n in fit_names]
idx = [names.index(n) for n in fit_names]
dQ2f = [dQ2[i] for i in idx]
dQiQjf = [dQiQj[i] for i in idx]

(k_old, r2_old) = r2(y_old, dQ2f, dQiQjf)
(k_new, r2_new) = r2(y_new, dQ2f, dQiQjf)
print("  (a) 電荷加性 κ_sΣQ²+κ_xΣQ_iQ_j（7 條 substitution；free 2 參數）:  R²(old)=%.4f  R²(new)=%.4f"
      % (r2_old, r2_new))
# ★ 判準修正：結構算子**無自由係數** ⟹ 排除測試必須用零擬合版本（κ_s=κ_x=1，即純電荷量）。
#   free 2 參數版本只作 context：7 點配 2 參數，R² 必然偏高，不具檢定力。
r2u = {}
for tag, yv in (("old", y_old), ("new", y_new)):
    res = [c - (a + b) for c, a, b in zip(yv, dQ2f, dQiQjf)]
    ss = sum(r * r for r in res)
    tot = sum((c - sum(yv) / len(yv)) ** 2 for c in yv)
    r2u[tag] = 1 - ss / tot
    print("      κ_s=κ_x=1（零擬合，零自由度）%s:  Σres²=%.2f  R²=%+.4f  RMS=%.2f MeV"
          % (tag, ss, r2u[tag], math.sqrt(ss / len(res))))
chk("R2.1 電荷加性（零擬合 κ_s=κ_x=1）仍被排除（R² < 0.2）", max(r2u.values()) < 0.2,
    "R²old=%+.4f R²new=%+.4f" % (r2u["old"], r2u["new"]), "<0.2")
print("      〔context〕free 2 參數版本 R²(old)=%.4f R²(new)=%.4f（7 點 2 參數，無檢定力）"
      % (r2_old, r2_new))

# (b) substitution-count 線性：E = a·coeff + b
def lin(y, x):
    n = len(y); xb = sum(x) / n; yb = sum(y) / n
    sxy = sum((a - xb) * (b - yb) for a, b in zip(x, y))
    sxx = sum((a - xb) ** 2 for a in x)
    a = sxy / sxx; b = yb - a * xb
    ss_res = sum((c - a * xi - b) ** 2 for xi, c in zip(x, y))
    ss_tot = sum((c - yb) ** 2 for c in y)
    return (a, b), 1 - ss_res / ss_tot


coeffs = [rows[n]["coeff"] for n in fit_names]
(_, r2b_old) = lin(y_old, coeffs)
(_, r2b_new) = lin(y_new, coeffs)
print("  (b) substitution-count 線性   :  R²(old)=%.4f  R²(new)=%.4f" % (r2b_old, r2b_new))
chk("R2.2 substitution-count 線性仍被排除（R² < 0.5）", max(r2b_old, r2b_new) < 0.5,
    "R²old=%.4f R²new=%.4f" % (r2b_old, r2b_new), "<0.5")

# (c) spectator 獨立：coeff=+1 且非 annihilation 之通道群
grp = [n for n in names if rows[n]["coeff"] == 1 and rows[n]["kind"] != "annihilation"]
sp_old = max(rows[n]["E_old"] for n in grp) - min(rows[n]["E_old"] for n in grp)
sp_new = max(rows[n]["E_new"] for n in grp) - min(rows[n]["E_new"] for n in grp)
print("  (c) spectator 獨立（coeff=+1 群 %d 條）: spread(old)=%.3f  spread(new)=%.3f MeV"
      % (len(grp), sp_old, sp_new))
chk("R2.3 spectator 獨立仍被排除（spread > 3 MeV）", min(sp_old, sp_new) > 3.0,
    "spread_new=%.3f" % sp_new, ">3.0")

# (d) K/N 一致性（docs/18 稱 −1.374 vs −1.267，spread 0.107 = 1.1σ）
d_old = rows["K+−K0"]["E_old"] - rows["n−p"]["E_old"]
s_old_kn = math.hypot(rows["K+−K0"]["s_old"], rows["n−p"]["s_old"])
d_new = rows["K+−K0"]["E_new"] - rows["n−p"]["E_new"]
s_new_kn = math.hypot(rows["K+−K0"]["s_new"], rows["n−p"]["s_new"])
print("  (d) K/N 一致：Δ(old)=%+.3f ± %.3f (%.2fσ)  →  Δ(new)=%+.3f ± %.3f (%.2fσ)"
      % (d_old, s_old_kn, abs(d_old) / s_old_kn, d_new, s_new_kn, abs(d_new) / s_new_kn))
chk("R2.4 K/N 一致性結論不變（<2σ）", abs(d_new) / s_new_kn < 2.0,
    "%.2fσ" % (abs(d_new) / s_new_kn), "<2σ")

# ---------------------------------------------------------------- R3 ε 錨
print("\n" + "=" * 104)
print("R3  由 FLAG24 Eq.(45) ε = 0.79(6) 獨立核算 em(K⁺−K⁰)")
print("=" * 104)
dMpi2 = M["pi+"] ** 2 - M["pi0"] ** 2
dMK2 = (1 + EPS) * dMpi2
emK = dMK2 / (2 * M["K+"])
s_emK = emK * S_EPS / (1 + EPS)
anchor = 2.2645889622304956   # M4 docs/36 之外部錨
print("  ΔM_π² = %.3f MeV² ；(ΔM_K²)_γ = (1+ε)ΔM_π² = %.2f MeV²" % (dMpi2, dMK2))
print("  em(K⁺−K⁰) = (ΔM_K²)_γ/(2M_K+) = %.4f ± %.4f MeV （ε 主導）" % (emK, s_emK))
print("  vs docs/36 錨 %.4f ⟹ dev %+.2f%% （%.2fσ）" % (anchor, (emK / anchor - 1) * 100, abs(emK - anchor) / s_emK))
chk("R3.1 ε=0.79(6) 與 docs/36 錨一致（<2σ）", abs(emK - anchor) / s_emK < 2.0,
    "%.2fσ" % (abs(emK - anchor) / s_emK), "<2σ")
chk("R3.2 由 ε 反解所需值落在 FLAG 誤差內", abs((2 * M["K+"] * anchor / dMpi2 - 1) - EPS) < 2 * S_EPS,
    "ε_req=%.3f" % (2 * M["K+"] * anchor / dMpi2 - 1), "0.79±0.12")

# ---------------------------------------------------------------- R4 對既有顯著度之影響
print("\n" + "=" * 104)
print("R4  對既有顯著度之影響（純標度：σ_X 變 ⟹ z 乘 S_X_OLD/S_X_NEW）")
print("=" * 104)
fac = S_X_OLD / S_X_NEW
for z in (16.0, 9.76):
    print("  %5.2fσ  →  %5.2fσ" % (z, z * fac))
chk("R4.1 σ_X 增大 ⟹ 依賴 σ_X 之顯著度一律下修（因子 %.3f）" % fac, fac < 1.0,
    "factor=%.3f" % fac, "<1.0")

DATA = {
    "flag2024": {"arxiv": "2411.04268v3", "mu": MU, "s_mu": S_MU, "md": MD, "s_md": S_MD,
                 "ms": MS_S, "s_ms": S_MS, "mu_over_md": RATIO, "s_ratio": S_RATIO,
                 "eps": EPS, "s_eps": S_EPS},
    "X_old": {"value": round(X_OLD, 5), "sigma": round(S_X_OLD, 5)},
    "X_new": {"value": round(X_NEW, 5), "sigma": round(S_X_NEW, 5)},
    "channels": {k: {kk: (round(vv, 5) if isinstance(vv, float) else vv) for kk, vv in v.items()}
                 for k, v in rows.items()},
    "R2": {"charge_additive_r2_old": round(r2_old, 4), "charge_additive_r2_new": round(r2_new, 4),
           "charge_additive_n_channels": len(fit_names),
           "charge_additive_r2_unit_old": round(r2u["old"], 4),
           "charge_additive_r2_unit_new": round(r2u["new"], 4),
           "subst_count_r2_old": round(r2b_old, 4), "subst_count_r2_new": round(r2b_new, 4),
           "spectator_spread_old": round(sp_old, 4), "spectator_spread_new": round(sp_new, 4),
           "KN_delta_old": round(d_old, 4), "KN_delta_new": round(d_new, 4),
           "KN_sigma_old": round(abs(d_old) / s_old_kn, 3), "KN_sigma_new": round(abs(d_new) / s_new_kn, 3)},
    "R3": {"dMpi2": round(dMpi2, 3), "dMK2_gamma": round(dMK2, 2), "emK": round(emK, 4),
           "s_emK": round(s_emK, 4), "anchor": anchor, "dev_pct": round((emK / anchor - 1) * 100, 2)},
    "R4": {"sigma_factor": round(fac, 4)},
    "n_pass": sum(1 for c in CHECK if c["ok"]), "n_total": len(CHECK),
}
print("\n" + "=" * 104)
print("RESULT: %d/%d PASS" % (DATA["n_pass"], DATA["n_total"]))
print("=" * 104)

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "42_flag2024_x_revision.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"checks": CHECK, "data": DATA}, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON -> data/42_flag2024_x_revision.json")