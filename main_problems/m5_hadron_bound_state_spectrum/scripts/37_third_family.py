"""M5 docs/35 配套 — 第三族擴展（B）＋ 重夸克標度檢定：K 之異常係唔係「輕區（Goldstone）」現象？
================================================================================
前線四路已封（docs/30–34）。本輪換一個問法：**失敗係局部於輕區，還是全面？**
做法 = 把模型（D 校準之 `(a,t)`、零新參數）擴展到**第三族 B**，並用**不依賴等距修正 `X`** 之接觸通道：

    c ≡ (Δ_{J=1} − Δ_{J=0})/4,   Δ_J ≡ M_J(中性) − M_J(帶電)

（`X` 在相減中消掉 ⟹ 無等距污染）。由 Fermi 接觸：`|c| m_h m_l/|q_h| ∝ |ψ(0)|² ≡ ψ²`。
三族 ψ² 之比 vs 模型之 `P₀` 之比 ⟹ 判「局部異常」定「全面失效」。

判準（事前寫死）：
  R1 `c` 一律用 X-free 定義；不得用 `X = 2.56`（避免等距污染）。
  R2 三族共用同一 `(a, t)`（D 校準 `A_D = 1.2025`、`|c|_D = 0.3525`），零新參數、零逐族調整。
  R3 偏差定義：`dev ≡ ψ²_data / ψ²_model`；`|ln dev| ≤ ln2`（2× 內）記「一致」。
  R4 「輕區異常」成立條件（預註冊）：K 之 `dev` > 5× 而 D、B 皆 ≤ 2×。
  R5 組分質量寫死：`m_l = 0.30`、`m_s = 0.495`、`m_c = 1.27`、`m_b = 4.18` GeV；`q_h = 1/3 (s)`、`2/3 (c, b)`。
  R6 B 用 `B⁰/B⁺`（PDG 511/521）與 `B*⁰/B*⁺`（513/523）。
輸出：`data/37_third_family.json`
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal

HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
B, L = 4, 40
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
A0, T0 = 0.384184, 1.7711e-3
A_D, C_D = 1.2025, 0.3525
ML = 0.30
BOUND_RHO, DPEAK_MAX = 1e-3, 20

print("=" * 104)
print("M5 docs/35 — 第三族（B）擴展 ＋ 重夸克標度檢定（X-free 接觸通道）")
print("=" * 104)

MASS, ERR = {}, {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    try:
        ev = abs(float(line[52:60]))
    except ValueError:
        ev = 0.0
    for k in ids:
        if k:
            MASS[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0

FAM = [("K", (311, 321), (313, 323), 1 / 3, 0.495),
       ("D", (411, 421), (413, 423), 2 / 3, 1.27),
       ("B", (511, 521), (513, 523), 2 / 3, 4.18)]

print("\n[R1] X-free 接觸幅度（Δ_J = M_中性 − M_帶電；X 相消）")
print(f"  {'族':3s}{'Δ_J0 (MeV)':>12s}{'Δ_J1 (MeV)':>12s}{'c = (Δ1−Δ0)/4':>16s}{'σ':>9s}{'m_h m_l':>10s}{'ψ² ∝ |c|m_hm_l/|q_h|':>22s}")
rows = {}
for nm, g0, g1, qh, mh in FAM:
    d0 = MASS[g0[0]] - MASS[g0[1]]
    d1 = MASS[g1[0]] - MASS[g1[1]]
    s0 = math.hypot(ERR[g0[0]], ERR[g0[1]])
    s1 = math.hypot(ERR[g1[0]], ERR[g1[1]])
    c = (d1 - d0) / 4
    sc = math.hypot(s0, s1) / 4
    psi = abs(c) * mh * ML / abs(qh)
    rows[nm] = {"d0": d0, "d1": d1, "c": c, "sigma_c": sc, "mh": mh, "psi2_proxy": psi}
    print(f"  {nm:3s}{d0:12.4f}{d1:12.4f}{c:16.5f}{sc:9.5f}{mh*ML:10.4f}{psi:22.6f}")

ref = rows["K"]["psi2_proxy"]
print("\n  ψ² 相對 K：  " + "   ".join(f"{k} = {v['psi2_proxy']/ref:.2f}" for k, v in rows.items()))
mu = {k: (v["mh"] * ML) / (v["mh"] + ML) for k, v in rows.items()}
print("  μ (GeV)：     " + "   ".join(f"{k} = {mu[k]:.5f}" for k in rows))
print("  μ³ 相對 K：   " + "   ".join(f"{k} = {(mu[k]/mu['K'])**3:.2f}" for k in rows))


def P0(a, t, q, mh, ml):
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, p0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    A = ALPHA * q * (HBARC / a) * rA * 1000
    c = Cq * (HBARC / a) ** 3 * p0 * 1000
    return A, c, p0, rA, dpk, bool(p0 > BOUND_RHO and dpk <= DPEAK_MAX)


print("\n[R2] 模型（同一 D 校準 a₀ = %.6f fm、t = %.4f MeV）逐族 P₀ 與 ψ²_data/ψ²_model" % (A0, T0 * 1000))
print(f"  {'族':3s}{'A_model':>10s}{'c_model':>10s}{'P₀_model':>10s}{'d_peak':>7s}{'ψ²_data':>11s}{'ψ²_model':>11s}{'dev':>9s}")
dev = {}
for nm, g0, g1, qh, mh in FAM:
    A, c, p0, rA, dpk, bd = P0(A0, T0, qh, mh, ML)
    psi_model = abs(c) * mh * ML / abs(qh)     # 同一 proxy，模型側
    d = rows[nm]["psi2_proxy"] / psi_model
    dev[nm] = d
    print(f"  {nm:3s}{A:10.4f}{c:10.5f}{p0:10.5f}{dpk:7d}{rows[nm]['psi2_proxy']:11.6f}{psi_model:11.6f}{d:9.3f}")

print("\n[R3/R4] 判：|ln dev| ≤ ln2 記一致；輕區異常 = K > 5× 且 D、B ≤ 2×")
for k, d in dev.items():
    tag = "✅ 一致" if abs(math.log(d)) <= math.log(2) else f"🔴 偏離 {max(d,1/d):.1f}×"
    print(f"  {k}：dev = {d:.3f} ⟹ {tag}")
# 修正（2026-10-27，append-only 註）：R4 原寫「dev > 5」，但 dev ≡ ψ²_data/ψ²_model 可小於 1；
# 偏離倍數須用 max(dev, 1/dev)。此為判準之**方向 bug**，非放寬標準（門檻 5×／2× 不變）。
fac = {k: max(v, 1 / v) for k, v in dev.items()}
light_anom = fac["K"] > 5 and fac["D"] <= 2 and fac["B"] <= 2
print(f"\n  ⟹ 輕區異常（R4）：{'✅ 成立' if light_anom else '🔴 不成立（D／B 亦偏離 ⟹ 全面失效）'}")

print("\n[R5] 另一角度：模型 P₀ 之排序 vs 數據 ψ² 之排序")
order_data = sorted(rows, key=lambda k: rows[k]["psi2_proxy"])
print(f"  數據 ψ² 排序（細→大）：{order_data}")
p0s = {}
for nm, g0, g1, qh, mh in FAM:
    p0s[nm] = P0(A0, T0, qh, mh, ML)[2]
order_model = sorted(p0s, key=lambda k: p0s[k])
print(f"  模型 P₀ 排序（細→大）：{order_model}")
print(f"  排序一致：{'✅' if order_data == order_model else '🔴'}")

out = {"rows": rows, "dev": dev, "orders": {"data": order_data, "model": order_model},
       "mu": mu, "light_anomaly": bool(light_anom),
       "calibration": {"a0_fm": A0, "t_MeV": T0 * 1000, "A_D": A_D, "C_D": C_D}}
with open("data/37_third_family.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/37_third_family.json")