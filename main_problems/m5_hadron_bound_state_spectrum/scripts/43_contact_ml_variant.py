"""M5 docs/41 配套 — 接觸項 m_l 依賴：v3（heavy-quark 極限 p = 0）全鏈重校準
================================================================================
觸發：`docs/38` 判 ΔI=0 超精細失敗（`c(D_s)/c(D)` 170σ、`c(B_s)/c(B)` 16σ），反解需 `c ∝ m_l^{+0.30}`，
而 v2 用 `C_q ∝ 1/(m_h m_l)`（`p = −1`）。本輪定 `p = 0`（heavy-quark symmetry：輕夸克相對論化，
其質量不入接觸係數），命名 **v3**，並**全鏈重校準**（不得局部補丁）。

定義：
  v2：`C_q = (8π/3)·αq/(16·m_h·m_l)`   —— 文獻 Fermi contact（非相對論）
  v3：`C_q = (8π/3)·αq/(16·m_h·m_ref)` —— `m_ref = 0.30 GeV`（D 校準之輕夸克參考值，固定）
  兩者均用量綱完整（`C_q·(ħc/a)³` 為能量）⟹ `m_ref` 係**參考標度**，非新自由度（由 v2 之 `m_l,D` 取值固定）。

判準（事前寫死）：
  V1 校準目標一律 D 之兩數（`A_D = 1.2025`、`|c_D| = 0.3525` MeV），與 `docs/30` 同。
  V2 預言：ΔI=1 通道 `B`（`A_B = −2.4825`、`c_B = −0.0775`）、`K`（`A_K = 1.3335`、`c_K = −0.0135`）；
     ΔI=0 通道 `c(D_s)/c(D)`、`c(B_s)/c(B)`（PDG 直讀，同電荷 `E_hfs`）。
  V3 判優：v3 勝 ⟺ ①ΔI=0 兩比值皆 `|dev| ≤ 12%`；②ΔI=1 之 `B` dev 不劣於 v2（不得以 ΔI=1 換 ΔI=0）。
  V4 報所有 z（PDG 誤差），但**不**以 z 判優（K 側 σ 牆 `docs/37`）。
  V5 不得改 `q`、`m_h`、`X`、`α`、`b`、`ρ`；只改 `C_q` 之 `m_l` 依賴。
輸出：`data/43_contact_ml_variant.json`
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
A_D, C_D = 1.2025, 0.3525
M_REF = 0.30
BOUND_RHO, DPEAK_MAX = 1e-3, 20

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


def amps(a, t, q, mh, ml, variant):
    g = ALPHA * q * HBARC / a / t
    denom = 16 * mh * (ml if variant == "v2" else M_REF)
    Cq = (8 * np.pi / 3) * ALPHA * q / denom
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    _, evec = eigh_tridiagonal(diag, OFF, select="i", select_range=(0, 0))
    v = evec[:, 0] if evec.ndim == 2 else evec
    w = v ** 2
    w /= w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    dpk = int(np.argmax(w))
    return (ALPHA * q * (HBARC / a) * rA * 1000,
            Cq * (HBARC / a) ** 3 * P0 * 1000, P0,
            bool(P0 > BOUND_RHO and dpk <= DPEAK_MAX))


def bisect(f, lo, hi, n=200):
    flo = f(lo)
    if f(hi) * flo > 0:
        return None
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-13:
            break
    return 0.5 * (lo + hi)


def calib(variant, ngrid=320):
    def t_solve(a):
        x = bisect(lambda xx: amps(a, float(np.exp(xx)), 2 / 3, 1.27, 0.30, variant)[0] - A_D,
                   np.log(1e-14), np.log(1e4))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), ngrid):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = amps(a, t, 2 / 3, 1.27, 0.30, variant)[1] - C_D
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: amps(float(np.exp(xx)), t_solve(float(np.exp(xx))),
                                       2 / 3, 1.27, 0.30, variant)[1] - C_D, prev[0], x)
            if r is not None:
                aa = float(np.exp(r))
                return aa, t_solve(aa)
        prev = (x, k)
    return None


print("=" * 100)
print("M5 docs/41 — 接觸項 m_l 依賴：v3（p = 0）全鏈重校準")
print("=" * 100)

data = {
    "B": {"A": -2.4825, "c": -0.0775, "sA": 0.23354, "sc": 0.07554},
    "K": {"A": 1.33350, "c": -0.0135, "sA": 0.26341, "sc": 0.08216},
}
hfs = {"D": MASS[423] - MASS[421], "D_s": MASS[433] - MASS[431],
       "B": MASS[513] - MASS[511], "B_s": MASS[533] - MASS[531]}
print("\n[V1] D 校準")
res = {}
for variant in ("v2", "v3"):
    sol = calib(variant)
    a0, t0 = sol
    AD, cD, pD, bd = amps(a0, t0, 2 / 3, 1.27, 0.30, variant)
    print(f"  {variant}: a₀ = {a0:.6f} fm、t = {t0*1000:.6f} MeV、A_D = {AD:.6f}、"
          f"|c_D| = {cD:.6f}、P₀_D = {pD:.6f}（束縛 {bd}）")
    res[variant] = {"a0_fm": a0, "t_MeV": t0 * 1000, "P0_D": pD}
    print(f"     ΔI=1 預言")
    # docs/37 之約定：q = |q_h| = 2/3 (c, b)、1/3 (s)；c 報幅值（符號由定理另給）
    for fam, q, mh, ml in (("B", 2 / 3, 4.18, 0.30), ("K", 1 / 3, 0.495, 0.30)):
        A, c, P0, bd2 = amps(a0, t0, q, mh, ml, variant)
        c = abs(c)
        d = data[fam]
        devc = c / abs(d["c"]) - 1
        zc = (c - abs(d["c"])) / d["sc"]
        print(f"       {fam}: c = {c:8.5f}（dev {devc*100:+7.2f}%、z = {zc:+.2f}σ）"
              f"  [A 不比分：含 −X = −2.56 offset]")
        res[variant][f"devc_{fam}"] = devc
        res[variant][f"z_{fam}"] = zc
    print(f"     ΔI=0 預言（比值）")
    for pair, hi, lo, (qh, mhh, mlh), (ql, mhl, mll) in (
            ("D_s/D", "D_s", "D", (2 / 3, 1.27, 0.495), (2 / 3, 1.27, 0.30)),
            ("B_s/B", "B_s", "B", (2 / 3, 4.18, 0.495), (2 / 3, 4.18, 0.30))):
        star_hi, star_lo = (433, 423) if hi == "D_s" else (533, 513)
        ch = amps(a0, t0, qh, mhh, mlh, variant)[1]
        cl = amps(a0, t0, ql, mhl, mll, variant)[1]
        rd = hfs[hi] / hfs[lo]
        dev = (ch / cl) / rd - 1
        print(f"       {pair}: 模型 {ch/cl:.4f} vs 數據 {rd:.4f} ⟹ dev {dev*100:+6.2f}%  [z 見 scripts/40]")
        res[variant][f"dev_{pair}"] = dev
        res[variant][f"R_{pair}"] = ch / cl

v2, v3 = res["v2"], res["v3"]
improve = (abs(v3["dev_D_s/D"]) <= 0.12 and abs(v3["dev_B_s/B"]) <= 0.12)
no_loss = abs(v3["devc_B"]) <= abs(v2["devc_B"]) + 1e-9
print(f"\n[V3] v3 勝出條件：ΔI=0 兩比值 ≤12% → {'✅' if improve else '🔴'}；"
      f"ΔI=1 之 B 不劣化 → {'✅' if no_loss else '🔴'}")
verdict = "v3_adopt" if (improve and no_loss) else "v2_keep"
print(f"[總結] {verdict}")
out = {"results": res, "hfs": hfs, "m_ref": M_REF, "verdict": verdict,
       "criteria": {"improve": bool(improve), "no_loss": bool(no_loss)}}
with open("data/43_contact_ml_variant.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("-> data/43_contact_ml_variant.json")
