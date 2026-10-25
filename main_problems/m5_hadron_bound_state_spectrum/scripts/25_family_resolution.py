"""M5 docs/23 配套 — D5 終局：族間分辨力定理（r = A/|c| 只依賴夸克質量 ⟹ 同族必須同 r）
================================================================================================
`docs/22` 把殘餘 1/6.14 歸屬到 A3 之層指派。本輪換一個**與層指派無關**嘅刀刃：

`docs/14` §四.1 恒等式：`A/|c|` 與電荷無關，只由 `(a, t, m_h, m_l)` 決定。
而 `(a, t)` 係全域標度 ⟹ **同一夸克內容（同 `m_h`、同 `m_l`）之所有族，模型必給同一個 `r ≡ A/|c|`**。

數據（`data/08_sampleout_channels.json`，6 族；`A = (E1+3E2)/4`、`c = (E1−E2)/4`）：

| 族 | 夸克內容 | `r` 數據 |
|:--|:--|--:|
| K vs K* | (s, ud) | **+98.778** |
| Σ vs Σ(1385) [−/0] | (s, ud) | **+3.877** |
| Σ vs Σ(1385) [0/+] | (s, ud) | **−1.814** |
| D vs D* | (c, ud) | **+3.411** |
| Σc(2455) vs Σc(2520) | (c, ud) | **+47.065** |
| Ξc vs Ξc′ | (c, s) | **−2.695** |

判準（寫死）：
  Z1 族分類只按夸克內容（`m_h`, `m_l`），不得引入自旋／徑向激發之額外構造。
  Z2 模型 `r` 由既有校準（D 兩數）定 ⟹ 同族預言**唯一**，不得逐族調參。
  Z3 命中 = 同族內所有族之 `r` 與模型值差 ≤ 2σ（σ 由 E1／E2 之 σ 傳播）。
  Z4 若同族內 `r` 數據本身跨 > 3× ⟹ 模型（`r` 為質量函數）**結構上不可能**同時命中 ⟹ 記定理。
  Z5 禁：不得為救任何一族而改 `m_h`／`m_l`（`m_l` 全域 0.30；`m_s = 0.495`、`m_c = 1.27` 為既有輸入）。

輸出：`data/25_family_resolution.json`
"""
import json
from fractions import Fraction as F

import numpy as np

B, L = 4, 40
HBARC, ALPHA = 0.1973269804, 1 / 137.035999084
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = float(F(2581, 352)), float(F(1447, 352))

FAMILIES = [
    ("K vs K*", "s_ud", 1.374, 1.320, 0.0960728889958036, 0.34122719704032967),
    ("Σ vs Σ(1385) [−/0]", "s_ud", 2.247, 0.940, 0.10125709851659784, 1.121978609421766),
    ("Σ vs Σ(1385) [0/+]", "s_ud", 0.712, -1.690, 0.11963277142990543, 1.0603942662990968),
    ("D vs D*", "c_ud", 2.260, 0.850, 0.1176265276202609, 0.1176265276202609),
    ("Σc(2455) vs Σc(2520)", "c_ud", -3.880, -3.570, 0.2771930735065362, 0.7397540131692426),
    ("Ξc vs Ξc′", "c_s", 0.170, -2.060, 0.3743474322070341, 0.7133274143056608),
]
CLASS = {"s_ud": (M_S, M_L, Q_S, G0), "c_ud": (M_C, M_L, Q_C, G7), "c_s": (M_C, M_S, Q_C, G7)}


def shape(a, t, q, mh, ml):
    """v2（2026-10-25 審計修正）：加回 scripts/15／19 之束縛守門（bound = ⟨B⟩>1e-3 ∧ d_peak≤20）
    —— v1 缺此守門，s_ud／c_s 之解為尾端帶邊態（d_peak 36–38）而被當作束縛態使用 ⟹ 舊 r 值作廢。"""
    g = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    return (rA, P0,
            ALPHA * q * (HBARC / a) * rA * 1000,
            Cq * (HBARC / a) ** 3 * P0 * 1000,
            int(np.argmax(w)), bool(rA > 1e-3 and int(np.argmax(w)) <= 20))


def bisect(f, lo, hi, n=300):
    flo = f(lo)
    if flo * f(hi) > 0:
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


def calib(kappa):
    chd = kappa * G7 * Q_C
    target = 1.2025 - chd
    if target <= 0:
        return None

    def t_solve(a):
        x = bisect(lambda xx: shape(a, float(np.exp(xx)), Q_C, M_C, M_L)[2] - target, np.log(1e-12), np.log(1e3))
        return float(np.exp(x)) if x is not None else None

    prev = None
    for x in np.linspace(np.log(1e-3), np.log(3.0), 400):
        a = float(np.exp(x))
        t = t_solve(a)
        if t is None:
            prev = None
            continue
        k = shape(a, t, Q_C, M_C, M_L)[3] - 0.3525
        if prev is not None and prev[1] * k <= 0:
            r = bisect(lambda xx: shape(float(np.exp(xx)), t_solve(float(np.exp(xx))), Q_C, M_C, M_L)[3] - 0.3525, prev[0], x)
            if r is not None:
                return float(np.exp(r)), t_solve(float(np.exp(r)))
        prev = (x, k)
    return None


print("=" * 104)
print("M5 docs/23 — D5：族間分辨力定理（r = A/|c| 為質量函數 ⟹ 同族必須同 r）")
print("=" * 104)
print(f"\n[Z1] 數據逐族 r（精確由 E1／E2 反解）")
print(f"     {'族':26s}{'內容':>7s}{'A':>10s}{'c':>10s}{'r':>11s}{'σ(r)':>9s}")
rows = []
for name, cls, E1, E2, s1, s2 in FAMILIES:
    A, c = (E1 + 3 * E2) / 4, (E1 - E2) / 4
    sr = abs(A / c) * np.sqrt((s1 / (E1 - E2)) ** 2 + (s2 / (E1 - E2)) ** 2) if c != 0 else np.inf
    rows.append({"family": name, "class": cls, "A": A, "c": c, "r": A / c, "sigma_r": sr})
    print(f"     {name:26s}{cls:>7s}{A:10.4f}{c:10.4f}{A/c:11.3f}{sr:9.3f}")

print(f"\n[Z4] 同族內 r 之跨度（模型對同族只可給一個值）")
for cls in ("s_ud", "c_ud"):
    sub = [r for r in rows if r["class"] == cls]
    mags = [abs(r["r"]) for r in sub]
    if len(sub) > 1:
        print(f"     {cls}：r = {[round(r['r'],3) for r in sub]} ⟹ |r| 跨度 = {max(mags)/min(mags):.2f}×"
              f" ⟹ {'🔴 > 3× ⟹ 模型結構上不可能（Z4）' if max(mags)/min(mags) > 3 else '✅'}")

print(f"\n[Z2] 模型逐類 r（κ = 0.2150 校準點；(a, t) 由 D 兩數定 ⟹ 各類唯一）")
sols = {}
for kp in (0.2150, 0.26986):
    sol = calib(kp)
    if sol is None:
        print(f"     κ = {kp:.5f}：D 校準無解 ⟹ N/A")
        continue
    a, t = sol
    line = {}
    print(f"     κ = {kp:.5f}（a₀ = {a:.6f} fm、t = {t*1000:.4f} MeV）：")
    for cls, (mh, ml, q, gch) in CLASS.items():
        rA, P0, A, c, dpk, bd = shape(a, t, q, mh, ml)
        ch = kp * gch * q
        r = (A + ch) / c if bd else float("nan")
        line[cls] = r
        tag = "✅ 束縛" if bd else "🔴 非束縛（帶邊，d_peak=%d）⟹ r 未定義" % dpk
        rr = f"{r:9.3f}" if bd else "      N/A"
        print(f"        {cls:6s}（m_h = {mh:.3f}、m_l = {ml:.3f}）⟹ ⟨B⟩ = {rA:.4f}、P₀ = {P0:.4f}、"
              f"d_peak = {dpk:2d} {tag}、A = {A+ch:8.4f}、c = {c:8.4f} ⟹ **r = {rr}**")
    sols[str(kp)] = line

print(f"\n[Z3] 對照（用 |r| 與其 σ 傳播；禁以點估計之比作判）")
zs = {}
for kp, line in sols.items():
    print(f"     κ = {kp}：")
    chi = {}
    for r in rows:
        if line[r["class"]] != line[r["class"]]:      # nan ⟹ 非束縛，模型無預言
            print(f"        {r['family']:26s} 模型 r = N/A（非束縛）⟹ 不入 χ²")
            continue
        pred = abs(line[r["class"]])
        obs, sd = abs(r["r"]), r["sigma_r"]
        z = (pred - obs) / sd if sd not in (0, float("inf")) else 0.0
        chi.setdefault(r["class"], []).append((r["family"], z))
        print(f"        {r['family']:26s} 模型 |r| = {pred:9.3f}、數據 |r| = {obs:8.3f} ± {sd:8.3f}"
              f" ⟹ z = {z:+8.2f}")
    for cls, lst in chi.items():
        c2 = sum(z * z for _, z in lst)
        print(f"        {cls:6s} 類 χ² = {c2:10.2f}（{len(lst)} 族，0 自由參數）⟹ "
              f"{'🔴 排除（%.1fσ）' % np.sqrt(c2) if np.sqrt(c2) > 5 else '🟡 未排除'}")
    zs[kp] = {cls: sum(z * z for _, z in lst) for cls, lst in chi.items()}

out = {"data_r": rows, "model_r": sols, "class_chi2": zs,
       "verdict": "🔴 v2 審計修正：加束縛守門後 s_ud／c_s 於 D 校準點（兩 κ 皆）為非束縛帶邊態（d_peak 35–38）⟹ 模型對該兩類無預言；v1 之 13.403／73.963 及 χ² 14.7σ／73.2σ 作廢。結構定理（r 只依 m_h·m_l）不受影響"}
with open("data/25_family_resolution.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/25_family_resolution.json")