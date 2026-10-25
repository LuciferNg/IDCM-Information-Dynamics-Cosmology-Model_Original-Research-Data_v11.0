"""M4 docs/26 配套 — em 承載窗定理（零參數）：8 通道電荷代數上界 vs 實測幅度
=========================================================================================
問題：M4 之 Q2／Q3／Q8 曾以「群論 no-go ＋ 四類形式擬合排除」判死；M5（`docs/23`）再提供
      「族分辨力定理」（同夸克內容必同 r）。本輪唔加任何新族、唔用任何擬合參數，只做一件事：

      對每個通道，由**電荷代數**給出 em 能量差之**嚴格上界**：
          |ΔE_em| ≤ (α ħc / a) · Q^max ,
          Q^max ≡ |Δ(Σ_i q_i²)| + |Δ(Σ_{i<j} q_i q_j)|
      （因 self 重疊與 pair 重疊皆 ≤ 1，且各項同號相加取絕對值上界；見 §定理陳述）
      ⟹ 該通道若真由 em 承載，必要條件為
          a ≤ a_i^max ≡ α ħc Q_i^max / A_i        （A_i = |E_i| 實測）
      ⟹ 全部 8 通道一致之**充要必要條件** = a ≤ min_i a_i^max ≡ **em 承載窗**。

驗證紀律：
  W1 電荷代數用 Fraction 精確（禁浮點比大小）。
  W2 通道夸克內容按 PDG 價夸克寫死（π⁺−π⁰ 屬湮滅型，另列並標記「代數上界不適用」）。
  W3 A_i 用 M4 `docs/18` 之 E 值（`E = Δm − coeff·X`），σ 一併傳播至 a_i^max。
  W4 唔准為放寬窗口而改通道集或夸克內容。
輸出：`data/26_em_window.json`
"""
import json
from fractions import Fraction as F

HBARC_MEVFM = 197.3269804
ALPHA = 1 / 137.035999084
ALPHA_HC = ALPHA * HBARC_MEVFM                       # MeV·fm

QU = F(2, 3)
QD = F(-1, 3)
QU_BAR = -QU
QD_BAR = -QD
QS = F(-1, 3)
QS_BAR = F(1, 3)
QC = F(2, 3)

# 通道：名稱、初態價夸克、終態價夸克（u↔d 方向：初 → 終）、E ± σ（docs/18）
CH = [
    ("n−p",              "udd",        "uud",         -1.2667, 0.0940),
    ("K⁺−K⁰",            "u sbar",     "d sbar",      -1.3740, 0.0963),
    ("Σ⁻−Σ⁰",            "dds",        "uds",         +2.2470, 0.1010),
    ("Σ⁰−Σ⁺",            "uds",        "uus",         +0.7120, 0.1196),
    ("Σ⁻−Σ⁺",            "dds",        "uus",         +2.9590, 0.2028),
    ("Ξ⁻−Ξ⁰",            "dss",        "uss",         +4.2900, 0.2301),
    ("D⁺−D⁰",            "c dbar",     "c ubar",      +2.2600, 0.1172),
    ("π⁺−π⁰（湮滅型）",  "u dbar",     "(u ubar − d dbar)/√2", +2.0336, 0.0940),
]
QMAP = {"u": QU, "d": QD, "s": QS, "c": QC, "ubar": QU_BAR, "dbar": QD_BAR,
        "sbar": QS_BAR, "cbar": -QC}


def parse(spec):
    """價夸克串 → 電荷列表；支援 'udd'／'u sbar'／'c dbar' 兩種寫法"""
    if "(" in spec:
        return None
    toks = []
    for tok in spec.split():
        if tok in QMAP:
            toks.append(tok)
        else:                                  # 逐字母（u/d/s/c）
            toks.extend(list(tok))
    return [QMAP[t] for t in toks]


def charges(qs):
    """回傳 (Σq², Σ_{i<j} q_i q_j)，精確分數"""
    s2 = sum(q * q for q in qs)
    pair = sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs)))
    return s2, pair


print("=" * 100)
print("em 承載窗定理（零參數）：8 通道之電荷代數上界 vs 實測幅度")
print(f"α ħc = {ALPHA_HC:.9f} MeV·fm")
print("=" * 100)
print(f"{'通道':<18}{'ΔΣq²':>8}{'ΔΣq_iq_j':>11}{'Q^max':>8}{'|E| ± σ':>14}{'a_i^max [fm]':>14}{'±':>9}")
rows = []
for name, ini, fin, E, sig in CH:
    qi, qf = parse(ini), parse(fin)
    if qi is None or qf is None:
        print(f"{name:<18}{'—':>8}{'—':>11}{'—':>8}{abs(E):>9.4f} ± {sig:<5}{'不適用':>12}")
        rows.append({"name": name, "type": "annihilation", "E": E, "sigma": sig,
                     "note": "混合電荷態（非單一 u↔d 置換）⟹ 代數上界不適用"})
        continue
    s2i, pi = charges(qi)
    s2f, pf = charges(qf)
    ds2, dpr = s2f - s2i, pf - pi
    qmax = abs(ds2) + abs(dpr)
    amax = float(F(1) / 1) * ALPHA_HC * float(qmax) / abs(E)
    damax = amax * sig / abs(E)
    print(f"{name:<18}{str(ds2):>8}{str(dpr):>11}{str(qmax):>8}{abs(E):>9.4f} ± {sig:<5.4}{amax:>12.4f}{damax:>9.4f}")
    rows.append({"name": name, "type": "flip", "dSigmaQ2": str(ds2), "dSigmaPair": str(dpr),
                 "Qmax": str(qmax), "E": E, "sigma": sig, "a_max": amax, "a_max_sigma": damax})

win = min(r["a_max"] for r in rows if r["type"] == "flip")
win_lo = min(r["a_max"] - r["a_max_sigma"] for r in rows if r["type"] == "flip")
win_hi = min(r["a_max"] + r["a_max_sigma"] for r in rows if r["type"] == "flip")
print("-" * 100)
print(f"⟹ em 承載窗：a ≤ {win:.4f} fm（σ 帶 [{win_lo:.4f}, {win_hi:.4f}]）；"
      f"限制通道 = {min((r for r in rows if r['type']=='flip'), key=lambda r: r['a_max'])['name']}")

print("\n對照 M5 兩個校準點（M5 docs/13、docs/19/23；已於 M4 docs/25 獨立驗證）")
for tag, a in (("κ=0 主校準（docs/13）", 0.377422), ("κ=0.215 天花板點（docs/23）", 0.243274)):
    print(f"  ── {tag}: a = {a:.6f} fm ⟹ {'🔴 窗外（em 結構上不可能承載全部通道）' if a > win else '🟡 窗內'}")
    for r in rows:
        if r["type"] != "flip":
            continue
        sat = abs(r["E"]) * a / (ALPHA_HC * float(F(r["Qmax"])))
        print(f"     {r['name']:<10} 所需飽和度 O = A·a/(αħc Q^max) = {sat:6.3f}"
              f" {'🔴 >1 不可達' if sat > 1 else '✅ ≤1'}")

print("\n單一 a 之下各通道所需 ⟨B⟩ 飽和度之離散度（窗內取 a = 0.15 fm）")
for a in (0.10, 0.15, 0.20):
    sats = [abs(r["E"]) * a / (ALPHA_HC * float(F(r["Qmax"]))) for r in rows if r["type"] == "flip"]
    print(f"  a = {a:.2f} fm：O ∈ [{min(sats):.3f}, {max(sats):.3f}] ⟹ 跨 {max(sats)/min(sats):.2f}×"
          f"（同一 (a,t) 兩體解須同時給出此跨度；M5 A3 已證 0/5 標度族命中）")

json.dump({"alpha_hc": ALPHA_HC, "rows": rows,
           "window": {"a_max_fm": win, "band": [win_lo, win_hi],
                      "limiting": min((r for r in rows if r["type"] == "flip"),
                                      key=lambda r: r["a_max"])["name"]}},
          open("data/26_em_window.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/26_em_window.json")