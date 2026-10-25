"""M4 docs/28 配套 — 強 IB 矩陣元 ⟨O⟩ ≠ 1：docs/27 之 14.6σ 異常之定量終結
=========================================================================================
外部錨（**權威來源，已下載原文核對**）：
  Borsányi et al. (BMW), PRL 111, 252001 (2013) = arXiv:1306.2287v2, **TABLE I**（MeV）：
      X      ΔM_X             ΔQED M_X          ΔQCD M_X
      N    −0.68(39)(36)    +1.59(30)(35)     −2.28(25)(7)
      Σ    −7.84(87)(72)    +0.08(12)(34)     −7.67(79)(105)
      Ξ    −7.16(76)(47)    −1.29(15)(8)      −5.87(76)(43)
  其符號約定：`ΔM_X ≡ M(高電荷) − M(低電荷)`（N: p−n；Σ: Σ⁺−Σ⁻；Ξ: Ξ⁰−Ξ⁻）。
  其分解：`ΔM_X = A_X α + B_X ΔM²`（`ΔM² ≡ M²_ūu − M²_dd`，由 kaon `ΔM²_K` 校準）。
  ⟹ `ΔQCD M_X` = 夸克質量差部分；`ΔQED M_X` = 電磁部分。

本輪定義（M4 取向，與 docs/18 之 `Δm` 同向）：
  `Δm_M4(ch) = M(第一態) − M(第二態)`；`coeff` = u↔d 置換數
  `⟨O⟩_ch ≡ |ΔQCD| / (coeff · X)`，`X = 2.56 ± 0.094`（FLAG24 `m_d − m_u`）
  ⟹ `em_M4(ch) = Δm_M4 − coeff · X · ⟨O⟩_ch`（等效 `= −ΔQED`，取向見下）

驗證紀律：
  M1 外部值逐位抄自 TABLE I（含兩段誤差，quadrature 合成）。
  M2 取向逐通道明列，不得靠記憶：`ΔM_X = M(高) − M(低)` ⟹ `em_M4 = −ΔQED M_X`。
  M3 三個 `⟨O⟩` 中，**`⟨O⟩_K` 不可得**：BMW 以 `ΔM²_K` 作輸入校準 `δm` ⟹ K 通道係**入**而非**出**，不得當證據。
  M4 重新計算 docs/26 之承載窗（把 `|em|` 換成含 `⟨O⟩` 之值），並報告窗之變化。
輸出：`data/28_matrix_element.json`
"""
import json
from fractions import Fraction as F

HBARC, ALPHA = 197.3269804, 1 / 137.035999084
AH = ALPHA * HBARC
X, SX = 2.56, 0.094
QF, QS_, QN = 2 / 3, 1.0, 1.0            # docs/26 之 Q^max（n−p 2/3、Σ 系 1、Ξ 1）

# BMW TABLE I（stat, sys）；取向 = M(高電荷) − M(低電荷)
BMW = {
    "N":  dict(total=(-0.68, 0.39, 0.36), qed=(+1.59, 0.30, 0.35), qcd=(-2.28, 0.25, 0.07)),
    "Σ":  dict(total=(-7.84, 0.87, 0.72), qed=(+0.08, 0.12, 0.34), qcd=(-7.67, 0.79, 1.05)),
    "Ξ":  dict(total=(-7.16, 0.76, 0.47), qed=(-1.29, 0.15, 0.08), qcd=(-5.87, 0.76, 0.43)),
}
# 實驗總分裂（PDG），同取向：M(高電荷) − M(低電荷)
EXP = {"N": -1.29333, "Σ": -8.0790, "Ξ": -6.8500}

# M4 通道：對應 BMW 之 X、coeff（u↔d 置換數）、Δm_M4、Q^max、E's σ
CH = [
    ("n−p",     "N", 1, +1.29330, QF, 0.0940),
    ("Σ⁻−Σ⁺",   "Σ", 2, +8.07900, QS_, 0.2028),
    ("Σ⁻−Σ⁰",   "Σ", 1, +4.80700, QS_, 0.1010),
    ("Σ⁰−Σ⁺",   "Σ", 1, +3.27200, QF, 0.1196),
    ("Ξ⁻−Ξ⁰",   "Ξ", 1, +6.85000, QN, 0.2301),
]


def comb(t):
    v, a, b = t
    return v, (a * a + b * b) ** 0.5


print("=" * 104)
print("⟨O⟩ ≠ 1：以 BMW PRL111,252001 (arXiv:1306.2287v2) TABLE I 定量終結 docs/27 之 14.6σ")
print(f"α ħc = {AH:.9f} MeV·fm ；X = {X} ± {SX} MeV（FLAG24 m_d − m_u）")
print("=" * 104)
print(f"{'X':<4}{'ΔM_X(exp)':>12}{'ΔQED':>16}{'ΔQCD':>18}{'⟨O⟩ = |ΔQCD|/(coeff·X)':>26}{'BMW 自洽殘差':>14}")
O, rows = {}, []
for k, d in BMW.items():
    tq, sq = comb(d["qcd"])
    tq_, sq_ = comb(d["qed"])
    nflip = {"N": 1, "Σ": 2, "Ξ": 1}[k]
    o = abs(tq) / (nflip * X)
    so = o * ((sq / abs(tq)) ** 2 + (SX / X) ** 2) ** 0.5
    O[k] = (o, so)
    resid = EXP[k] - comb(d["total"])[0]
    print(f"{k:<4}{EXP[k]:>12.4f}{tq_:>9.2f} ± {sq_:<5.2f}{tq:>12.2f} ± {sq:<6.2f}"
          f"{o:>17.3f} ± {so:<6.3f}{resid:>+9.3f}")
    rows.append({"sector": k, "exp": EXP[k], "qed": [tq_, sq_], "qcd": [tq, sq],
                 "O": o, "O_sigma": so, "bmw_vs_exp_residual": resid, "n_flip": nflip})

print("\n【判讀 1】⟨O⟩ 是否為 1？（若 ⟨O⟩ ≡ 1 則 docs/27 之算法成立）")
for k, (o, so) in O.items():
    print(f"  {k}：⟨O⟩ = {o:.3f} ± {so:.3f} ⟹ 距 1 有 {abs(o-1)/so:.1f}σ ⟹ "
          f"{'🔴 顯著 ≠ 1' if abs(o-1)/so > 3 else ('🟡 偏離' if abs(o-1)/so > 1 else '✅ 與 1 相容')}")

print("\n【判讀 2】以 ⟨O⟩ 重算 em（= Δm_M4 − coeff·X·⟨O⟩），並與 docs/27 之符號定理預言對照")
SIGN_PRED = {"n−p": -1, "Σ⁰−Σ⁺": -1}          # docs/27 有判別力之通道（同號）；其餘無預言
print(f"{'通道':<9}{'coeff':>6}{'Δm_M4':>9}{'⟨O⟩':>8}{'em(M4)':>10}{'±':>7}"
      f"{'lattice −ΔQED':>14}{'em 符':>7}{'Coulomb 預言':>13}{'判':>7}")
for name, sec, nf, dm, q, sE in CH:
    o, so = O[sec]
    em = dm - nf * X * o
    sem = ((nf * X * so) ** 2 + (nf * o * SX) ** 2) ** 0.5
    qed = -BMW[sec]["qed"][0]                  # M4 取向
    pred = SIGN_PRED.get(name, 0)
    if pred:
        verdict = "PASS" if (em > 0) == (pred > 0) else "FAIL"
    else:
        verdict = "無預言"
    print(f"{name:<9}{nf:>6}{dm:>9.4f}{o:>8.3f}{em:>+10.4f}{sem:>7.3f}"
          f"{qed:>+14.3f}{('+' if em>0 else '−'):>7}"
          f"{(('+' if pred>0 else '−') if pred else '—'):>13}{verdict:>7}")
    rows.append({"channel": name, "coeff": nf, "dm_M4": dm, "O": o, "em": em, "em_sigma": sem,
                 "qed_M4_orientation": qed, "coulomb_pred": pred, "verdict": verdict})

print("\n  ⟹ `em(M4)` 與 `−ΔQED` 之差 = BMW 自身與實驗之張力（N −0.61／Σ −0.24／Ξ +0.31），非 M4 之殘餘。")

print("\n【判讀 3】`⟨O⟩_K` 不可由本表得出（BMW 以 ΔM²_K 作 δm 之輸入 ⟹ K 係入非出）")
o_need = 3.9340 / X
print(f"  docs/27 之 K 符號要求：`X·⟨O⟩_K > 3.9340` ⟹ `⟨O⟩_K > {o_need:.3f}`（此為必要條件，非量度）")

print("\n【判讀 4】重算 docs/26 之承載窗：a_i^max = α ħc · Q^max / |em_i|")
print(f"{'通道':<9}{'|em|':>9}{'Q^max':>8}{'a_i^max [fm]':>14}")
wins = []
for name, sec, nf, dm, q, sE in CH:
    o, _ = O[sec]
    em = abs(dm - nf * X * o)
    if em < 1e-9:
        print(f"{name:<9}{em:>9.4f}{q:>8.2f}{'∞':>14}")
        continue
    amax = AH * q / em
    wins.append((name, amax))
    print(f"{name:<9}{em:>9.4f}{q:>8.2f}{amax:>14.4f}")
wmin = min(w for _, w in wins)
print(f"  ⟹ 新窗 a ≤ {wmin:.4f} fm（舊窗 docs/26 = 0.3357 fm）")
for tag, a in (("M5 主校準 0.377422", 0.377422), ("M5 天花板點 0.243274", 0.243274)):
    print(f"     對照 {tag} fm ⟹ {'🔴 窗外' if a > wmin else '✅ 窗內（舊判 no-go 失效）'}")

json.dump({"BMW_ref": "PRL 111, 252001 (2013) / arXiv:1306.2287v2 TABLE I",
           "X": X, "X_sigma": SX, "rows": rows,
           "O": {k: {"value": v[0], "sigma": v[1]} for k, v in O.items()},
           "window_new_fm": wmin,
           "O_K_required_min": o_need},
          open("data/28_matrix_element.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/28_matrix_element.json")