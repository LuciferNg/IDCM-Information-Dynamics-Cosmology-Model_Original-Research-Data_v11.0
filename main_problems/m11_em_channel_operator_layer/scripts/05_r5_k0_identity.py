#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M11 scripts/05 — R5：Q6 `k₀` 身份 ＋ R4 判詞之 k₀ 依賴（靈敏度定理）

問題（Q6）：`k₀ = α ħc/a₀` 之 `a₀` 能否由 IDCM 譜零參數定出？
      若不能，R4 之判詞（尤其 T-N4 / D 之唯一性）在多大 `k₀` 範圍內成立？

做法：
  A 溯源：M11 用之 `ħc/a₀ = 0.514 GeV` vs M5 校準 `a₀ = 0.377422 fm`（⟹ 0.52277 GeV）之差
  B 解析：對每個（模型, 盲測通道）解出 `k₀^crit`（使箱約束剛好可行之 k₀）
  C 可行域：代入 corpus 之各 `a` 支（窗界 0.3357σ帶、M5 κ=0、κ=0.215），逐支重判
  D 零參數性：清點 v11 存貨內所有長度量級，檢有無 fm 級結構候選

判準（事前寫死）：
  R5-1 `k₀^crit = |em|/|box 界|`（解析；箱界用 Fraction 精確）
  R5-2 一模型在給定 `k₀` 之盲測判 = 三通道全部可行 ⟹ PASS_BLIND；否則 OUT
  R5-3 錨層判詞（z）與 `k₀` 無關 ⟹ 標明「硬判」
  R5-4 零參數性：存貨內若無 fm 級（0.1–1 fm）結構量 ⟹ 判「中間標度缺口」，禁以量級接近作匹配
禁：不得以「量級接近」（如 4.36 fm vs 0.38 fm）宣告推導。
"""

import json
import os
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D1 = json.load(open(os.path.join(HERE, "data", "01_r1_criteria.json"), encoding="utf-8"))
D4 = json.load(open(os.path.join(HERE, "data", "04_r4_necessity.json"), encoding="utf-8"))
OUT = os.path.join(HERE, "data", "05_r5_k0_identity.json")

CH9 = D1["C9"]
ROWS, C9 = D1["C2"]["rows"], CH9["channels"]
X = D1["C2"]["X"]
ALPHA_INV = 137.035999
HBARC = 197.3269804          # MeV·fm
BLIND = ["S0-S+", "S--S+", "D+-D0"]
BOX = {c: (F(C9[c]["box_min"]).limit_denominator(12), F(C9[c]["box_max"]).limit_denominator(12))
       for c in BLIND}


def k0_from_a(a_fm):
    return (HBARC / a_fm) / ALPHA_INV


def feas(em, k0, ch, tol=1e-12):
    """t = em/k0 是否落箱（箱界以 Fraction 精確載入，比較用 float）"""
    lo, hi = float(BOX[ch][0]), float(BOX[ch][1])
    t = em / k0
    return (lo - tol) <= t <= (hi + tol)


def main():
    log = []
    def emit(s=""):
        log.append(s); print(s)

    emit("M11 R5 — Q6：k₀ 身份 ／ R4 判詞之 k₀ 依賴")
    emit("=" * 104)
    emit(f"α⁻¹ = {ALPHA_INV}｜ħc = {HBARC} MeV·fm｜X = {X}｜k₀ = αħc/a₀")

    # ── A 溯源 ───────────────────────────────────────────────
    emit("")
    emit("A｜溯源：三支 ħc/a（沿用值 vs 現行正典）")
    a_m5 = 0.37742216816391858835          # M5 docs/13 重校（現行正典）
    a_v2 = 0.384184                        # M5 docs/28（舊 v2，被重校取代）
    a_old = HBARC / 514.0                  # 舊版 M11 所引（v2 之四捨五入）
    hc_a_m5, hc_a_v2, hc_a_old = HBARC / a_m5, HBARC / a_v2, HBARC / a_old
    k0_m5, k0_v2, k0_old = k0_from_a(a_m5), k0_from_a(a_v2), k0_from_a(a_old)
    emit(f"  ① 現行正典（M5 docs/13 重校）a₀ = {a_m5:.8f} fm ⟹ ħc/a₀ = {hc_a_m5:.4f} MeV ⟹ k₀ = {k0_m5:.6f} MeV  ← data/01 現用")
    emit(f"  ② M5 v2（docs/28，已取代）  a  = {a_v2:.6f} fm ⟹ ħc/a  = {hc_a_v2:.4f} MeV ⟹ k₀ = {k0_v2:.6f} MeV")
    emit(f"  ③ 舊版 M11（重算前）        ħc/a = 514.0 MeV           ⟹ a   = {a_old:.6f} fm ⟹ k₀ = {k0_old:.6f} MeV")
    emit(f"  ⟹ ③ = ② 之四捨五入（非獨立校準）；①↔③ 差 Δ(a) = {(a_old/a_m5-1)*100:+.2f}%、Δ(k₀) = {(k0_old/k0_m5-1)*100:+.2f}%")
    a_m11 = CH9["a0_fm"]
    emit(f"  核對：data/01 之 a₀ = {a_m11:.8f} fm、k₀ = {CH9['k0_MeV']:.6f} MeV ⟹ "
         f"{'與①一致 ✓' if abs(a_m11-a_m5) < 1e-12 else '與①不符 ✗'}")

    # ── B 解析 k₀^crit ──────────────────────────────────────
    emit("")
    emit("B｜R5-1 `k₀^crit`（|em|/|箱界|）：模型在 k₀ ≥ k₀^crit 時盲測可行")
    crit = {}
    for m in D4["models"]:
        lab = m["label"].split()[0]
        rows = []
        for ch in BLIND:
            v = m["verdicts"][ch]
            if v["O"] is None:
                rows.append((ch, None, None)); continue
            em, k0_used = v["em"], CH9["k0_MeV"]
            t = em / k0_used
            lo, hi = BOX[ch]
            # em/k0 需落 [lo,hi] ⟹ em<0：k0 ≥ |em|/|lo| ；em>0：k0 ≥ em/hi
            if em < 0:
                kc = abs(em) / float(abs(lo))
            else:
                kc = em / float(hi) if hi != 0 else float("inf")
            rows.append((ch, kc, t))
        worst = max([r[1] for r in rows if r[1] is not None], default=None)
        crit[lab] = worst
        emit(f"  {m['label']:34s}")
        for ch, kc, t in rows:
            if kc is None:
                emit(f"      {ch:8s} 無預言（新類別）")
            else:
                mark = ""
                if kc == max([r[1] for r in rows if r[1] is not None]):
                    mark = "  ← 約束通道"
                emit(f"      {ch:8s} t(k₀={k0_used:.4f})={t:+.4f}  k₀^crit = {kc:.4f} MeV{mark}")
        emit(f"      ⟹ 全盲測可行需 k₀ ≥ {worst:.4f} MeV  (a ≤ {HBARC/(worst*ALPHA_INV):.5f} fm)")
    kmax = max(v for v in crit.values() if v is not None)
    emit("")
    emit(f"  ★ 全局門檻 k₀^crit(all) = {kmax:.4f} MeV ⟺ a = {HBARC/(kmax*ALPHA_INV):.5f} fm")
    emit(f"    ⟹ 當 k₀ ≥ {kmax:.4f}（a ≤ {HBARC/(kmax*ALPHA_INV):.5f} fm）時，四模型（D 以其兩條可判通道計）**全部**通過盲測箱"
         f" ⟹ 箱之鑑別力歸零")

    # ── C 各 a 支重判 ────────────────────────────────────────
    emit("")
    emit("C｜corpus 各 a 支之重判（錨層硬判不變；盲測層隨 k₀ 變）")
    branches = [
        ("M5 κ=0 主校準（正典）", 0.37742216816391858835),
        ("M5 v2 docs/28（舊）", 0.384184),
        ("窗界（M4 docs/26）", 0.3357),
        ("M5 κ=0.215 天花板", 0.243274),
        ("窗內下限（示意 0.20）", 0.20),
        ("窗內更小（示意 0.10）", 0.10),
    ]
    verdict = {}
    for name, a in branches:
        k0 = k0_from_a(a)
        line = f"  a = {a:.6f} fm ⟹ k₀ = {k0:7.4f} MeV ｜"
        for m in D4["models"]:
            lab = m["label"].split()[0]
            ok = all(True if m["verdicts"][ch]["O"] is None else
                     feas(m["verdicts"][ch]["em"], k0, ch) for ch in BLIND)
            # 綜合 = 錨層（k₀ 無關）＋ 盲測層
            ga = m["grade_anchor"]
            if ga == "KILL":
                st = "錨死"                      # 硬判，與 k₀ 無關
            elif not ok:
                st = "盲測死"
            elif ga == "TENSION":
                st = "錨爭議"                    # 未達 3σ，不算死
            else:
                st = "存活"
            line += f"{lab}:{st}  "
            verdict.setdefault(name, {})[lab] = dict(blind_ok=bool(ok), anchor=ga,
                                                     alive=bool(ok and ga in ("SURVIVE", "SATURATED")))
        emit(line)

    # ── D 零參數性清點 ──────────────────────────────────────
    emit("")
    emit("D｜R5-4 存貨長度清點（檢有無 fm 級零參數候選）")
    inv = [
        ("R_world（世界環半徑）", 4.4e3 * 3.0857e19, "Mpc→fm", "10²⁹"),
        ("r_d = πR_world/88", 147.4386 * 3.0857e19, "Mpc→fm", "10²⁷"),
        ("λ_sync = L_ring/88", 147.4386 * 3.0857e19, "Mpc→fm", "10²⁷"),
        ("ħc/Λ_IDCM（Λ=45.28487 MeV）", HBARC / 45.28487, "fm", "10⁰"),
        ("ħc/v_SM（v=246.2196 GeV）", HBARC / 246.2196e3, "fm", "10⁻⁴"),
        ("ℓ_Pl", 1.616e-35, "m→fm", "10⁻²⁰"),
    ]
    for name, val, unit, mag in inv:
        fm = val * 1e15 if "m→fm" in unit else val
        emit(f"  {name:34s} ≈ {fm:.3e} fm  ({mag})")
    emit(f"  ⟹ 唯一落 fm 數量級者 = ħc/Λ_IDCM = {HBARC/45.28487:.3f} fm，距所需 a ≈ {a_m11:.3f} fm"
         f" 差 {HBARC/45.28487/a_m11:.2f}×")
    emit("  ⟹ 存貨內**無**任何 0.1–1 fm 之結構長度（除 ħc/Λ，差 1 個數量級以上且比值非結構整數）")
    emit("  ⟹ `a₀` 屬『中間標度缺口』（與 O38 `M₀` 同族）；`k₀` 不可零參數化（禁以量級接近作匹配，R5-4）")

    # ── E X 靈敏度（R5-5）：X 之 1σ 如何移 k₀^crit ──────────
    emit("")
    emit("E｜R5-5 `X = 2.56 ± 0.094` 對 k₀^crit 之 1σ 影響（em 含 X，故 k₀^crit ∝ |em|）")
    xsens = {}
    for m in D4["models"]:
        lab = m["label"].split()[0]
        vals = {}
        for sx in (-1, 0, +1):
            dX = 0.094 * sx
            need = 0.0
            for ch in BLIND:
                v = m["verdicts"][ch]
                if v["O"] is None:
                    continue
                row = ROWS[ch]
                em_s = row["dm"] - row["coeff"] * (2.56 + dX) * v["O"]
                lo = float(BOX[ch][0])
                kc = abs(em_s) / abs(lo) if em_s < 0 else (em_s / float(BOX[ch][1]) if float(BOX[ch][1]) else float("inf"))
                need = max(need, kc)
            vals[f"{sx:+d}"] = need
        xsens[lab] = vals
        emit(f"  {m['label']:34s} X−1σ: {vals['-1']:7.4f}  X: {vals['+0']:7.4f}  X+1σ: {vals['+1']:7.4f} MeV"
             f"  ⟹ a^crit: {HBARC/(vals['-1']*ALPHA_INV):.5f} / {HBARC/(vals['+0']*ALPHA_INV):.5f} / {HBARC/(vals['+1']*ALPHA_INV):.5f} fm")
    gc = xsens["C"]
    emit(f"  ★ 對 C 之排除判詞：X−1σ ⟹ a > {HBARC/(gc['-1']*ALPHA_INV):.5f} fm 才排除；"
         f"X+1σ ⟹ a > {HBARC/(gc['+1']*ALPHA_INV):.5f} fm")
    emit("    逐支裁決（C 被排除 ⟺ a > a^crit(X)）：")
    for nm, a in (("κ=0 正典", 0.37742216816391858835), ("v2 docs/28", 0.384184), ("窗界", 0.3357), ("κ=0.215 支", 0.243274)):
        st = []
        for sx in ("-1", "+0", "+1"):
            acrit = HBARC / (gc[sx] * ALPHA_INV)
            st.append(f"X{sx}:{'排除' if a > acrit else '存活'}(界{acrit:.4f})")
        emit(f"      a={a:.6f} ({nm:11s}) ⟹ " + "  ".join(st))

    res = dict(alpha_inv=ALPHA_INV, hbarc=HBARC, x_sensitivity=xsens,
               provenance=dict(canonical=dict(a=a_m5, hc_over_a=hc_a_m5, k0=k0_m5, src="M5 docs/13 重校"),
                               v2=dict(a=a_v2, hc_over_a=hc_a_v2, k0=k0_v2, src="M5 docs/28（已取代）"),
                               m11_before=dict(a=a_old, hc_over_a=hc_a_old, k0=k0_old, src="舊版 M11：514.0 MeV"),
                               data01=dict(a=a_m11, k0=CH9["k0_MeV"]),
                               delta_a_pct=(a_old / a_m5 - 1) * 100,
                               delta_k0_pct=(k0_old / k0_m5 - 1) * 100),
               k0_crit=crit, k0_crit_all=kmax, a_crit_all=HBARC / (kmax * ALPHA_INV),
               branches={n: dict(a=a, k0=k0_from_a(a), verdict=v) for (n, a), v in
                         zip(branches, [verdict[n] for n, _ in branches])},
               inventory_scales=[dict(name=n, fm=(v * 1e15 if "m→fm" in u else v), mag=m)
                                 for n, v, u, m in inv],
               zero_parameter=dict(available=False,
                                   note="存貨無 0.1–1 fm 結構長度；ħc/Λ_IDCM=4.36 fm 差 11× 且比值非結構整數"))
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    emit("")
    emit(f"[JSON] {OUT}")
    emit(f"[LINES] {len(log)}")


if __name__ == "__main__":
    main()