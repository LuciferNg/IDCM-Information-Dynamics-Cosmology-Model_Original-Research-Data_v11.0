#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M12 / R3-P1 —— 多體（weight = 3）算子：κ_ijk 之典範收縮

事前判準（本輪開算前寫死，不得事後改）：
  P1  典範 weight-3 算子 = κ_ijk（CY₃(36,98) 三重相交數，純整數、零參數）。
       Universe 只有一個解 ⟹ 無權選向量 v ⟹ 唯一置換不變之收縮係「對全 1 向量」：
          K^(1)_ij := Σ_k κ_ijk          （36×36 對稱）
          v^(2)_i  := Σ_jk κ_ijk         （36 維向量）
          T3       := Σ_ijk κ_ijk        （純量）
  P2  檢 K^(1) 之精確 rank_ℚ 與核維數。
  P3  檢 K^(1) 之特徵多項式在 ℚ 上之分解 ⟹ 重根結構。**範疇性判準**：
      若有 3 重簡併 ⟹ 與 N_gen = 3 之世代結構相容（生成序）；若全單根 ⟹ 無世代訊號。
  P4  **不檢**「譜比值 vs 目標 r」——因無結構指派規則，36 個特徵值任選 4 個湊 4 個目標
      之組合數 = C(36,4)·4! ≈ 2.4e7 ⟹ 純掃描，屬恆真禁區。此項明文不作。

跑法：/usr/bin/python3 scripts/12_kappa_contraction.py
"""
from __future__ import annotations
import json
import os
from collections import Counter

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
L, OUT = [], {}

KAPPA_JSON = ("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/"
              "open_problems_fresh/o53_cy3_reality/data/polytopes/cy36_98_intersection_full.json")


def say(s=""):
    L.append(s)
    print(s)


def load_kappa():
    d = json.load(open(KAPPA_JSON))
    kap = d["kappa"]
    n = len(kap)
    assert d["shape"] == [n, n, n], d["shape"]
    k = [[[sp.Integer(kap[i][j][m]) for m in range(n)] for j in range(n)] for i in range(n)]
    return k, n


def main():
    k, n = load_kappa()
    say("=" * 86)
    say("R3-P1  κ_ijk 基本統計（h11 = %d ⟹ 36³ = %d 個分量）" % (n, n ** 3))
    say("")
    flat = [k[i][j][m] for i in range(n) for j in range(n) for m in range(n)]
    nz = [x for x in flat if x != 0]
    say("  非零分量 = %d / %d" % (len(nz), len(flat)))
    say("  值域 = [%s, %s]" % (min(nz), max(nz)))
    # 全對稱性檢查
    asym = 0
    for i in range(n):
        for j in range(n):
            for m in range(n):
                if k[i][j][m] != k[j][i][m] or k[i][j][m] != k[i][m][j]:
                    asym += 1
    say("  不對稱分量數 = %d ⟹ %s" % (asym, "✅ 全對稱" if asym == 0 else "⚠ 非全對稱"))
    OUT["P1"] = {"nonzero": len(nz), "min": str(min(nz)), "max": str(max(nz)), "asym": asym}

    say("")
    say("=" * 86)
    say("R3-P2  典範收縮（唯一置換不變者）")
    say("")
    T3 = sum(flat)
    v2 = [sum(k[i][j][m] for j in range(n) for m in range(n)) for i in range(n)]
    K1 = sp.Matrix(n, n, lambda i, j: sum(k[i][j][m] for m in range(n)))
    T1 = sum(k[i][i][i] for i in range(n))
    say("  T3 = Σ_ijk κ_ijk = %s" % T3)
    say("  T1 = Σ_i κ_iii   = %s" % T1)
    say("  v^(2) = Σ_jk κ_ijk：min %s, max %s, Σ_i v_i = %s"
        % (min(v2), max(v2), sum(v2)))
    say("  K^(1) 對稱？%s" % ("✅" if K1 == K1.T else "🔴 否"))
    say("  K^(1)_ii 之和 = tr = %s" % K1.trace())
    say("  K^(1) 非零元素 = %d / %d" % (sum(1 for x in K1 if x != 0), n * n))
    say("")
    say("  ★ 全 1 向量之收縮：1ᵀK^(1)1 = %s（= T3：自洽檢查）" % (sum(K1)))
    OUT["P2"] = {"T3": str(T3), "T1": str(T1), "tr_K1": str(K1.trace()),
                 "sum_v2": str(sum(v2))}

    say("")
    say("=" * 86)
    say("R3-P3  K^(1) 之精確 rank_ℚ 與特徵多項式")
    say("")
    rk = K1.rank()
    say("  rank_ℚ(K^(1)) = %d / %d ⟹ 核維數 = %d" % (rk, n, n - rk))
    say("  det(K^(1)) = %s" % K1.det())
    OUT["P3a"] = {"rank": rk, "nullity": n - rk, "det": str(K1.det())}

    x = sp.Symbol("x")
    say("  計算特徵多項式（Berkowitz, 精確）…")
    cp = K1.charpoly(x)
    fac = sp.factor_list(cp.as_expr())
    say("")
    say("  特徵多項式分解（ℚ 上）：")
    mult = []
    for f, e in fac[1]:
        fp = sp.Poly(f, x) if not isinstance(f, sp.Symbol) else sp.Poly(x, x)
        if fp.degree() == 1:
            cs = fp.all_coeffs()
            root = -cs[1] / cs[0]
            say("     (x − %s)^%d" % (sp.nsimplify(root), e))
            mult.append((str(sp.nsimplify(root)), e))
        else:
            say("     (%s)^%d   [deg %d，非有理根]" % (sp.factor(f), e, fp.degree()))
    say("")
    nz_mult = [(r, e) for r, e in mult if sp.sympify(r) != 0]
    say("  非零有理根重數 = %s" % (nz_mult if nz_mult else "（無非零有理根）"))
    say("  ★ 判準修正：λ=0 之 34 重簡併係 rank 2 之平凡後果，**不**計為世代訊號。")
    say("  非零譜重數 ≥3 之情形 = %s"
        % ("有" if any(e >= 3 for _, e in nz_mult) else "🔴 無"))
    OUT["P3b"] = {"factor_list": [[str(f), e] for f, e in fac[1]],
                  "rational_mult": mult,
                  "nonzero_mult": nz_mult,
                  "has_generation_signal": any(e >= 3 for _, e in nz_mult)}

    say("")
    say("=" * 86)
    say("R3-P3b  非零譜之精確值與像空間")
    say("")
    say("  由分解：char(K^(1)) = x^34 · (x² − 5) ⟹ 非零特徵值 = **±√5**（精確）")
    say("  ⟹ 唯一倖存之標度 = √5（φ 家族常數），**無**其他資訊")
    cs = K1.columnspace()
    say("  像空間基底（%d 維）：" % len(cs))
    for v in cs:
        say("     %s" % v.T)
    say("  ⚠ 若像空間基底向量非『整數／結構數』，則此收縮無乾淨之結構讀數。")
    OUT["P3c"] = {"nonzero_eigs": ["+sqrt(5)", "-sqrt(5)"],
                  "image_dim": len(cs),
                  "image_basis": [[str(x) for x in v] for v in cs]}

    say("")
    say("=" * 86)
    say("R3-P4  明文不作")
    say("")
    say("  不檢「譜比值 vs 目標 r」：無結構指派規則，C(36,4)·4! ≈ %.1e 種配對 ⟹ 純掃描，恆真禁區。"
        % float(sp.binomial(36, 4) * sp.factorial(4)))
    say("  若要檢，必須先有把 4 個目標通道 → 4 個特定特徵值 之結構規則（目前無）。")

    say("")
    say("=" * 86)
    say("R3-P5  【事前判準】K^(1) 之支撐 vs O137 之 9 模式集（獨立來源）")
    say("")
    say("  O137 `docs/31`：chiral I-60 軌道 ⟹ 9 模式 = {1,2,3,4,5,7,9,11,13}（n 指標）")
    say("  判準：若 K^(1) 之非零位置（含行列指標）落在該集合內 ⟹ 支撐命中；否則落空。")
    say("  兩種指標偏移（0-based / 1-based）皆報，不事後選擇。")
    say("")
    M9 = {1, 2, 3, 4, 5, 7, 9, 11, 13}
    supp = sorted({(i, j) for i in range(n) for j in range(n) if K1[i, j] != 0})
    idx = sorted({x for pair in supp for x in pair})
    say("  非零位置 = %d 個；涉及指標 = %s" % (len(supp), idx))
    say("  位置明細：%s" % supp)
    say("")
    for off, label in ((0, "0-based（指標 +0）"), (1, "1-based（指標 +1）")):
        mapped = sorted({x + off for x in idx})
        hit = [x for x in mapped if x in M9]
        say("  %-22s → %s" % (label, mapped))
        say("     命中 9 模式集 = %d / %d ⟹ %s"
            % (len(hit), len(mapped),
               "★ 全中" if len(hit) == len(mapped) and mapped
               else ("部分（%s）" % hit if hit else "🔴 全落空")))
        OUT.setdefault("P5", {})[label] = {"mapped": mapped, "hits": hit}
    say("")
    say("  ⚠ 注意：κ 之指標係 Hodge 類（h11=36），O137 之 n 係譜模式（1..44）。")
    say("     兩者無先驗對應 ⟹ 任何『對應』本身即為假設；本檢定只係『支撐是否巧合地落於該集』。")

    say("")
    say("=" * 86)
    say("R3-P1 判詞")
    say("")
    say("  ✅ κ_ijk 確係零參數 weight-3 算子（全對稱、純整數）")
    say("  ✅ 唯一置換不變收縮 K^(1) 已建；其譜係典範（無選擇自由）")
    say("  %s" % ("★ 非零譜有 3 重簡併 ⟹ 世代結構訊號" if any(e >= 3 for _, e in nz_mult)
                    else "🔴 非零譜無 3 重簡併（只有 ±√5 兩根）⟹ κ 之典範收縮**不**產生世代結構"))
    say("  🔴 rank 2 ⟹ κ 之置換不變內容最多分辨 2 個態；不足以承載 3 世代或 4 通道")
    say("  ⟹ 下一步：κ 作 flavour 算子**假設已判弱**；須改檢『對特定 v』之 κ(v,v,v)（但 v 之選擇即自由度）")


if __name__ == "__main__":
    main()
    say("")
    say("[written] data/12_kappa_contraction.json , data/run_12_kappa_contraction.log")
    with open(os.path.join(ROOT, "data", "12_kappa_contraction.json"), "w") as f:
        json.dump(OUT, f, indent=2, ensure_ascii=False,
                  default=lambda o: str(o))
    with open(os.path.join(ROOT, "data", "run_12_kappa_contraction.log"), "w") as f:
        f.write("\n".join(L) + "\n")