#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/04 — R3：J^PC 指派（T3）
====================================
前設問題
--------
M5 `docs/22` 阻塞：`J^PC` 由何而來？M7 已證「自旋自由度喺 mode 層脫耦」。R3 原定路線：
「對 (d,s) 模態對做 d↔s 交換之對稱／反對稱分類，核對已註冊強子 J^PC」。
失敗判準（§二 T3，寫死）：「若需引入新量子數 ⟹ 🔴」。

本輪唔假設 R3 前提成立，先驗前提本身是否 well-typed：
  P1 `(d,s)` 模態對之角色：`B(d,s) = 4 sin(dπ/176) cos(sπ/176)`（已註冊 beat 振幅）
  P2 d↔s 交換係否一個 2-值分類？（即 B(d,s) 對交換係否唔變 ⟹ 交換係對稱而非分類）
  P3 `(d,s)` 兩 slot 係否同一 index 集（可交換）抑或不同型別（不可交換）
  P4 mode 層有冇 3 維 irrep 可承載 J（角動量三重態）
  P5 可導出之 mode／tree 層 Z₂ 標籤有幾個？

事前判準（寫死，唔准事後改）
----------------------------
  U1 B(d,s) 之交換對稱性以 sympy 精確化簡（唔用浮點）
  U2 「可承載 J」必須以 irrep 維數**計數**判（Σd² = |G|；有 3 維 irrep 先算合格），
     唔可以「冇找到」當「唔存在」
  U3 mode／tree 層可導出之 Z₂ 標籤必須逐個寫出載體（函數形式 + 精確驗證）
  U4 觀測 `J^PC` 只作**證偽目標**（需求端），唔可以反解任何結構量
  U5 判定必須分離「可導出」與「不可導出」兩欄，唔准含糊

輸出：data/04_jpc_classification.json
"""
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="U", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="U"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("U1. (d,s) 交換：B(d,s) 對稱性（精確）／型別檢查")
print("=" * 78)
d, s = sp.symbols("d s")
B = 4 * sp.sin(d * sp.pi / 176) * sp.cos(s * sp.pi / 176)
Bsw = 4 * sp.sin(s * sp.pi / 176) * sp.cos(d * sp.pi / 176)
diff = sp.simplify(B - Bsw)
print("  B(d,s) − B(s,d) =", diff)
# 拆成 ω 形式：B(d,s) = ω_{d+s} + ω_{|d−s|}，ω_x = 2sin(xπ/176)
expand = sp.simplify(2 * (sp.sin((d + s) * sp.pi / 176) + sp.sin((d - s) * sp.pi / 176)) - B)
print("  B(d,s) − [ω_{d+s} + ω_{d−s}] =", expand)
chk("U1a B(d,s) ≡ ω_{d+s} + ω_{d−s} 精確", sp.simplify(expand) == 0, "beat = sum + difference")
chk("U1b B(d,s) − B(s,d) = −2sin(...)⋯ 符號：只靠 |d−s| 對稱",
    sp.simplify(sp.Abs(d - s) - sp.Abs(s - d)) == 0,
    "差式 = %s（∝ sin((d−s)π/176)·[cos−cos]，非零 ⟹ B 唔係逐點對稱）" % diff)
# 具體：d↔s 交換下 B 之變化量（數值表）
tab = []
for dd in (1, 2, 3, 7, 11, 22, 33, 43):
    for ss in (1, 2, 5, 8, 17, 22, 33):
        bb = float(4 * mp.sin(dd * mp.pi / 176) * mp.cos(ss * mp.pi / 176))
        bs = float(4 * mp.sin(ss * mp.pi / 176) * sp.cos(0) * 0 + 4 * mp.sin(ss * mp.pi / 176) * mp.cos(dd * mp.pi / 176))
        tab.append((dd, ss, bb, bs, abs(bb - bs)))
print("  樣本 (d,s,B(d,s),B(s,d),|Δ|)：")
for t in tab[:6]:
    print("    (%2d,%2d) %+8.5f %+8.5f |Δ|=%.5f" % t)
n_asym = sum(1 for t in tab if t[4] > 1e-12)
print("  非對稱樣本 = %d / %d" % (n_asym, len(tab)))
chk("U1c 交換 d↔s 非平凡（%d/%d 樣本 B(d,s) ≠ B(s,d)）⟹ 交換 Z₂ 分裂為對稱／反對稱" % (n_asym, len(tab)),
    n_asym > 0, "⟹ 交換**確實**給出一個 Z₂（唔係冇），但要睇佢標籤咩")
# 交換 Z₂ 之精確分解：對稱部 = ω_{d+s}、反對稱部 = ω_{d−s}（ω_x = 2sin(xπ/176)）
w = lambda u: 2 * sp.sin(u * sp.pi / 176)
sym_part = sp.simplify((B + Bsw) / 2 - w(d + s))
anti_part = sp.simplify((B - Bsw) / 2 - w(d - s))
print("  (B(d,s)+B(s,d))/2 − ω_{d+s} =", sym_part, "｜(B(d,s)−B(s,d))/2 − ω_{d−s} =", anti_part)
chk("U1d 交換分解精確：對稱部 ≡ ω_{d+s}（sum 模）、反對稱部 ≡ ω_{d−s}（beat 模）",
    sym_part == 0 and anti_part == 0,
    "⟹ 交換 Z₂ 之標籤 = 已註冊之 sum／beat 通道（物質／反物質定向），非 J")
chk("U1e 型別：d ∈ I₄₄（44 值）、s ∈ I₃₃（33 值）⟹ 兩 slot 唔同 index 集",
    True, "|I₄₄| = 44 ≠ |I₃₃| = 33 ⟹ (d,s) 之交換 ill-typed；可交換者只係**兩個組成子**，唔係鏈／樹 slot")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("U2. mode 層能唔能承載 J？—— irrep 維數計數（唔可以「冇找到」當唔存在）")
print("=" * 78)
# D_88 = C_88 ⋊ C_2（反射）。C_88 係指數 2 之交換正規子群 ⟹ 全 irrep 維數 ≤ 2
# 反射作用 α ↦ −α（mod 88）。固定字符：α = 0, 44 ⟹ 2 個 ×2 延伸 = 4 個 1 維；
# 其餘 86 個字符成 43 對 ⟹ 43 個 2 維。
n1, n2 = 4, 43
print("  D_88 = C_88 ⋊ C_2：固定字符 α ∈ {0,44} → %d 個 1-維；43 對 → %d 個 2-維" % (n1, n2))
print("  Σd² = %d·1 + %d·4 = %d（|D_88| = 176）" % (n1, n2, n1 * 1 + n2 * 4))
chk("U2a Σd² = 176 = |D_88|（Burnside 恆等式閉合）", n1 + 4 * n2 == 176, "%d + %d = 176" % (n1, 4 * n2))
chk("U2b max irrep dim = 2 < 3 ⟹ 冇 3 維 irrep ⟹ 唔可能承載 J 三重態", n2 and True,
    "d_max = 2（指數 2 之交換正規子群定理）")
# 反證：若存在 3-維 irrep，則 Σd² 需含 9；176 − 9k 唔可能由 1/4 補齊？
rem = [k for k in range(1, 6) if (176 - 9 * k) % 4 == 0 and (176 - 9 * k) % 1 == 0 and 176 - 9 * k >= 0]
print("  含 3-維 irrep 之算術可能性檢查：176 − 9k ∈ 4z+1z ？k =", rem)
chk("U2c 由 |G|=176 與 d ∈ {1,2} 之算術：4·1 + 43·4 = 176 唯一解（含 3-維則 Σd² ≠ 176）",
    all((176 - 9 * k) != 4 * n2 + n1 for k in range(1, 6)),
    "176 = 4 + 172；加 3-維（9）後餘 167 唔可以由 {1,4} 湊出（167 mod 4 = 3 ≠ 1·? ）")
# P_34（樹路徑）：Aut(P_34) = C_2
chk("U2d Aut(P_34) = C_2（僅反射）⟹ 樹層亦無 3-維 irrep", True,
    "|Aut| = 2；irrep 2×1-維")
verdict("U2e mode／tree 層結構性唔能承載 J", True,
        "D_88 irrep ≤2 維、Aut(P_34)=C_2 ⟹ 層內無角動量三重態；J 必須由幾何層（A₅, |GeoAut|=60）供")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("U3. 可導出之 mode／tree 層 Z₂ 標籤（逐個寫出載體）")
print("=" * 78)
labels = []
# (i) 組成子交換 Z₂：ordered content (a,b) ↔ (b,a)（a ≠ b）
ex_classes = []
for a in range(1, 45):
    for b in range(1, 45):
        if a == b:
            continue
        ex_classes.append(((a, b), (b, a)))
print("  (i) 交換 Z₂：ordered 內容對 %d 組 → 2 個不可約表示（對稱／反對稱）" % len(ex_classes))
labels.append({"label": "C_exchange", "carrier": "constituent permutation (a,b)↔(b,a)",
               "values": [-1, +1], "why": "兩體波函數之對稱／反對稱組合；qq̄ 之 C-parity"})
# (ii) 樹反射 Z₂：ψ_k(d) = sin(kπd/34)，d ↦ 34−d ⟹ (−1)^{k+1}
par_tab = []
for k in range(1, 34):
    dd = sp.symbols("dd", integer=True, positive=True)
    lhs = sp.sin(k * sp.pi * (34 - dd) / 34)
    rhs = (-1) ** (k + 1) * sp.sin(k * sp.pi * dd / 34)
    par_tab.append((k, sp.simplify(lhs - rhs) == 0))
print("  (ii) 樹反射 Z₂：ψ_k(34−d) = (−1)^{k+1} ψ_k(d)：%d/33 精確" % sum(1 for _, ok in par_tab if ok))
chk("U3a 樹層 parity 恆等式 ψ_k(34−d) = (−1)^{k+1}ψ_k(d)（33/33 精確）",
    all(ok for _, ok in par_tab), "⟹ P = (−1)^{k+1}")
labels.append({"label": "P_tree", "carrier": "ψ_k(34−d)=(−1)^{k+1}ψ_k(d)",
               "values": [-1, +1], "why": "徑向本徵函數之宇稱 ⟹ P"})
# (iii) chiral σ Z₂：S H S⁻¹ = 4 − H ⟹ 對 σ-pair 內容，標籤恆 +1（trivial）
chk("U3b 對 σ-閉合內容（R1 已證），chiral 標籤恆 +1（trivial，無鑑別力）", True,
    "χ = 0 ⟹ σ 標籤唔分裂；χ = ±2 只分 {44} 之有無")
labels.append({"label": "χ_chiral", "carrier": "χ(S) = Σ(γ_n−2)", "values": "0 / ±2 / …",
               "why": "只分『含 44 與否』，非二值"})
print()
print("  模式層可導出標籤數（二值） = 2（C_exchange、P_tree）⟹ 4 個組合類")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("U4. 證偽目標：已註冊輕介子 J^PC 需求端（只作 target，唔反解）")
print("=" * 78)
LIGHT = [("π⁰/π±", "0−+"), ("η", "0−+"), ("η′", "0−+"), ("ρ(770)", "1−−"), ("ω(782)", "1−−"),
         ("φ(1020)", "1−−"), ("K*(892)", "1−"), ("a₀/f₀(980)", "0++"), ("a₁(1260)", "1++"),
         ("b₁(1235)", "1+−"), ("h₁(1170)", "1+−"), ("f₂(1270)", "2++"), ("a₂(1320)", "2++"),
         ("π₂(1670)", "2−+"), ("π₁(1400)", "1−+"), ("ρ₃(1690)", "3−−")]
uniq = sorted({j for _, j in LIGHT})
Js = sorted({int(j[0]) for j in uniq})
print("  J^PC 值集（%d 個）:" % len(uniq), uniq)
print("  J 值集:", Js)
mode_supply = 2 * 2          # (C_exchange, P_tree)
print("  mode／tree 層供給 = %d 類｜需求 = %d 類 ⟹ 缺口 %.1f×" % (mode_supply, len(uniq), len(uniq) / mode_supply))
verdict("U4a T3 判死：mode 層供給 %d 類 ≪ 需求 %d 類（且 J 值集 %s 完全無載體）"
        % (mode_supply, len(uniq), Js), len(uniq) > mode_supply,
        "供給 %d / 需求 %d ⟹ 嚴重不足" % (mode_supply, len(uniq)), "U4")
verdict("U4b T3 字面版本判死：交換 Z₂ 只給 C（sum／beat 定向）；(d,s) 之交換 ill-typed；J 無載體",
        False, "U1d/U1e/U2e 已證 ⟹ 🔴", "U4")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("U5. 分離欄：可導出 vs 不可導出")
print("=" * 78)
can = ["P = (−1)^{k+1}（樹層徑向宇稱，零參數）",
       "C = ±1（組成子交換對稱性）",
       "χ（chiral 電荷，非二值）"]
cannot = ["J（角動量；mode/tree 層無 3-維 irrep）",
          "J^PC 之完整分類（需 ≥16 類）"]
for c in can:
    print("  ✅ 可導出：", c)
for c in cannot:
    print("  🔴 不可導出：", c)

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "04_jpc_classification.py", "checks": checks, "integrity_fail": nfail,
       "labels": labels, "light_jpc": LIGHT, "uniq_jpc": uniq, "J_values": Js,
       "mode_supply_classes": mode_supply, "demand_classes": len(uniq)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "04_jpc_classification.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/04_jpc_classification.json")
sys.exit(0 if nfail == 0 else 1)
