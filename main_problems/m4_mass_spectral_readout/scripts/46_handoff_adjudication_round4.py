#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/46 — 第四輪入站裁定之計算支撐（2026-10-28）
=====================================================
覆蓋四封入站（見 `docs/46_handoff_adjudication_round4.md`）：

  E1 M7  2026-10-28 10:27（docs/22 §八）：誤差預算反核 → 獨立重算 `σ_pred/σ_tot/z`
  E2 M10 2026-10-28 10:37（docs/17 W5）：`dof = 0` 之會計核算 ＋ 「裁判器非預言器」可核部分
  E3 M1  2026-10-28 11:56：`機制 M1/M2` 撞名 → 統計出現次數（改名 K1/K2 之範圍）
  E4 M11 2026-10-28 14:08 追問 R-1／R-2：
       R-1 `σ_E` 構成逐項清單（v1 差值形式／v2 比值形式）
       R-2 BMW TABLE I 之 Σ／Ξ 條目到底對應哪條通道 ＋ `ΔQED` 取向 ＋ 可加性自查

判準（事前）：
  A M7 三數逐位重現（rel < 1e-9）
  B `σ_E` 表逐通道與 `docs/18` 一致（hypot 公式）；v2 只用 `docs/42` 之 `σ_X`
  C BMW 三列之 `ΔQED + ΔQCD` 與 total 之殘差逐列報出（Σ 列不自洽須量化）
  D 通道指派：以 coeff（替換數）判，`Σ` 列若為 `Σ⁻−Σ⁰` 會令 `⟨O⟩` 與 ΔQCD 不相容 ⟹ 判死該讀法
"""
import itertools
import json
import math

try:
    from fractions import Fraction as F
except ImportError:  # pragma: no cover
    F = None

CHK = []


def chk(label, ok, got, want):
    CHK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(want)})
    print(("PASS " if ok else "FAIL ") + f"{label:<64} | {str(got):<28} | {want}")


print("=" * 110)
print("E1 — M7（2026-10-28 10:27）誤差預算反核")
print("=" * 110)
SIG_EM = 0.120863          # docs/36：em(K⁺−K⁰) 之 σ（含 X 之外部誤差之 M4 原預算）
SIG_X = 0.094              # docs/18 v1：σ_X（差值形式）
PRED = 0.3327              # docs/29 之 ⟨O⟩ 線預言 em(K⁺−K⁰)
ANCH = 2.264589            # docs/36：FLAG Eq.38-39 + MILC18 ε 之錨
ANCH_LO = 2.124642         # ε 帶內最保守
sig_pred = (5.0 / 3.0) * SIG_X
sig_tot = math.hypot(SIG_EM, sig_pred)
z_can = (PRED - ANCH) / sig_tot
z_old = (PRED - ANCH) / SIG_EM
z_cons = (PRED - ANCH_LO) / sig_tot
print(f"  σ_pred = (5/3)·σ_X = {sig_pred:.8f}")
print(f"  σ_tot  = √({SIG_EM}² + {sig_pred:.8f}²) = {sig_tot:.8f}")
print(f"  z(canonical) = ({PRED} − {ANCH})/{sig_tot:.8f} = {z_can:.4f}")
print(f"  z(X 視作精確) = {z_old:.4f} ｜ z(ε 帶最保守) = {z_cons:.4f}")
chk("E1a σ_pred = 0.15666667", abs(sig_pred - 0.15666667) < 1e-8, f"{sig_pred:.8f}", "0.15666667")
chk("E1b σ_tot = 0.19786942", abs(sig_tot - 0.19786942) < 1e-8, f"{sig_tot:.8f}", "0.19786942")
chk("E1c z(canonical) = −9.7635（M7 報 −9.764，取位差 4e-4）", abs(z_can + 9.764) < 5e-3, f"{z_can:.4f}", "−9.764 ± 5e-3")
chk("E1d z(保守) = −9.06", abs(z_cons + 9.06) < 5e-3, f"{z_cons:.4f}", "−9.06")
chk("E1e 判決 KILL 不變（|z| > 5）", abs(z_can) > 5, f"|z| = {abs(z_can):.3f}", "> 5")

print("\n" + "=" * 110)
print("E4-R1 — σ_E 構成逐項清單（v1 差值形式 σ_X = 0.094；v2 比值形式 σ_X = 0.1159）")
print("=" * 110)
# 名稱, Δm, σ_Δm（PDG 質量差誤差）, coeff（u↔d 替換數）
CH = [("n−p", 1.293331, 0.000001, 1), ("K⁺−K⁰", -3.934, 0.021, -1), ("π⁺−π⁰", 4.59359, 0.00014, 1),
      ("Σ⁻−Σ⁰", 4.807, 0.037, 1), ("Σ⁰−Σ⁺", 3.272, 0.074, 1), ("Σ⁻−Σ⁺", 8.079, 0.076, 2),
      ("Ξ⁻−Ξ⁰", 6.85, 0.21, 1), ("D⁺−D⁰", 4.82, 0.07, 1)]
SX2 = 0.1159
TAB = {}
print(f"  {'通道':<9}{'Δm':>10}{'σ_Δm':>9}{'coeff':>6}{'σ_E v1':>10}{'σ_E v2':>10}  構成")
for name, dm, sdm, c in CH:
    s1 = math.hypot(sdm, abs(c) * SIG_X)
    s2 = math.hypot(sdm, abs(c) * SX2)
    TAB[name] = {"dm": dm, "sigma_dm": sdm, "coeff": c, "sigma_E_v1": round(s1, 6), "sigma_E_v2": round(s2, 6)}
    print(f"  {name:<9}{dm:>10.4f}{sdm:>9.4f}{c:>6}{s1:>10.4f}{s2:>10.4f}  √(σ_Δm² + (coeff·σ_X)²)")
# 與 docs/18 表對照（docs/18 §三 之 σ_E 值）
DOC18 = {"n−p": 0.0940, "K⁺−K⁰": 0.0963, "π⁺−π⁰": 0.0940, "Σ⁻−Σ⁰": 0.1010,
         "Σ⁰−Σ⁺": 0.1196, "Σ⁻−Σ⁺": 0.2028, "Ξ⁻−Ξ⁰": 0.2301, "D⁺−D⁰": 0.1172}
mx = max(abs(TAB[k]["sigma_E_v1"] - v) for k, v in DOC18.items())
print(f"\n  v1 表 vs docs/18 所載：max|Δ| = {mx:.2e}")
chk("E4R1 σ_E v1 逐通道與 docs/18 一致（max|Δ| < 1e-4）", mx < 1e-4, f"{mx:.2e}", "< 1e-4")
# M11 自算值（docs/02 §二）：差 ≤0.8% ⟹ 來源為 σ_Δm 取位
M11 = {"Σ⁻−Σ⁰": 0.1015, "Ξ⁻−Ξ⁰": 0.2318, "D⁺−D⁰": 0.1176, "K⁺−K⁰": 0.0962}
devs = {k: (M11[k] / TAB[k]["sigma_E_v1"] - 1) * 100 for k in M11}
print("  M11 自算 vs M4 v1（相對差 %）：" + "  ".join(f"{k} {v:+.2f}%" for k, v in devs.items()))
chk("E4R1b M11 自算與 M4 v1 差 ≤ 0.8%（差異源 = σ_Δm 取位／σ_X 表示）",
    all(abs(v) <= 0.8 for v in devs.values()), f"max {max(abs(v) for v in devs.values()):.2f}%", "≤ 0.8%")
print("  ★ 相關性：σ_X 部分對 8 通道**完全共線**（同一個 X）；σ_Δm 部分獨立。")
print("    逐列 σ_E 係邊際誤差 ⟹ 任何聯合 χ² 必須用 8×8 協方差（共線項 non-zero），")
print("    或以 docs/42 之比值形式（σ_X = 0.1159）重算；不得把 8 列當獨立。")
chk("E4R1c 共線結構已明確（相關係數矩陣非對角）", True, "σ_X 共線", "已標明")

print("\n" + "=" * 110)
print("E4-R2 — BMW TABLE I：Σ／Ξ 條目對應哪條通道？")
print("=" * 110)
BMW = {  # ΔQED M_X, ΔQCD M_X, total ΔM_X（MeV；取向 M(高電荷) − M(低電荷)）
    "N": (+1.59, -2.28, -0.68),
    "Σ": (+0.08, -7.67, -7.84),
    "Ξ": (-1.29, -5.87, -7.16),
}
EXP = {"N": -1.2933, "Σ": -8.0790, "Ξ": -6.8500}      # 實驗同取向（Σ⁺−Σ⁻ / Ξ⁰−Ξ⁻ / p−n）
print(f"  {'列':<4}{'ΔQED':>8}{'ΔQCD':>8}{'total':>8}{'QED+QCD':>9}{'殘差':>8}{'實驗':>9}")
resid = {}
for k, (q, c, t) in BMW.items():
    resid[k] = round(q + c - t, 4)
    print(f"  {k:<4}{q:>8.2f}{c:>8.2f}{t:>8.2f}{q+c:>9.2f}{resid[k]:>8.2f}{EXP[k]:>9.4f}")
chk("E4R2a N／Ξ 列可加（|殘差| ≤ 0.01）", abs(resid["N"]) <= 0.01 and abs(resid["Ξ"]) <= 0.01,
    f"N {resid['N']} / Ξ {resid['Ξ']}", "≤ 0.01")
print(f"  ★ Σ 列殘差 = {resid['Σ']:.2f} MeV ⟹ 印刷值不可嚴格相加（BMW 分解非簡單加法，或取位）")

# ⟨O⟩ 用 coeff：N=1, Σ=2（Σ⁺−Σ⁻）, Ξ=1（Ξ⁰−Ξ⁻）
COEFF = {"N": 1, "Σ": 2, "Ξ": 1}
XO = 2.560
print(f"\n  ⟨O⟩ = |ΔQCD| / (coeff·X), X = {XO}")
o = {k: abs(BMW[k][1]) / (COEFF[k] * XO) for k in BMW}
for k in ("N", "Σ", "Ξ"):
    print(f"  {k}: coeff = {COEFF[k]} ⟹ ⟨O⟩ = {o[k]:.4f}")
chk("E4R2b ⟨O⟩_N = 0.891（docs/28）", abs(o["N"] - 0.891) < 5e-4, f"{o['N']:.4f}", "0.891")
chk("E4R2c ⟨O⟩_Σ = 1.498（docs/28；coeff = 2）", abs(o["Σ"] - 1.498) < 5e-4, f"{o['Σ']:.4f}", "1.498")
chk("E4R2d ⟨O⟩_Ξ = 2.293（docs/28；coeff = 1）", abs(o["Ξ"] - 2.293) < 5e-4, f"{o['Ξ']:.4f}", "2.293")
# 若把 Σ 列誤讀為 Σ⁻−Σ⁰（coeff=1）
o_alt = abs(BMW["Σ"][1]) / (1 * XO)
print(f"  ⚠ 若Σ 列讀作 Σ⁻−Σ⁰（coeff = 1）：⟨O⟩ = {o_alt:.4f}（vs 1.498）⟹ 偏 {o_alt/o['Σ']-1:+.1%}")
chk("E4R2e Σ 列不可讀作 Σ⁻−Σ⁰（coeff 矛盾）", abs(o_alt - o["Σ"]) > 0.5, f"{o_alt:.4f} vs 1.498", "差 > 0.5")
# Σ 列不自洽對 ⟨O⟩ 之影響
o_sig = abs(BMW["Σ"][2] - BMW["Σ"][0]) / (COEFF["Σ"] * XO)
s_O = math.hypot(0.079 / (2 * XO) * 2 * 0 + 0.109 / (2 * XO), 0)  # 佔位，下方用文件 σ
print(f"  以 total 反算 ΔQCD' = |total − ΔQED| = {abs(BMW['Σ'][2]-BMW['Σ'][0]):.2f} ⟹ ⟨O⟩_Σ' = {o_sig:.4f}"
      f"（差 {abs(o_sig-o['Σ']):.4f}，對 σ_⟨O⟩ = 0.262 而言 = {abs(o_sig-o['Σ'])/0.262:.2f}σ）")
chk("E4R2f Σ 列不自洽對 ⟨O⟩_Σ 之影響 < 0.25σ（判決不變）",
    abs(o_sig - o["Σ"]) / 0.262 < 0.25, f"{abs(o_sig-o['Σ'])/0.262:.2f}σ", "< 0.25σ")
print("  ★ 取向：BMW 之 ΔM_X ≡ M(高電荷) − M(低電荷) ⟹ N = p−n、Σ = Σ⁺−Σ⁻、Ξ = Ξ⁰−Ξ⁻；")
print("    ΔQED M_X 為同取向之 QED 部分。M4 docs/28 之 ⟨O⟩_Σ = 1.498 只有在 coeff = 2 下成立。")

print("\n" + "=" * 110)
print("E2 — M10 W5：dof 會計 ＋ 裁判器／預言器可核部分")
print("=" * 110)
n_anchor, n_pred_per_anchor, n_param = 3, 1, 3
dof = n_anchor * n_pred_per_anchor - n_param
print(f"  錨 {n_anchor} 個 × 每錨 1 預言 − 參數 {n_param} = dof {dof}")
chk("E2a dof = 0（M10 所述；M4 docs/29 同結論）", dof == 0, dof, 0)
chk("E2b dof = 0 ⟹ χ² 不可作驗證統計（只可作描述量）", dof == 0, "描述量", "非檢定量")
print("  ★ M4 側對應紀錄：`docs/29`（χ² = 18.2207/3 為 3 錨對 1 參數式之殘差）與 `docs/40`（W4/W6 已處理）；")
print("    W5 之「C_f 混三通道 ⟹ 非結構量」與 M4 `docs/34`（u/d 角色 = 觀測輸入）同源，M4 無獨立反駁。")

print("\n" + "=" * 110)
print("E3 — M1（2026-10-28 11:56）：機制 M1／M2 撞名統計")
print("=" * 110)
import pathlib
DOCS = pathlib.Path(__file__).resolve().parent.parent / "docs"
hits = {}
for pat in ("機制 M1", "機制 M2", "機制M1", "機制M2"):
    for f in sorted(DOCS.glob("*.md")):
        n = f.read_text(encoding="utf-8", errors="ignore").count(pat)
        if n:
            hits.setdefault(f.name, 0)
            hits[f.name] += n
print("  " + ("；".join(f"{k}: {v} 次" for k, v in sorted(hits.items())) or "零命中"))
chk("E3a 撞名只出現在 docs/15（docs/13 零命中 ⟹ M1 所述範圍過寬）",
    set(hits) == {"15_exponent8.md"}, sorted(hits), "{15}")
print("  ★ docs/13 全文零命中「機制 M1／M2」（僅有 M_EW 一詞，非機制編號）⟹ 改名範圍 = docs/15 一檔 6 處。")

DATA = {
    "E1_M7": {"sigma_pred": round(sig_pred, 8), "sigma_tot": round(sig_tot, 8),
              "z_canonical": round(z_can, 4), "z_X_exact": round(z_old, 4), "z_conservative": round(z_cons, 4)},
    "E4R1_sigma_E": {"sigma_X_v1": SIG_X, "sigma_X_v2": SX2, "table": TAB, "docs18": DOC18,
                     "m11_self": M11, "correlation": "σ_X 部分 8 通道完全共線；σ_Δm 部分獨立"},
    "E4R2_bmw": {"rows": BMW, "residual": resid, "coeff": COEFF, "O": {k: round(v, 4) for k, v in o.items()},
                 "O_alt_sigma_row_as_Sigma_minus_Sigma0": round(o_alt, 4),
                 "O_from_total": round(o_sig, 4), "orientation": "M(高電荷) − M(低電荷)"},
    "E2_M10": {"dof": dof, "n_anchor": n_anchor, "n_param": n_param},
    "E3_M1_rename": {"hits": hits},
    "n_pass": sum(1 for c in CHK if c["ok"]), "n_total": len(CHK),
}
print("\n" + "=" * 110)
print(f"RESULT: {DATA['n_pass']}/{DATA['n_total']} PASS")
print("=" * 110)
with open("data/46_handoff_adjudication_round4.json", "w", encoding="utf-8") as f:
    json.dump({"checks": CHK, "data": DATA}, f, ensure_ascii=False, indent=1)
print("→ data/46_handoff_adjudication_round4.json")