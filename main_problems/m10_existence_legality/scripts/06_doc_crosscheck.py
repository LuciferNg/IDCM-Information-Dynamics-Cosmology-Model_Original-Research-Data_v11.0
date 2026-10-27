#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/06 — 文檔數字交叉核（docs ↔ data，單向：data 為權威）

前設問題:
  M5 docs/11 與 M7 scripts/04 之做法 —— 文檔引述之每個數字必須可回溯到 data/*.json。
  本腳本逐條核 docs/01–05 之關鍵數字；**期望值一律由 JSON 欄位即場計算**，
  腳本內不寫死任何數值（除格式），以免「腳本抄文檔、文檔抄腳本」之循環。

判準（事前寫死）:
  C1 每條 claim 之格式化字串必須在對應 doc 出現（Unicode 負號正規化後比對）
  C2 由 JSON 重算之關鍵合成量（x_n、佔比）與 doc 一致
  C3 doc 中所有 64-hex 封印必須是 data 內之封印（不得有手寫假 hash）
  C4 每份 doc 至少被核 5 條（防止空核）

執行：/usr/bin/python3 scripts/06_doc_crosscheck.py
"""
import json
import os
import re
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DOCS, DATA = os.path.join(ROOT, "docs"), os.path.join(ROOT, "data")
PASS, FAIL = [], []


def norm(t):
    return t.replace("\u2212", "-").replace("\u00d7", "x").replace("\u2009", "")


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


def e_(v, dp):
    """科學記號，指數去前導零，用於比對 doc 之「e−5」寫法。"""
    m, ex = f"{v:.{dp}e}".split("e")
    return f"{m}e{int(ex)}"


J = {}
for f in ("01_criteria_freeze", "02_sign_pattern", "03_precision_wall",
          "04_topology_gate", "05_blind_prediction", "07_environment",
          "08_nuclide_chart", "09_sphaleron"):
    with open(os.path.join(DATA, f + ".json"), encoding="utf-8") as fh:
        J[f] = json.load(fh)

DOC = {}
for f in ("01_criteria_freeze", "02_three_layer_classification",
          "03_precision_wall_correction", "04_topology_gate", "05_preregistration",
          "06_environment_completeness", "07_nuclide_chart",
          "08_sphaleron_relaxation"):
    with open(os.path.join(DOCS, f + ".md"), encoding="utf-8") as fh:
        DOC[f] = norm(fh.read())

print("=" * 78)
print("M10 — 文檔數字交叉核（docs ← data，單向）")
print("=" * 78)

CLAIMS = []


def cl(doc, s, why):
    CLAIMS.append((doc, norm(s), why))


# ---------------- docs/01 ← 01_criteria_freeze
j1 = J["01_criteria_freeze"]
cl("01_criteria_freeze", j1["C5_seal"]["sha256"][:16], "判據封印前綴")
q = j1["C2_quadrants"]
cl("01_criteria_freeze", f"{q['Q1 閘✅ 能✅']['x_MeV']:+.4f}", "Q1 n→peν x")
cl("01_criteria_freeze", f"{q['Q2 閘✅ 能❌']['x_MeV']:+.4f}", "Q2 ²H→pn x")
cl("01_criteria_freeze", f"{q['Q3 閘❌ 能✅']['x_MeV']:+.4f}", "Q3 p→π⁰e⁺ x")
cl("01_criteria_freeze", f"{q['Q4 閘❌ 能❌']['x_MeV']:+.1f}", "Q4 ppp̄ x")
e = j1["C3_environment"]
cl("01_criteria_freeze", e_(e["E_cmb_peak_MeV"], 3), "CMB 峰值能")
cl("01_criteria_freeze", e_(e["T_dissociation_K"], 3).replace("e", "e"), "解離溫度")
cl("01_criteria_freeze",
   f"{float(F(e['E_bind_2H_MeV']))/e['E_cmb_peak_MeV']:.2e}".split("e")[0],
   "門檻/CMB 量級")

# ---------------- docs/02 ← 02_sign_pattern
j2 = J["02_sign_pattern"]
for r in j2["A_reproduction"]["rows"]:
    cl("02_three_layer_classification", r["x_str"], f"重現閘 x {r['channel']}")
    cl("02_three_layer_classification", r["ratio_str"], f"重現閘 |x|/m {r['channel']}")
b = j2["B_three_layers"]
cl("02_three_layer_classification", f"{b['t']['ctau_over_fm']:.3f}", "t 之 cτ/1fm")
for key in ("b", "c", "s", "d"):
    cl("02_three_layer_classification", f"{b[key]['ctau_over_fm']:.2e}", f"{key} 之 cτ")
d = j2["D_single_bit"]
cl("02_three_layer_classification", f"{d['X_real']:+.4f}", "X_real")
cl("02_three_layer_classification", f"{d['x_p_swap']:+.4f}", "互換 x_p")
emc = j2["E_mc_null"]
cl("02_three_layer_classification", f"{emc['q_n_illegal_within_10pct_error']:.6f}", "E.1 q")
cl("02_three_layer_classification", f"{emc['q_p_legal_unconstrained_random']:.5f}", "E.2 q")
f6 = j2["F_idcm_side"]
cl("02_three_layer_classification", f"{f6['gamma3_over_gamma2']:.11f}", "γ₃/γ₂")
cl("02_three_layer_classification", f"{f6['z_vs_data']:+.2f}", "z vs 數據")
cl("02_three_layer_classification", f"{f6['X_pred_MeV']:.4f}", "X_pred")
cl("02_three_layer_classification", f"{f6['em_needed_pred']:.4f}", "所需 em（γ 路線）")
cl("02_three_layer_classification", f"{f6['em_needed_flag']:.4f}", "所需 em（FLAG 路線）")
cl("02_three_layer_classification", f"{f6['A_max_v2']}", "A_max v2")

# ---------------- docs/03 ← 03_precision_wall
j3 = J["03_precision_wall"]
for r in j3["G1_channels"]:
    cl("03_precision_wall_correction", f"{r['cancel_factor']:.3e}", f"相消倍數 {r['channel']}")
    cl("03_precision_wall_correction", f"{r['need_rel_precision']:.2e}", f"判號精度 {r['channel']}")
g2 = j3["G2_two_precisions"]
for k, v in g2["delta_star"].items():
    cl("03_precision_wall_correction", e_(v, 1), f"δ★ {k}")
g4 = j3["G4_neutron_budget"]
cl("03_precision_wall_correction", f"{g4['share_X']*100:.0f}%", "X 佔 |x| 327%")
cl("03_precision_wall_correction", f"{g4['share_em']*100:.0f}%", "em 佔 |x| 162%")
cl("03_precision_wall_correction", f"{g4['share_me']*100:.0f}%", "m_e 佔 |x| 65%")
cl("03_precision_wall_correction", f"{g4['X_precision_needed']*100:.1f}%", "X 所需 30.6%")
cl("03_precision_wall_correction", f"{g4['em_precision_needed']*100:.1f}%", "em 所需 61.8%")
cl("03_precision_wall_correction", f"{g4['A_max_v2']}", "A_max v2")
cl("03_precision_wall_correction", f"{g4['em_needed_gamma']:.4f}", "γ 路線所需 em")
cl("03_precision_wall_correction", f"{100*(g4['em_needed_gamma']/g4['A_max_v2']-1):.2f}%", "超界 11.76%")
cl("03_precision_wall_correction", f"{100*(g4['em_needed']/g4['A_max_v2']-1):.2f}%", "超界 1.38%")

# ---------------- docs/04 ← 04_topology_gate
j4 = J["04_topology_gate"]
h1 = j4["H1_spectral_no_go"]
cl("04_topology_gate", f"{h1['gamma11']-h1['gamma9']:.10f}", "γ₁₁−γ₉")
cl("04_topology_gate", f"{h1['gamma9']:.10f}", "γ₉")
cl("04_topology_gate", f"{h1['gamma11']:.10f}", "γ₁₁")
h3 = j4["H3"]
cl("04_topology_gate", f"{h3['N_d']} mod {h3['N_gen']} = {h3['N_d_mod_N_gen']}", "88 mod 3")
cl("04_topology_gate", f"{h3['anomaly_coefficient']} ≠ 0", "反常係數 3")

# ---------------- docs/05 ← 05_blind_prediction
j5 = J["05_blind_prediction"]
t = j5["J_t3H"]
cl("05_preregistration", f"{t['frozen']:.3f}", "³H 凍結值")
cl("05_preregistration", f"{t['correct']:.4f}", "³H 正確值")
cl("05_preregistration", f"{abs(t['frozen']-t['correct']):.4f}", "Δ 0.5301")
cl("05_preregistration", f"{abs(t['x_frozen']):.4f}", "x 凍結 −9.0119")
cl("05_preregistration", f"{abs(t['x_correct']):.4f}", "x 更正 8.4818")
cl("05_preregistration", f"{abs(t['x_correct'])/t['correct']:.3e}", "更正後 |x|/m")
cl("05_preregistration", j5["L_seal"], "預註冊封印")
cl("05_preregistration", f"{j5['K_min_rel']:.3e}", "外推族最小 |x|/m")
cl("05_preregistration", f"{j5['K_out_of_sample'][0]['rel']:.2e}", "Λ→pπ⁻ 之 |x|/m")


# ---------------- docs/06 ← 07_environment
j6 = J["07_environment"]
cl("06_environment_completeness", f"{j6['M1']['T_D_MeV_idcm']*1e3:.2f}", "T_D (IDCM η)")
cl("06_environment_completeness", f"{j6['M1']['y_idcm']:.3f}", "y = Q/kT")
cl("06_environment_completeness", f"{j6['M1']['ratio_to_BBN']:.3f}", "T_D/T_BBN")
cl("06_environment_completeness", f"{j6['M2']['rho_th_gcm3']:.3e}", "ρ_th")
cl("06_environment_completeness", f"{j6['M2']['gap_MeV']:.2f}", "W3 缺口")
cl("06_environment_completeness", f"{j6['M2']['mu_drip_MeV']:.2f}", "μ_e(drip)")
cl("06_environment_completeness", f"{j6['M3']['hw_c_1e12G_keV']:.3f}", "ħω_c @1e12G")
cl("06_environment_completeness", f"{j6['M3']['hw_c_1e15G_MeV']:.2f}", "ħω_c @1e15G")
cl("06_environment_completeness", f"{j6['M3']['B_for_078MeV_G']:.2e}", "B(0.78MeV)")
cl("06_environment_completeness", f"{j6['M3']['B_for_3H_Q_G']:.2e}", "B(³H Q)")
cl("06_environment_completeness", f"{j6['T_3H']['rel']:.3e}", "³H |x|/m")
cl("06_environment_completeness", f"{j6['T_3H']['ratio']:.1f}", "δ★/rel")

# ---------------- docs/07 ← 08_nuclide_chart
j7 = J["08_nuclide_chart"]
cl("07_nuclide_chart", str(j7["n_nuclides"]), "核素條數")
cl("07_nuclide_chart", str(j7["crosscheck"]["compared"]), "交叉核可比條數")
cl("07_nuclide_chart", f"{j7['crosscheck']['worst_keV']:.2f}", "最大偏差 keV")
for st in j7["stages"]:
    cl("07_nuclide_chart", str(st["n_legal"]), f"合法數 {st['channels'][-1]}")
cl("07_nuclide_chart", str(j7["stages"][-1]["A_max"]), "全表 A_max")
for st in j7["stages_measured_only"]:
    cl("07_nuclide_chart", str(st["n_legal"]), f"實測合法數 {st['channels'][-1]}")
cl("07_nuclide_chart", str(j7["stages_measured_only"][-1]["A_max"]), "實測 A_max")
cl("07_nuclide_chart", f"{100*j7['sf_contraction_measured']:.1f}%", "SF 收縮率")
cl("07_nuclide_chart", str(j7["tier_bins"]["III (<1e-4)"]), "Tier III 條數")
cl("07_nuclide_chart", str(j7["tier_bins"]["II (1e-4…1e-2)"]), "Tier II 條數")
cl("07_nuclide_chart", f"{j7['tier3_median']:.2e}", "Tier III 中位數")
cl("07_nuclide_chart",
   str(j7["stages_measured_only"][2]["n_legal"] - j7["stages_measured_only"][-1]["n_legal"]),
   "漏 SF 之誤判數")

# ---------------- docs/08 ← 09_sphaleron
j8 = J["09_sphaleron"]
cl("08_sphaleron_relaxation", f"{j8['E_sph_GeV']['1.9']:.0f}", "E_sph (B=1.9)")
cl("08_sphaleron_relaxation", f"{j8['E_sph_GeV']['1.5']/1e3:.2f}", "E_sph (B=1.5)")
cl("08_sphaleron_relaxation", f"{j8['E_sph_GeV']['2.5']/1e3:.2f}", "E_sph (B=2.5)")
tab8 = {r["T_GeV"]: r["supp"] for r in j8["suppression_table"]}
for T in (159, 140, 120, 100, 50, 20):
    cl("08_sphaleron_relaxation", f"{tab8[T]:.2e}", f"抑制 @{T} GeV")
cl("08_sphaleron_relaxation", f"{j8['T_for_1e30_supp']:.1f}", "T(1e−30)")
cl("08_sphaleron_relaxation", f"{j8['B_minus_L_required']:.4e}", "(B−L) 需求")
cl("08_sphaleron_relaxation", f"{j8['S_inst']:.1f}", "S_inst")
cl("08_sphaleron_relaxation", f"{j8['inst_rate']:.2e}", "e^-2S")

print(f"\n[C1] 逐條核 {len(CLAIMS)} 條 claim：")
bad = [(doc, s, why) for doc, s, why in CLAIMS if s not in DOC[doc]]
for doc, s, why in bad:
    print(f"    ❌ {doc}: 「{s}」({why})")
ck("C1 全部 claim 字串在 doc 出現", not bad, f"{len(CLAIMS)-len(bad)}/{len(CLAIMS)}")

print("\n[C2] 由 JSON 重算之合成量：")
re_x = g4["X"] - g4["em_needed"] - g4["m_e"]
ck("C2.1 x_n 可由 X − em − m_e 重算", abs(re_x - g4["x_n"]) < 1e-9,
   f"{re_x:.5f} vs JSON {g4['x_n']:.5f}")
ok_x = f"{re_x:.4f}".lstrip("+") in DOC["03_precision_wall_correction"]
ck("C2.2 重算之 x_n 在 doc 出現", ok_x, f"「{re_x:.4f}」")
ck("C2.3 三佔比之和 ≈ 5.44（= (X+em+m_e)/|x|）", abs(
    g4["share_X"] + g4["share_em"] + g4["share_me"] -
    (g4["X"] + g4["em_needed"] + g4["m_e"]) / g4["x_n"]) < 1e-9,
   "恆等式 self-consistent（佔比係相對 |x| 之和，非 100%）")

print("\n[C3] 封印完備性：")
seals = {j1["C5_seal"]["sha256"], j5["L_seal"]}
found = set(re.findall(r"\b[0-9a-f]{64}\b", " ".join(DOC.values())))
ck("C3.1 doc 中 64-hex 皆為 data 內封印", found <= seals,
   f"doc {len(found)} 個 / data {len(seals)} 個；多出 {sorted(found-seals)}")
ck("C3.2 兩個封印互異（判據 ≠ 預言）", j1["C5_seal"]["sha256"] != j5["L_seal"],
   f"判據 {j1['C5_seal']['sha256'][:12]}… / 預言 {j5['L_seal'][:12]}…")

print("\n[C4] 覆蓋度：")
per = {}
for doc, s, why in CLAIMS:
    per[doc] = per.get(doc, 0) + 1
for doc in sorted(per):
    print(f"    {doc:<34} {per[doc]:>3} 條")
ck("C4 每份 doc 至少 5 條", all(v >= 5 for v in per.values()) and len(per) == 8,
   f"{len(per)} 份 doc，最少 {min(per.values())} 條")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL:", FAIL)
print("=" * 78)

out = {"script": "06_doc_crosscheck.py", "role": "docs ↔ data 交叉核",
       "claims": len(CLAIMS), "per_doc": per,
       "seals_data": sorted(seals), "seals_in_docs": sorted(found),
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL),
       "failures": FAIL}
with open(os.path.join(DATA, "06_doc_crosscheck.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/06_doc_crosscheck.json 已寫")
sys.exit(0 if not FAIL else 1)