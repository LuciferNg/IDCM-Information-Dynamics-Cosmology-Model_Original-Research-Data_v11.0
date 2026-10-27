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
          "08_nuclide_chart", "09_sphaleron", "10_channel_completion",
          "11_falsification_nubase", "12_definition_d5", "13_primordial_test",
          "14_barrier_gate"):
    with open(os.path.join(DATA, f + ".json"), encoding="utf-8") as fh:
        J[f] = json.load(fh)

DOC = {}
for f in ("01_criteria_freeze", "02_three_layer_classification",
          "03_precision_wall_correction", "04_topology_gate", "05_preregistration",
          "06_environment_completeness", "07_nuclide_chart",
          "08_sphaleron_relaxation", "09_channel_completion",
          "10_falsification_nubase", "11_definition_d5",
          "12_primordial_label", "13_barrier_gate"):
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


# ---------------- docs/09 ← 10_channel_completion
j9c = J["10_channel_completion"]
cl("09_channel_completion", f"{j9c['n_measured']}", "實測池條數")
cl("09_channel_completion", f"{j9c['baseline']['n_legal']}", "基線合法數")
cl("09_channel_completion", f"{j9c['baseline']['A_max']}", "基線 A_max")
cl("09_channel_completion", f"{j9c['final']['n_legal']}", "𝒞_v2 合法數")
cl("09_channel_completion", f"{j9c['final']['A_max']}", "𝒞_v2 A_max")
cl("09_channel_completion",
   f"{100*(1-j9c['final']['n_legal']/j9c['baseline']['n_legal']):.1f}%", "總收縮率")
for r in j9c["marginal"]:
    cl("09_channel_completion", f"{r['legal_after']}", f"加入 {r['added']} 後之合法數")
    cl("09_channel_completion", f"{r['killed']}", f"{r['added']} 之殺數")
    cl("09_channel_completion", f"{100*r['contraction']:.1f}%", f"{r['added']} 之收縮率")
kc = j9c["kill_count_on_baseline"]
for name in ("SF_opt", "cluster", "2β⁻", "2β⁺/ECEC"):
    cl("09_channel_completion", f"{100*kc[name]/j9c['baseline']['n_legal']:.1f}%",
       f"{name} 之獨立殺傷率")
cl("09_channel_completion", f"{j9c['kill_count_on_baseline']['SF_opt']}", "SF_opt 獨立殺數")
cl("09_channel_completion", f"{j9c['kill_count_on_baseline']['cluster']}", "cluster 獨立殺數")
cl("09_channel_completion", f"{j9c['kill_count_on_baseline']['2β⁻']}", "2β⁻ 獨立殺數")
cl("09_channel_completion", f"{j9c['kill_count_on_baseline']['2β⁺/ECEC']}", "2β⁺ 獨立殺數")

# ---------------- docs/10 ← 11_falsification_nubase
j10c = J["11_falsification_nubase"]
tb = j10c["test_B"]
cl("10_falsification_nubase", f"{j10c['n_nubase']}", "NUBASE 基態條數")
cl("10_falsification_nubase", f"{j10c['n_parse_fail']}", "NUBASE 解析失敗數")
cl("10_falsification_nubase", f"{j10c['anchors']['3H']:.2f}", "³H T½ (yr)")
cl("10_falsification_nubase", e_(j10c["anchors"]["209Bi"], 3), "²⁰⁹Bi T½ (yr)")
cl("10_falsification_nubase", f"{j10c['anchors']['14C']/1000:.2f}", "¹⁴C T½ (ky)")
cl("10_falsification_nubase", f"{j10c['P_size']}", "P 條數")
cl("10_falsification_nubase", f"{tb['S_size']}", "S(1e15) 條數")
cl("10_falsification_nubase", f"{tb['SP_size']}", "S\\P 條數")
cl("10_falsification_nubase", f"{tb['tier3']}", "Tier III 條數")
cl("10_falsification_nubase", f"{100*tb['tier3_frac']:.1f}%", "Tier III 佔比")
cl("10_falsification_nubase", f"{tb['n_non_tier3']}", "Tier II 例外條數")
cl("10_falsification_nubase", f"{100*(1-tb['tier3_frac']):.1f}%", "Tier II 佔比")
cl("10_falsification_nubase", f"{tb['non_tier3_A_range'][0]}", "例外 A 下限")
cl("10_falsification_nubase", f"{tb['non_tier3_A_range'][1]}", "例外 A 上限")
cl("10_falsification_nubase", f"{tb['SP_stbl']}", "S\\P 中 stbl 條數")
cl("10_falsification_nubase", f"{tb['SP_finite']}", "S\\P 中有限 T½ 條數")
cl("10_falsification_nubase", e_(tb["max_rel"], 2), "失效半徑 max rel")
cl("10_falsification_nubase", f"{tb['n_tier1']}", "Tier I 條數")
for nm in ("α", "SF_opt", "2β⁺/ECEC", "2β⁻", "β⁻", "β⁺/EC"):
    cl("10_falsification_nubase", f"{tb['SP_channel_dist'][nm]}", f"S\\P 之 {nm} 條數")
cl("10_falsification_nubase", f"{tb['III_channel_dist']['SF_opt']}", "Tier III 內 SF_opt 條數")
_nt = [t for t in tb["items"] if t["tier"] != "III"]
cl("10_falsification_nubase", e_(min(t["rel"] for t in _nt), 2), "例外 rel 下限")
cl("10_falsification_nubase", e_(max(t["rel"] for t in _nt), 2), "例外 rel 上限")
tc = j10c["test_C"]
cl("10_falsification_nubase", f"{tc['n_ill']}", "不合法集（有限 T½）條數")
cl("10_falsification_nubase", e_(tc["min_yr"], 1), "最短 T½")
cl("10_falsification_nubase", e_(tc["max_yr"], 2), "最長 T½")
cl("10_falsification_nubase", e_(tc["span"], 1), "T½ 跨度")


# ---------------- docs/11 ← 12_definition_d5 (+ 11)
j12 = J["12_definition_d5"]
cl("11_definition_d5", j12["seal"], "D5 封印（逐字）")
cl("11_definition_d5", f"{j10c['P_size']}", "合法集條數")
cl("11_definition_d5", f"{tb['S_size']}", "S(1e15) 條數")
cl("11_definition_d5", f"{tb['SP_size']}", "S\\P 條數")
cl("11_definition_d5", f"{tb['tier3']}", "Tier III 條數")
cl("11_definition_d5", f"{100*tb['tier3_frac']:.1f}%", "Tier III 佔比")
cl("11_definition_d5", f"{tb['n_non_tier3']}", "Tier II 例外條數")
cl("11_definition_d5", f"{tb['non_tier3_A_range'][0]}", "例外 A 下限")
cl("11_definition_d5", f"{tb['non_tier3_A_range'][1]}", "例外 A 上限")
cl("11_definition_d5", f"{tb['max_rel']:.2e}", "失效半徑（doc 格式）")
cl("11_definition_d5", f"{j10c['test_C']['span']:.1e}", "時間盲跨度（doc 格式）")
cl("11_definition_d5", f"{len(j9c['final']['channels'])}", "𝒞_v2 通道類數")
cl("11_definition_d5", f"{j9c['final']['n_legal']}", "𝒞_v2 合法集條數")


# ---------------- docs/12 ← 13_primordial_test
j13 = J["13_primordial_test"]
cl("12_primordial_label", f"{j13['n_nat']}", "天然存在集條數")
cl("12_primordial_label", f"{j13['n_P']}", "P 條數")
cl("12_primordial_label", f"{j13['n_P_and_nat']}", "P∩Nat 條數")
cl("12_primordial_label", f"{j13['n_Nat_not_P']}", "Nat\\P 條數")
cl("12_primordial_label", f"{j13['tier_dist']['III']}", "Nat\\P 之 Tier III")
cl("12_primordial_label", f"{j13['tier_dist']['II']}", "Nat\\P 之 Tier II")
cl("12_primordial_label", f"{j13['max_rel']:.2e}", "失效半徑（丰度標籤）")
cl("12_primordial_label", f"{j13['n_mono']}", "單一同位素元素個數")
for r in j13["mono"]:
    cl("12_primordial_label", f"{r['rel']:.2e}", f"單一同位素 Z={r['Z']} 之 rel")


# ---------------- docs/13 ← 14_barrier_gate
j14 = J["14_barrier_gate"]
cl("13_barrier_gate", f"{j14['ldm']['a_s']:.2f}", "a_s")
cl("13_barrier_gate", f"{j14['ldm']['kappa']:.4f}", "κ")
cl("13_barrier_gate", f"{j14['ldm']['coeff_Bf']:.2f}", "B_f 係數")
cl("13_barrier_gate", f"{j14['n_sf_observed']}", "觀測 SF 條數")
cl("13_barrier_gate", f"{j14['n_band']}", "屏障帶條數")
cl("13_barrier_gate", f"{j14['n_P']}", "𝒞_v2 合法集")
cl("13_barrier_gate", f"{j14['G1']['n_violations']} 條違反", "G1 違反數")
_v = max(j14["G1"]["violations"], key=lambda r: r["Q"] or 0)
cl("13_barrier_gate", f"{_v['Q']:.1f}", "最大 Q_SF")
cl("13_barrier_gate", f"{_v['Bf']:.1f}", "該核素 B_f")
cl("13_barrier_gate", f"{_v['ratio']:.2f}", "該核素 Q/B_f")
g2 = j14["G2"]
cl("13_barrier_gate", f"{g2['bf_sf_min']:.1f}", "SF 集 B_f min")
cl("13_barrier_gate", f"{g2['bf_sf_med']:.1f}", "SF 集 B_f 中位")
cl("13_barrier_gate", f"{g2['bf_sf_max']:.1f}", "SF 集 B_f max")
cl("13_barrier_gate", f"{g2['bf_band_min']:.1f}", "屏障帶 B_f min")
cl("13_barrier_gate", f"{g2['bf_band_med']:.1f}", "屏障帶 B_f 中位")
cl("13_barrier_gate", f"{g2['bf_band_max']:.1f}", "屏障帶 B_f max")
cl("13_barrier_gate", f"{g2['gap']:.1f}", "兩集 B_f 間隙")
cl("13_barrier_gate", f"{g2['B_star']:.2f}", "B* 門檻")
g3 = j14["G3"]
cl("13_barrier_gate", f"{g3['P_gated']}", "P_gated 條數")
cl("13_barrier_gate", f"{g3['n_nonpermanent']}", "P_gated 中非永久者")
cl("13_barrier_gate", f"{100*g3['frac']:.1f}%", "非永久者佔比")
cl("13_barrier_gate", f"{g3['SP_gated']}", "S\\P_gated 條數")
cl("13_barrier_gate", f"{100*g3['tier3_frac_gated']:.1f}%", "gated 之 Tier III 佔比")

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
seals = {j1["C5_seal"]["sha256"], j5["L_seal"], J["12_definition_d5"]["seal"]}
found = set(re.findall(r"\b[0-9a-f]{64}\b", " ".join(DOC.values())))
ck("C3.1 doc 中 64-hex 皆為 data 內封印", found <= seals,
   f"doc {len(found)} 個 / data {len(seals)} 個；多出 {sorted(found-seals)}")
ck("C3.2 封印互異（判據 ≠ 預言；D5 亦互異）",
   len({j1["C5_seal"]["sha256"], j5["L_seal"], J["12_definition_d5"]["seal"]}) == 3,
   f"判據 {j1['C5_seal']['sha256'][:12]}… / 預言 {j5['L_seal'][:12]}…")

print("\n[C4] 覆蓋度：")
per = {}
for doc, s, why in CLAIMS:
    per[doc] = per.get(doc, 0) + 1
for doc in sorted(per):
    print(f"    {doc:<34} {per[doc]:>3} 條")
ck("C4 每份 doc 至少 5 條", all(v >= 5 for v in per.values()) and len(per) == 13,
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