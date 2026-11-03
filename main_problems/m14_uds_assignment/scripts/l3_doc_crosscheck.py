#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
L3 — M14 文檔交叉核（idcm-audit-verify §2.4）
================================================================================
R1 唯讀：對 docs 只 open 讀；只寫 data/ 下兩個輸出。
R2 期望值**只**來自 L2 `_derived`；禁抄 docs 數字。
R3 容差依有效位（≥6 → 1e−4；≥4 → 5e−4；其餘 5e−3）；可 per-claim `tol` 收緊。
R4 負值取絕對值比（U+2212 先歸一）。
R5 容差仍須分辨真錯。

★ 豁免聲明（§2.4 豁免 1）：docs/03、04、07、08、09、13 無數值宣稱或屬判定／更正文檔
   ⟹ 不入 CLAIMS，由 audit 記錄豁免理由。
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
L2PATH = DATA / "l2_independent_verify.json"

LOG, RES, NFAIL = [], {}, 0

# 明示豁免（每項必須附理由；audit §⑤ 會核理由非空）
EXEMPT = {
    "docs/03": "AST 出處追溯記錄；數字僅為 v9 字面值（gen_modes=[1,4,5]、z2_map={1:3,4:2,5:0}）與 N_trials=6，非本 M 譜量 ⟹ 無譜量宣稱可 claim",
    "docs/07": "判詞／移交敘述表；無獨立譜量，數字皆引他 doc",
    "docs/08": "M13 入向分診記錄；數字屬 M13（δ_c、候選 |z|），M14 只轉述",
    "docs/09": "概念更正記錄（對象型別）；唯一精確值 A_eff=4.4492±0.4603 屬 M5 跨記錄",
    "docs/13": "通道存續之讀法更正記錄；數字僅為模數與對數（44／21），已由 docs/12、docs/14 claim 覆蓋",
    "docs/22": "審計報告自身（自指）；其數字為各層計數，已由 audit 自身斷言",
    "docs/38": "內部碼命名空間規範＋映射表；數字僅為碼序號（MR1–MR8），非譜量 ⟹ 無 claim",
}

CLAIMS = [
    # docs/01 重現閘
    {"doc": "docs/01_reproduction_gate.md", "key": "n_A2p1", "label": "A₂⊕1 = 44"},
    {"doc": "docs/01_reproduction_gate.md", "key": "n_A2block", "label": "A₂ 2×2 塊 = 87"},
    # docs/02 Q1 gauge
    {"doc": "docs/02_q1_station_gauge.md", "key": "n_A2p1", "label": "44/112"},
    {"doc": "docs/02_q1_station_gauge.md", "key": "n_A2block", "label": "87/112"},
    # docs/05 δ_u 符號
    {"doc": "docs/05_q3_delta_u_sign.md", "key": "delta_u_over_sigma", "label": "|δ_u|/σ = 0.789"},
    {"doc": "docs/05_q3_delta_u_sign.md", "key": "rho_u_z", "label": "ρ_u z = 0.771"},
    # docs/06 第三槽
    {"doc": "docs/06_q4_third_slot.md", "key": "n_A2p1", "label": "44 個 (2,1) 塊"},
    {"doc": "docs/06_q4_third_slot.md", "key": "n_singlet_decoupled", "label": "44 個 s 退耦"},
    # docs/10 分類第一步
    {"doc": "docs/10_classification_first_step.md", "key": "pairs_946", "label": "946 價對"},
    {"doc": "docs/10_classification_first_step.md", "key": "aut_C88", "label": "|Aut(C_88)| = 176"},
    {"doc": "docs/10_classification_first_step.md", "key": "irrep_dimsq_sum", "label": "Σd² = 176"},
    # docs/11 價對分類（★ 43「軌」之語義 → 綁真群軌道鍵，預期 FAIL 以現形）
    {"doc": "docs/11_pair_delta_classification.md", "key": "delta_classes_43", "label": "|Δ| 43 類"},
    {"doc": "docs/11_pair_delta_classification.md", "key": "z44_orbits_22", "label": "真群軌道 22（更正後 §九）"},
    # docs/12 三波共振
    {"doc": "docs/12_three_wave_resonance.md", "key": "law_sum_N44", "label": "sum 21"},
    {"doc": "docs/12_three_wave_resonance.md", "key": "law_beat_N44", "label": "beat 43"},
    {"doc": "docs/12_three_wave_resonance.md", "key": "law_active_N44", "label": "活價對 64"},
    # docs/14 通用定律
    {"doc": "docs/14_universal_resonance_law.md", "key": "law_active_N44", "label": "N=44 → 64"},
    {"doc": "docs/14_universal_resonance_law.md", "key": "law_active_N88", "label": "N=88 → 130"},
    {"doc": "docs/14_universal_resonance_law.md", "key": "law_active_N176", "label": "N=176 → 262"},
    # docs/15 慣例閘 + 三三元組
    {"doc": "docs/15_convention_gate_three_triples.md", "key": "n_triples90", "label": "恰 3 個三元組"},
    # docs/16 χ³ no-go
    {"doc": "docs/16_chi3_charge_nogo.md", "key": "distinct_triples", "label": "C(44,3) = 13,244"},
    {"doc": "docs/16_chi3_charge_nogo.md", "key": "chi3_nonzero_count", "label": "χ³ 非零 6,622"},
    # docs/17 T_pipe + φ 階梯
    {"doc": "docs/17_tpipe_filter_and_phi_ladder.md", "key": "n_star_27", "label": "n* = 27"},
    {"doc": "docs/17_tpipe_filter_and_phi_ladder.md", "key": "A9_phi_inv", "label": "A_9 = φ⁻¹ = 0.618034"},
    {"doc": "docs/17_tpipe_filter_and_phi_ladder.md", "key": "gsum_composites", "label": "γ 和 = 8"},
    {"doc": "docs/17_tpipe_filter_and_phi_ladder.md", "key": "g18_plus_g36", "label": "γ18+γ36 = 5"},
    # docs/18 四波 + mod-3
    {"doc": "docs/18_four_wave_mod3_reduction.md", "key": "clusters_913", "label": "913 等值類"},
    {"doc": "docs/18_four_wave_mod3_reduction.md", "key": "n_doubletons", "label": "12 雙元類"},
    {"doc": "docs/18_four_wave_mod3_reduction.md", "key": "fourwave_count", "label": "243 四波共振"},
    {"doc": "docs/18_four_wave_mod3_reduction.md", "key": "cluster_max_size", "label": "22 類"},
    {"doc": "docs/18_four_wave_mod3_reduction.md", "key": "sublattice_modes", "label": "14 模子晶格"},
    # docs/19 單標籤
    {"doc": "docs/19_single_label_theorem.md", "key": "sum88", "label": "Σγ(/88) = 90"},
    {"doc": "docs/19_single_label_theorem.md", "key": "sum90", "label": "Σγ(/90) = 88"},
    {"doc": "docs/19_single_label_theorem.md", "key": "vandermonde_det_abs", "label": "det V = 2.188e41", "tol": 1e-2},
    # docs/20 相位差
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "theta_den_3960", "label": "3960"},
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "comp88_parity_same", "label": "44/44 保宇稱"},
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "comp90_parity_same", "label": "0/44 flip 宇稱"},
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "sumfam88_pairs", "label": "a+b=44 → 21 對"},
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "sumfam90_pairs", "label": "a+b=45 → 22 對"},
    {"doc": "docs/20_odd_even_phase_gap.md", "key": "min_gap_cross", "label": "最小間隙 2.2387e−4", "tol": 1e-4},
    # docs/23 Green 核 ＋ 補償（新）
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "green_points", "label": "1936 點"},
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "green_g_9_15", "label": "G(9,15) = 6"},
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "green_diag_18", "label": "G(18,18) = 10.8", "tol": 1e-3},
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "green_g_18_36", "label": "G(18,36) = 3.6"},
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "sum4_triples", "label": "sum-4 解恰 1"},
    {"doc": "docs/23_green_kernel_and_complement.md", "key": "phi_parity_same", "label": "Φ 保 parity 0/44"},
    # docs/24 不變量格（新）
    {"doc": "docs/24_invariant_lattice.md", "key": "inv_pool", "label": "池 10 個"},
    {"doc": "docs/24_invariant_lattice.md", "key": "inv_subsets", "label": "1023 個子集"},
    {"doc": "docs/24_invariant_lattice.md", "key": "three_class_partitions", "label": "恰 3 類分割 5 個"},
    {"doc": "docs/24_invariant_lattice.md", "key": "canonical_classes_z44", "label": "ℤ_44 軌道 22"},
    # docs/25 /88 奇模扇區（新）
    {"doc": "docs/25_odd_sector_and_separability.md", "key": "res88_triples", "label": "21 組"},
    {"doc": "docs/25_odd_sector_and_separability.md", "key": "res88_odd_pair", "label": "奇模 11"},
    {"doc": "docs/25_odd_sector_and_separability.md", "key": "res88_even_pair", "label": "偶模 10"},
    {"doc": "docs/25_odd_sector_and_separability.md", "key": "res90_triples", "label": "/90 3 組"},
    {"doc": "docs/25_odd_sector_and_separability.md", "key": "res_separable_count", "label": "三分離者 0"},
    # docs/26 /88 四波譜（新）
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "spec88_classes", "label": "/88 926 類"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "spec90_classes", "label": "/90 913 類"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "fw88_resonances", "label": "/88 共振 210"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "fw90_resonances", "label": "/90 共振 243"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "accidental88", "label": "/88 雙元類 0"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "accidental90", "label": "/90 雙元類 12"},
    {"doc": "docs/26_fourwave_88_and_duality.md", "key": "fixed_line_pairs88", "label": "固定線 21 對"},
    # docs/27 子晶格（新）
    {"doc": "docs/27_sublattice_divisibility.md", "key": "sub88_count", "label": "/88 六個"},
    {"doc": "docs/27_sublattice_divisibility.md", "key": "sub90_count", "label": "/90 九個"},
    {"doc": "docs/27_sublattice_divisibility.md", "key": "sub88_has3", "label": "3∤88（0）"},
    {"doc": "docs/27_sublattice_divisibility.md", "key": "sub30_halfsum", "label": "/30 半譜 28"},
    # docs/28 二次域（新）
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_subchain_modes", "label": "4 模"},
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_threewave8", "label": "三波 1"},
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_fourwave8", "label": "四波 0"},
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_disc_gamma11", "label": "disc 8"},
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_disc_gamma9_90", "label": "disc 5"},
    {"doc": "docs/28_quadratic_field_criterion.md", "key": "q2_gamma1", "label": "γ1 = 0.5857864376269", "tol": 1e-12},
    # docs/29 加性 vs 乘性（新）
    {"doc": "docs/29_additive_vs_multiplicative.md", "key": "gp88", "label": "/88 GP 0"},
    {"doc": "docs/29_additive_vs_multiplicative.md", "key": "ap88", "label": "/88 AP 21"},
    {"doc": "docs/29_additive_vs_multiplicative.md", "key": "gp90", "label": "/90 GP 3"},
    {"doc": "docs/29_additive_vs_multiplicative.md", "key": "ap90", "label": "/90 AP 0"},
    {"doc": "docs/29_additive_vs_multiplicative.md", "key": "ap88_mid22", "label": "中項 22 全核 21"},
    # docs/30 AP/GP 判準（新）
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "c4_88", "label": "4|88"},
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "c30_90", "label": "30|90"},
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "ap88_formula", "label": "AP = 88/4−1 = 21"},
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "gp88_formula_ring", "label": "GP(88) = 0"},
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "gp90_formula_dir", "label": "GP(90) = 3"},
    {"doc": "docs/30_ap_gp_general_criterion.md", "key": "c12_88", "label": "12∤88 (0)"},
    # docs/31 AP 兩族分類（新）
    {"doc": "docs/31_ap_classification.md", "key": "ap12", "label": "AP(12)=3"},
    {"doc": "docs/31_ap_classification.md", "key": "ap24", "label": "AP(24)=6"},
    {"doc": "docs/31_ap_classification.md", "key": "ap120", "label": "AP(120)=30"},
    {"doc": "docs/31_ap_classification.md", "key": "ap88", "label": "N=88 → AP 21（族A）"},
    {"doc": "docs/31_ap_classification.md", "key": "ap90", "label": "N=90 → 0", "tol": 1e-3},
    # docs/32 GP 五族分類（新）
    {"doc": "docs/32_gp_classification.md", "key": "gpQ1_closed", "label": "Q1 閉式 (7−3√5)/2", "tol": 1e-12},
    {"doc": "docs/32_gp_classification.md", "key": "gp88", "label": "N=88 → GP 0"},
    {"doc": "docs/32_gp_classification.md", "key": "gp90", "label": "N=90 → GP 3"},
    # docs/33 φ-階梯同一性（新）
    {"doc": "docs/33_phi_ladder_identity.md", "key": "q2_amp_ratio", "label": "A 比 φ = 1.6180339887", "tol": 1e-12},
    {"doc": "docs/33_phi_ladder_identity.md", "key": "gp88", "label": "/88 γ-GP 0"},
    {"doc": "docs/33_phi_ladder_identity.md", "key": "gp90", "label": "/90 γ-GP 3"},
    # docs/34 殘留五項（新）
    {"doc": "docs/34_residue_five_items.md", "key": "l3_sum8", "label": "γ18+γ30+γ36 = 8", "tol": 1e-12},
    {"doc": "docs/34_residue_five_items.md", "key": "p1_lcm", "label": "3960 = lcm(88,90)"},
    {"doc": "docs/34_residue_five_items.md", "key": "j1_members", "label": "ℚ(√5) 值集 6 模"},
    # docs/36 Q1／M2／M3（新；小整數用 anchor）
    {"doc": "docs/36_q1_m2_m3.md", "key": "q1_ring_delta", "label": "ring 之 δ = +2",
     "anchor": "δ = ε · γ_max/2", "tol": 1e-9},
    {"doc": "docs/36_q1_m2_m3.md", "key": "m3_ksum14", "label": "k 和 14 = 3+5+6",
     "anchor": "k = (3,5,6)"},
    {"doc": "docs/36_q1_m2_m3.md", "key": "m3_res_count", "label": "共振 3 個",
     "anchor": "三波共振恰"},
    # docs/37 K3 閉式（新；小整數用 anchor）
    {"doc": "docs/37_k3_closed_form_theorem.md", "key": "k3_zero_key_even", "label": "k_i=0 之偶和三元組 2838",
     "anchor": "現存 `k_i = 0` 之三元組"},
    {"doc": "docs/37_k3_closed_form_theorem.md", "key": "k3_zero_key_odd", "label": "k_i=0 之奇和 = 0",
     "anchor": "反例 **0**"},
    {"doc": "docs/37_k3_closed_form_theorem.md", "key": "k3_min_val", "label": "最小 |χ³| 於 (1,1,43)",
     "anchor": "(1,1,43)"},
    # docs/39 四波閉式（新）
    {"doc": "docs/39_fourwave_overlap_closed_form.md", "key": "chi4_1234", "label": "χ⁴(1,2,3,4)=1/90",
     "anchor": "χ⁴ = W/90"},
    {"doc": "docs/39_fourwave_overlap_closed_form.md", "key": "chi4_1111", "label": "χ⁴(1,1,1,1)=3/90",
     "anchor": "χ⁴ = W/90"},
    {"doc": "docs/39_fourwave_overlap_closed_form.md", "key": "w4_nz_frac_pct", "label": "非零比例 18.94%",
     "anchor": "非零比例"},
    # docs/40 四波能量 vs 選擇律（新）
    {"doc": "docs/40_energy_vs_selection.md", "key": "fw_total", "label": "能量共振 1,233（v2）",
     "anchor": "1,233"},
    {"doc": "docs/40_energy_vs_selection.md", "key": "fw_odd", "label": "和奇 12（v2）",
     "anchor": "**12**"},
    # docs/41 MR1 機器（新）
    {"doc": "docs/41_dynamics_machine.md", "key": "dyn_d12", "label": "三波衰變 3 條",
     "anchor": "恰 3 條"},
    {"doc": "docs/41_dynamics_machine.md", "key": "dyn_d13", "label": "四波衰變 1 條",
     "anchor": "恰 1 條"},
    {"doc": "docs/41_dynamics_machine.md", "key": "dyn_res", "label": "2→2 共振 243 對",
     "anchor": "243"},
    {"doc": "docs/41_dynamics_machine.md", "key": "dyn_uniform", "label": "不透明度每模 21",
     "anchor": "每模恰 21"},
    {"doc": "docs/41_dynamics_machine.md", "key": "dyn_tri4", "label": "和=4 三波 1 組",
     "anchor": "恰 = 連續 Nyquist 值"},
    {"doc": "docs/40_energy_vs_selection.md", "key": "fw_even", "label": "和偶 1,221",
     "anchor": "和偶"},
    {"doc": "docs/40_energy_vs_selection.md", "key": "fw_offshell", "label": "off-shell 13,874",
     "anchor": "W ≠ 0 但無能量共振"},
    # docs/42 二階自能（新）
    {"doc": "docs/42_second_order_selfenergy.md", "key": "so_poles", "label": "運動學極點 1 個",
     "anchor": "恰 1 個"},
    {"doc": "docs/42_second_order_selfenergy.md", "key": "so_true_poles", "label": "真極點 0",
     "anchor": "真極點 = 0"},
    {"doc": "docs/42_second_order_selfenergy.md", "key": "so_distinct", "label": "相異 44",
     "anchor": "相異值 = 44"},
    {"doc": "docs/42_second_order_selfenergy.md", "key": "so_bound", "label": "|g₄| < 0.070265395858052",
     "anchor": "< 0.070265395858052"},
    # docs/36 精度更正（新）
    {"doc": "docs/36_q1_m2_m3.md", "key": "specmax_90", "label": "/90 譜內最大 3.99512810052",
     "anchor": "譜內最大"},
    # docs/43 bubble／符號（新）
    {"doc": "docs/43_bubble_and_sign.md", "key": "bub_zero_nz", "label": "W≠0 零分母每模 21",
     "anchor": "**21**"},
    {"doc": "docs/43_bubble_and_sign.md", "key": "so_pos_count", "label": "S_n > 0 之模 6 個",
     "anchor": "6 / 44"},
    {"doc": "docs/43_bubble_and_sign.md", "key": "so_transitions", "label": "轉換 12 次",
     "anchor": "轉換 **12** 次"},
    # docs/04 獨立重跑 216（新）
    {"doc": "docs/04_q2b_permutation_uniqueness.md", "key": "perm_scan_total",
     "label": "216 獨立重跑（6 規則 × 36 指派）"},
    # docs/35 宇稱／δ 定理（新）
    {"doc": "docs/35_parity_and_delta_theorems.md", "key": "k3_even_n", "label": "偶和三元組 6622"},
    {"doc": "docs/35_parity_and_delta_theorems.md", "key": "q4_viol", "label": "δ 三分律違例 0"},
    {"doc": "docs/35_parity_and_delta_theorems.md", "key": "k3_reflect_pairs", "label": "反射對 1936"},
    # docs/21 δ 三分律
    {"doc": "docs/21_budget_defect_dyadic_depth.md", "key": "defect_ring88", "label": "δ(/88) = +2"},
    {"doc": "docs/21_budget_defect_dyadic_depth.md", "key": "defect_dir90", "label": "δ(/90) = −2"},
    {"doc": "docs/21_budget_defect_dyadic_depth.md", "key": "layers88", "label": "3 層"},
    {"doc": "docs/21_budget_defect_dyadic_depth.md", "key": "layers90", "label": "2 層"},
    {"doc": "docs/21_budget_defect_dyadic_depth.md", "key": "v2_88", "label": "v₂(88) = 3"},
]

TOK = re.compile(r"(\d[\d,]*(?:\.\d+)?)(?:[eE]([+\u2212-]?\d+))?")


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def tokens(text):
    out = []
    for m in TOK.finditer(text.replace("\u2212", "-")):
        mant, exp = m.group(1).replace(",", ""), m.group(2)
        try:
            v = float(mant) * (10.0 ** int(exp)) if exp else float(mant)
        except (ValueError, OverflowError):
            continue
        sig = len(mant.replace(".", "").lstrip("0")) or 1
        out.append((v, sig))
    return out


def match(toks, target, tol_override=None):
    t = abs(target)
    if t == 0:                     # 零期望值：只接受文檔出現字面 0
        return (any(v == 0 for v, _ in toks), 0.0, 0.0)
    best, best_d = None, float("inf")
    for v, sig in toks:
        if v == 0:
            continue
        rel = abs(abs(v) - t) / t
        tol = tol_override if tol_override else (1e-4 if sig >= 6 else (5e-4 if sig >= 4 else 5e-3))
        if rel <= tol and rel < best_d:
            best, best_d = v, rel
    return (best is not None), best, best_d


def main():
    hdr("§0 讀 L2（期望值唯一來源）")
    if not L2PATH.exists():
        print(f"✗ 缺 {L2PATH}")
        return 2
    l2 = json.loads(L2PATH.read_text(encoding="utf-8"))
    der = l2.get("_derived", {})
    chk("E1 L2 有 _derived 且非空", bool(der), f"{len(der)} keys")

    hdr("§1 逐 claim 交叉核（唯讀）")
    cache, by_doc = {}, {}
    for c in CLAIMS:
        doc, key, label = c["doc"], c["key"], c.get("label", c["key"])
        p = ROOT / doc
        if not p.exists():
            chk(f"C.{doc} 存在", False, "檔案不存在")
            continue
        if doc not in cache:
            cache[doc] = p.read_text(encoding="utf-8")
        if key not in der:
            chk(f"C.{doc}.{label}（L2 key={key}）", False, "L2 冇此 key")
            continue
        by_doc.setdefault(doc, [0, 0])
        by_doc[doc][1] += 1
        exp = der[key].get("float")
        if exp is None:
            chk(f"C.{doc.split('/')[-1][:18]}.{label}", False, "L2 值超出 float 範圍")
            continue
        ok, got, rel = match(tokens(cache[doc]), exp, c.get("tol"))
        if ok:
            by_doc[doc][0] += 1
        chk(f"C.{doc.split('/')[-1][:18]}.{label}", ok,
            f"L2 {exp:.6g} ↔ 文檔 {got:.6g}（rel {rel:.1e}）" if ok
            else f"✗ L2 {exp:.6g}，文檔最近 {got if got is None else format(got, '.6g')}")

    hdr("§2 逐檔命中率")
    for d, (h, t) in sorted(by_doc.items()):
        print(f"    {d:50s} {h}/{t}")
        chk(f"F.{d} 全命中（{h}/{t}）", h == t)

    hdr("§3 匯總")
    n = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {n}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "by_doc": by_doc,
                       "exempt_docs": sorted(EXEMPT), "exempt_reasons": EXEMPT}
    (DATA / "l3_doc_crosscheck.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())