# docs/56 — M14 腳本與過程記錄清冊（scripts manifest）


> 建立：2026-11-25｜用途：滿足「全推導鏈腳本皆存於 scripts／過程與結果皆入文檔」之檢查
> 判準：任何在推導鏈中被執行並產生文檔數字之腳本，必須（i）存於 `scripts/`（含 `scripts/oneoff/`）且（ii）在文檔指明其產物。


## §一 L1 正典腳本（50 支，02–48、51–53 與審計器）

```
支數 = 51（含 23_audit_verify 審計器；L1 集 = 其餘 50）
01_reproduction_gate.py 02_q1_station_gauge.py 03_q2a_label_provenance.py 04_q2b_permutation_uniqueness.py 05_q3_delta_u_sign.py 06_q4_third_slot.py 07_m13_incoming_reply.py 08_classification_partitions.py 09_pair_galois_orbits.py 10_symmetry_group_irreps.py 11_pair_physical_classification.py 12_delta_support.py 13_channel_fate_classification.py 14_three_wave_resonance.py 15_general_resonance_law.py 16_convention_and_three_triples.py 17_chi3_charge_coefficient.py 18_triples_vs_tpipe.py 19_four_wave_and_mod3.py 20_single_label_theorem.py 21_odd_even_phase_gap.py 22_budget_defect_and_dyadic_depth.py 23_audit_verify.py 24_charge_kernel_and_complement.py 25_invariant_lattice.py 26_odd_sector_and_partition_separability.py 27_fourwave_88_and_duality.py 28_sublattice_divisibility.py 29_qsqrt2_sector_and_field_criterion.py 30_ladder_ap_vs_gp.py 31_general_ap_gp_criterion.py 32_ap_classification.py 33_gp_classification.py 34_phi_ladder_identity.py 35_residue_five_items.py 36_parity_and_delta_theorems.py 37_reverse_coverage.py 38_q1_m2_m3.py 39_k3_closed_form_proof.py 40_fourwave_overlap_theorem.py 41_energy_vs_selection.py 42_cross_section_machine.py 43_second_order_selfenergy.py 44_bubble_and_sign.py 45_third_order.py 46_fourth_order.py 47_quint_sector.py 48_fifth_order.py 51_precise_vertex.py 52_fourwave_operator.py 53_vertex_gpu.py
```

## §二 分層驗證

```
scripts/l2_independent_verify.py → data/run_l2_independent_verify.log
scripts/l3_doc_crosscheck.py   → data/run_l3_doc_crosscheck.log
```

## §三 oneoff 工具（8 支）

```
e7chain.py｜e7run.py｜e7runB.py｜gpu_op.py｜gpu_run.py｜matrix.py（E7／E8 鏈與 GPU 對照）
cs_sig_gate.py（CS = ∏ε 之閘）｜m0_resid_gate.py（M₀ 殘差族覆蓋＋null）
```

## §四 探針（scripts/oneoff/probes_m14／27 支）

```
【MR-A】probe_mrA_adjoint / adjoint2 / decode / decode2 / decode3 / final / rankcheck / raw (8)
   產物：data/mrA_gate.json；data/mrA_*．log
   文檔：docs/49 §五
   角色：多重度公約 ∏C vs raw 之伴隨檢定
【MR-B】probe_mrB_decomp / 48 / split / terms / law_L3 (5)
   產物：data/mrB_decomp.json、mrB_terms6.json、mrB_terms44.json、mrB_law44.json；data/mrB_lawL3_attempt.log
   文檔：docs/49 §六／§七／§八／§九／§十
   角色：rel 恆等式、四項分解、44 n 掃描、三層還原
【± 結構／相消】probe_vsign35 / vinter / sign_mech / sign_bound / a1_cancel / discord / mass (7)
   產物：data/mrB_vsign35.json、mrB_vinter35.json、mrB_sign_mech.json、mrB_sign_bound.json、mrB_a1_cancel.json、mrB_discord.json、mrB_mass35.json
   文檔：docs/49 §十一～§十六
   角色：w ∈ ℤ、路徑干涉、A2 = ‖w‖²、三路否證
【M₀ 標度】scripts/oneoff/m0_resid_gate.py（＋cs_sig_gate.py）
   產物：data/m0_resid_gate.json、m0_resid_gate.log；data/mrB_cs_sig.json
   文檔：docs/55；docs/49 §十二
   角色：殘差族覆蓋／null；CS = ∏ε
【E_7／E_8 鏈】e8_79p / e8_97b / e8_97c / e8b_97 / e8_tail / e8_tailB_C / e8_tailC (7)；oneoff/e7chain、e7run、e7runB、gpu_op、gpu_run、matrix
   產物：data/chain_matrix.json、data/e8_chain_clean.json、data/chain_logs/*；runtime p97C/p97D（ckpt）
   文檔：docs/53（矩陣）、docs/22 §二十七
   角色：四維矩陣 E1–E8；A/B 兩路徑
```

## §五 缺口（誠實聲明）

```
① 重計算資料（p97C/p97D、mrB/decomp48_all*.json、terms44 原表）為 runtime 產物，存於 /home/wsl/idcm_ckpt/e8，不入 git（*.npy 由 .gitignore 排除）
② scripts/oneoff/probes_m14 為執行時副本（原檔在 idcm_ckpt），兩者內容逐位相同（copy -p）
③ 所有文檔引用數字皆有 data/*.json|*.log 落地（37_reverse_coverage 未追溯 = 0）
```

## §六 更新（2026-11-26）：探針 27 → 45 支（新增 18）

```
新增探針與文檔映射（全部由本輪 §十七–§二十六 產生）：
  probe_rank1.py     → docs/49 §十七（秩一：sign(t2_r1)=sign(t2) 44/44；主導模份額）
  probe_tail5.py     → docs/49 §十七（N_eff(w²)、top-10 份額、近退化態數）
  probe_r1reason.py  → docs/49 §十八（H1 最近模 5/44 否；gap2/gap1；N_eff(1/den²)）
  probe_r1support.py → docs/49 §十八（H2 允許集最近模 14/44 否；支撐集 441,737／1,712,304）
  probe_degen5.py    → docs/49 §十九（EQ 簡併普查 815,171 distinct、最大重數 360）
  probe_offdelta.py  → docs/49 §十九（44 n 之 δ 分組：22 個精確值，dps 60）
  probe_delta_alg.py → docs/49 §二十（22 δ 之 PSLQ 有理判定：0/22）
  probe_gammadim.py  → docs/49 §二十（ℚ-維數 12 判定；360 重數簇 = 五餘弦零和）
  probe_oos_dim.py   → docs/49 §二十一（OOS-1：m88 dim40 預言 ≥40 ⟹ 實測 28 FAIL）
  probe_oos_dim2.py  → docs/49 §二十一（OOS-2：m45 20／m180 42 PASS）
  probe_oos_dim3.py  → docs/49 §二十二（OOS-3：R1／R2 關係數 ⟹ m60 18／m176 44 PASS）
  probe_wclass.py    → docs/49 §二十三（四腿 W 全域 667,920 簽名：5 值計數）
  probe_w3struct.py  → docs/49 §二十四（w=+3 之 88 簽名 = (k,k,k,k) ⊕ (k,k,45−k,45−k)）
  probe_wmask.py     → docs/49 §二十五（16 遮罩 census；popcount 律）
  probe_wclassify.py → docs/49 §二十五（重數型／反射對數分佈，輔證）
  probe_w17_check.py → docs/49 §二十六（12,000 路徑 w_rec ≡ W_table 零差 ⟹ §十一 值域更正）
庫／驅動（load-bearing，非探針）：
  e7chain.py（E₇ 鏈模組，被多支探針 import）｜e7run.py（E₇ 驅動，表 QUI／TRI 之來源）
⟹ 修正 §四 之計數：27 → 45（本節為 append-only 更新，§四 原文不動）
```

## §七 更新（2026-11-27 收尾）：探針實測 67 支＋全量索引（append-only；§四／§六 原文不動）

> 計數更正：§四 記 27 支、§六 記 45 支；本節實測 **67 支**（含本次補歸檔 1 支）。
> 原文保留為歷史；本節為生效值。

### 7.1 統計

```
probe 總數                67
有 json 產物              51（產物落 data/，命名前綴規則見 7.2）
無 json（只落 docs／stdout）16（見 7.4）
清冊 §四／§六 曾列名       28
```

### 7.2 命名前綴規則（避免再誤判「缺產物」）

```
probe 內寫 mrB/x.json      ⟹ data/ 落 mrB_x.json
probe 內寫 op_*.json       ⟹ data/ 落 op_*.json（無前綴）
probe 內寫 phys_*.json     ⟹ data/ 落 mrB_phys_*.json
probe 內寫 s2_/s3_/s4_     ⟹ data/ 落 mrB_s*_*.json
probe 內寫 r1_/r4_/r5_     ⟹ data/ 落 mrB_r*_*.json
```

### 7.3 新增探針索引（本輪 §三十五–§四十九 所用；產物 → 文檔節）

```
probe_r1_r5                      mrB_r1_masks_ordered／mrB_r5_distribution／mrB_dominance   → §三十五／§三十六
probe_r1_mask_algebra            mrB_r1_mask_algebra                                      → §三十九
probe_r4_support45               mrB_r4_support45／mrB_offdelta                           → §三十八
probe_phys_mass_ceiling          mrB_phys_mass_ceiling                                    → §四十
probe_phys_power_fit2            mrB_phys_power_fit                                       → §四十
probe_phys_power_null            mrB_phys_power_null                                      → §四十
probe_phys_lattice_membership    mrB_phys_lattice_membership                              → §四十一
probe_s2_lambda_trap             mrB_s2_lambda_trap                                       → §四十二
probe_s3_pair_orbits             mrB_s3_pair_orbits                                       → §四十三
probe_s3_assignment_dof          mrB_s3_assignment_dof                                    → §四十三
probe_s4_chunk_validate          mrB_s4_chunk_validate                                    → §四十四
probe_op_ring_identities         op_ring_identities                                       → §四十五
probe_op_sign_homogeneity        op_sign_homogeneity_proof                                → §四十五
probe_op_w_fourier               op_w_fourier                                             → §四十六
probe_op_fourwave_H              op_fourwave_H                                            → §四十七
probe_op_z_normalization         op_z_normalization                                       → §四十八
probe_op_dual_chain_nogo         op_dual_chain_nogo                                       → §四十九
```

### 7.4 無 json 產物之 16 支（誠實聲明）

```
庫／驅動（5）：e7chain（被多支 import）、e8_79p、e8_97b、e8_97c、e8b_97
舊探針（11）：probe_cs_vs_ident、probe_mrA_adjoint、probe_mrA_adjoint2、probe_mrA_decode、
              probe_mrA_decode2、probe_mrA_decode3、probe_mrA_final、probe_mrA_rankcheck、
              probe_mrA_raw、probe_w3struct、probe_wclassify
處置：其數值只落 docs（§十七–§三十四 各節）；重跑可重建但無 json 對賬 ⟹ 標【doc-only】
      （不補寫 json：事後補寫屬回溯污染，禁用）
```

### 7.5 歸檔動作（本節）

```
cp -p /home/wsl/idcm_ckpt/e8/probe_phys_power_fit.py scripts/oneoff/probes_m14/
  理由：ckpt 存在而 scripts/ 缺（①a 差集唯一項）；已被 probe_phys_power_fit2 取代 ⟹ 標【superseded】
  產物 mrB_phys_power_fit.json 由 fit2 產生（fit 為中間版本，保留作歷史）
複核：ckpt → scripts/ 差集 = 0
```
