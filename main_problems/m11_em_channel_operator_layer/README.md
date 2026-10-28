# M11 — em 通道算子層（U5 攻堅）：通道級電荷讀數算子之結構存在性

**pid** `M11` ｜ **created** 2026-10-28 ｜ **state** 🟡 partial ｜ **version** v11.0
**上游**：M4 `docs/19` U5（「em 通道之正確攻堅點 = 矩陣元；真正上游係一個未建立之結構層」）；
M4 `docs/37 §一 B4` 在版本紀律（2026-10-28 禁止反向推送）下懸空 ⟹ 本記錄承接。

## 問題（一句）
是否存在**零新參數**之結構算子 `Ô`，使 `E_ch = ⟨Ô⟩` 重現 7 通道 IB 質量差譜（跨度 5.66 MeV）？

## 首輪結果（`docs/01`，`scripts/01`，12 checks 10 PASS）
| 項 | 結果 |
|:--|:--|
| `E_ch` 重算 | ✅ 8/8 與 M4 `docs/18` 一致（max\|Δ\| = 3e-5 MeV） |
| ★ 定義 | 🔴 `E_ch` **唔係純 em**：K 通道 強 IB 殘餘 −3.639 > em +2.265 |
| QN 盲區 | ✅ rank = 2 精確；最大盲區 spread **5.664 MeV**（獨立再證 O124） |
| pairwise 環核（混合對象） | 🟡 5/7 符號對；K、Σ0−Σ+ **反號**；χ² = 1184 |
| pairwise 純 em（K，單位耦合、零擬合） | ✅ +2.0136 vs lattice +2.2646（**−11.1%**） |
| 強 IB 殘餘 `R/X` | ✅ 非普適，∈ [−1.32, +1.36]，spread 2.68（再證 M4 `docs/13`：單一 `⟨O⟩` 不可行） |

## 文檔
- `docs/01_em_operator_layer_round1.md` — 前設問題、事前判準、結果、發現 F1–F5、次輪計劃
- `01_problem_record.md` — 完整記錄（含時序索引）
- `scripts/01_em_operator_constraints.py` — 全部算術（`/usr/bin/python3`，`Fraction` 精確 rank）
- `data/01_em_operator_constraints.json` — checks + 數值
