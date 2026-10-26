# M7 — 強子模基展開：w_k(h) ⟹ ⟨O⟩ = w_1 − w_4

## 狀態: 🔴 closed / verdict **negative**（2026-10-26）

問「強子多夸克態在 `K_chain ⊗ K_tree` 基下之展開 `w_k(h)` 有無零參數結構式」。命題：味屬鏈模（`k∈{1,4,5}`）⟹ `ūu − d̄d ⟶ P_{k=1} − P_{k=4}` ⟹ `⟨O⟩_h = w_1(h) − w_4(h)`，並須零參數重現 BMW 三錨 `0.891 / 1.498 / 2.293`（`n_s = 0/1/2`，單調升）。

## 結果一句話

**`w_k(h)` 存在但退化為恒等**（整數定理 ⟹ `w_k` ≡ 價夸克數）；**`⟨O⟩ = w₁ − w₄` 雙重死**（定義層：v9 sector 表不含 u、d；數值層：sector 讀法 ⟹ `⟨O⟩≡2`、價夸克讀法 ⟹ `⟨O⟩≡1`，後者 `p = 3.96e−4`＝3.36σ）；
★ 最重要：**M2 已證「質量 = 加性泛函」⟹ IDCM 已註冊結構強制 `⟨O⟩ ≡ 1`，與 BMW 三錨衝突 3.36σ**，須上游裁定（M2 加性適用範圍／三體束縛，已 handoff）。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 ＋ §十 執行 timeline（8 輪） | ✅ |
| `final_result.md` | 判決摘要（T1 🟡／T2 🔴／T3 🔴；禁區；handoff） | ✅ |
| `SESSION_START.md` | 起手式（原文保留，頂部加完成標記） | ✅ |
| `docs/10_reproduction_gate.md` | R0 重現閘 **19/19**（含 mpmath 30 dps 打磨、weight 約定雙軌） | ✅ |
| `docs/11_definition_freeze.md` | R1 定義凍結（hash `883e91d5…92421`）＋ 語料 grep 證據 | 🔒 |
| `docs/12_r5_verdict.md` | R5 判死閘（K1–K8）＋ 4 項結構結果 ＋ 審計項 | 🔴 |
| `docs/13_file_deletion_audit.md` | 13:33 刪除事件審計：21 檔零缺失、6 檔離碟（全屬並行 session、全未 commit、6/6 副本在 scratch） | ✅ |
| `scripts/01_reproduction_gate.py` | 零 import 獨立重建兩體骨架 | ✅ |
| `scripts/02b_r5_mode_content_nogo.py` | 整數定理、χ²、非可加性、跨 M 核（v1 封存 `deprecated/`） | ✅ |
| `scripts/03_grep_evidence.py` | 概念詞證據 ＋ 凍結 hash | ✅ |
| `data/` | 3 JSON ＋ 全部 run log（含 v1 log 存證） | ✅ |

## 依賴

M4（`docs/28` 三錨／`docs/31` 密度插入否證／`docs/32` 發還原文）· M5（兩體徑向骨架、contact 項）· **M2（`docs/01` §2 加性質量泛函——本輪關鍵跨 M）** · O210（`𝒯`）· O216（色 = 嵌入層標籤）· v9 §D.1（味 sector 原文）

## 別名 ⚠️

調律者原話「開新m6」之時，**M6 = 尺度自由度定理（`m6_scale_dof_theorem`）已存在** ⟹ 本題開為 **M7**（`m7_hadron_model_basis`）。

## 註冊

pid `M7`（created 2026-10-26；closed 2026-10-26 verdict `negative`）
handoff ×3：→ **M2**（加性泛函適用範圍）· → **M5**（weight 約定 `N_dv²` vs `v²`）· → **M4**（`docs/29` p 值 3.96e−4；M1 後驗標記）

## 未閉合（下游）

1. 三體束縛 `E_bind^{(3)}`（唯一實質活口，需新題）
2. M2 加性泛函是否覆蓋複合態（決定 `⟨O⟩` 問題係數據讀法抑或結構缺口）
3. Q4：s/d、u/c 之 sector 內細分（v9／P10 側）