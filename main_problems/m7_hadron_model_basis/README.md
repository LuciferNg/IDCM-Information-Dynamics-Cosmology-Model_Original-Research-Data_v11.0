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
| `docs/14_rerun_verification.md` | 全 4 腳本重跑（JSON byte-identical）＋ 文檔數字 **88/88 PASS**；發現 3 個記錄層錯（已 append-only 更正） | ✅ |
| `docs/16_gap_closure.md` | 未推盡審計（R3／MC null／剛性／Q4／7-12）4 項全補；**撤回 `docs/12` K6「sector ⟹ ≡2」**；失配在斜率 4.19σ | ✅ |
| `docs/15_additivity_scope_verdict.md` | ① 加性適用範圍（容量／飽和攻擊）⟹ **NO-GO**；② **不獨立**（M2 唯一性禁止獨立束縛項）⟹ **單一 no-go**（§九 更正） | 🔴 |
| `scripts/01_reproduction_gate.py` | 零 import 獨立重建兩體骨架 | ✅ |
| `scripts/02b_r5_mode_content_nogo.py` | 整數定理、χ²、非可加性、跨 M 核（v1 封存 `deprecated/`） | ✅ |
| `scripts/03_grep_evidence.py` | 概念詞證據 ＋ 凍結 hash | ✅ |
| `scripts/04_doc_consistency.py` | 文檔↔輸出檢查器（**111** 條 claim；ulp 由文檔精度推，零可調容差） | ✅ |
| `scripts/06_gap_closure.py` | 未推盡補推 G1–G6（**12/12**）：MC null 10⁷、家族可辨識性、Q4 判死、剛性、7/11 重現 | ✅ |
| `scripts/05_additivity_scope.py` | ① 容量／強迫共享搜索 ＋ ② pairwise-only 量級核（**10/10**） | ✅ |
| `data/` | 5 JSON ＋ 5 run log（含 v1 log 存證） | ✅ |

## 依賴

M4（`docs/28` 三錨／`docs/31` 密度插入否證／`docs/32` 發還原文）· M5（兩體徑向骨架、contact 項）· **M2（`docs/01` §2 加性質量泛函——本輪關鍵跨 M）** · O210（`𝒯`）· O216（色 = 嵌入層標籤）· v9 §D.1（味 sector 原文）

## 別名 ⚠️

調律者原話「開新m6」之時，**M6 = 尺度自由度定理（`m6_scale_dof_theorem`）已存在** ⟹ 本題開為 **M7**（`m7_hadron_model_basis`）。

## 註冊

pid `M7`（created 2026-10-26；closed 2026-10-26 verdict `negative`）
handoff ×3：→ **M2**（加性泛函適用範圍）· → **M5**（weight 約定 `N_dv²` vs `v²`）· → **M4**（`docs/29` p 值 3.96e−4；M1 後驗標記）

## 已封（2026-10-26 第三輪）

1. ~~M2 加性泛函是否覆蓋複合態~~ ⟹ **已答：覆蓋**（每模容量 6 > `N_q` 3 ⟹ 永不飽和）→ `docs/15`
2. ~~三體束縛 `E_bind^{(3)}`~~ ⟹ **不獨立**（M2 唯一性窮盡 ⟹ 冇束縛能槽位；量級核降為旁證）→ `docs/15` §九

## 未閉合（下游／新結構需求，未開題）

1. ~~三體味奇頂點 `λ₃`~~ —— **撤回**（加入即放棄 M2 唯一性）
2. **連續夸克質量參數**（唯一活口） —— 現譜離散 ⟹ Feynman–Hellmann `⟨O⟩` 本身唔係 IDCM 內在量
3. Q4：s/d、u/c 之 sector 內細分（v9／P10 側）
4. M4：`⟨O⟩` 程序須重審（對應關係問題）