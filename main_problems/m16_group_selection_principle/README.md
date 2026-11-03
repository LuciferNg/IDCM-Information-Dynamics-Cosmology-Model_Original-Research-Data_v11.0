# M16 — 群選擇原理（Group Selection Principle）

**狀態**: 🟢 已閉合（2026-11-03）｜ **version** v11.0 ｜ **性質** 公理層缺口 → 已歸約至上游

---

## 定位

`M15` 已把宇宙預算分母由「凑數」升級為「A₅ 之 irrep 恆等式」：

```
Σd² = Σd + Σd(d−1)   ⟺   60 = 16 + 44      （d = [1,3,3,4,5]）
```

本 O 處理下一步：**為何係 A₅？**

## 結論

```
C(G) := [∏(相異非平凡 irrep 維數) = |G|] ∧ [Σ(相異非平凡 irrep 維數) = 12]
C(G) ⟺ G ≅ A₅          （覆蓋所有階；零自由度、零新增假設）
```

**A₅ = 條件導出，非公理輸入。** 缺口搬遷至上游（`MR1–MR3`，全部屬 `O171`）。

| 已閉 | 出處 |
|:--|:--|
| `C(G) ⟺ DNT = {3,4,5} ∧ \|G\| = 60`（算術引理） | `docs/01` |
| 92,570 群全域掃描 → 唯一 A₅ | `docs/02` |
| 修飾詞承重 ＋ look-elsewhere | `docs/02 §3–4` |
| `∏ = \|G\|` 之真內容 = 排除 abelian 膨脹（11 群族） | `docs/03 §3` |
| 判詞 ＋ R5 ＋ 發還 M3 | `docs/04` |

## 目錄

```
m16_group_selection_principle/
├── 01_problem_record.md     問題、缺口、參考集、推導時序
├── SESSION_START.md         承接指引（已閉合；殘項指向 O171）
├── README.md                本檔
├── final_result.md          🟢 判詞（一句話 ＋ 關鍵式 ＋ 審計 ＋ 已否證）
├── docs/01_r5_and_arithmetic_reduction.md
├── docs/02_global_scan_and_modifiers.md
├── docs/03_criterion_source_and_verdict.md
├── docs/04_verdict_r5_and_status.md
├── scripts/00_gap_named_degrees.g   具名群維數（GAP 權威）
├── scripts/01_arithmetic_reduction.py
├── scripts/02_gap_full_scan.g ｜ 02b_chunk_driver.py ｜ 02c_order_counts.g
├── scripts/03_variants_and_lookelsewhere.py
├── scripts/04_gap_alternative_principles.g
├── scripts/05_candidate_principles.py
├── scripts/gap_run.sh ｜ run_02_parallel.sh
└── data/                    全部 run 記錄（含 92,570 群掃描原始檔）
```

## 使用方式

1. 想知結論 → `final_result.md`；
2. 想核證明 → `docs/01`（算術）＋ `docs/03`（定理）；
3. 想核掃描／look-elsewhere → `docs/02` ＋ `data/run_02_scan_all_500.txt`；
4. 想接殘項 → `SESSION_START.md §四`（MR1–MR3，全屬 O171）。

## 環境備註

GAP 4.14.0 + SmallGrp 1.5.4 由 Debian deb 解包至
`/home/wsl/.hermes/cache/scratch/gapdebs/root`（無 root）；`scripts/gap_run.sh` 為啟動器。
掃描以 11 進程並行分片執行（`scripts/run_02_parallel.sh`），片段按其階之子區間切分。