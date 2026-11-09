# M17 — Session 啟動文件（CMB 大尺度異常）

> **狀態**: ⚪ 已重新登記、未開推（2026-11-09 重寫）。承接時請完整閱讀本檔與 `01_problem_record.md`。
> **重寫原因**: 舊版 §2.1「A3 冷斑 2.6σ 真信號」已被 v10 P4 兩次獨立撤回（`docs/48` 著色依賴、`docs/55` null 錯）。詳見 `01_problem_record §0`。

## 所需 Skills

1. `planck-cmb-map-analysis` / `idcm-cmb-analysis`（核心）— SMICA 讀取、ud_grade、平滑、mask、NESTED/RING
2. `falsifiable-check-design`（含 `references/optimization-statistic-null.md`）— 最優化統計量之 null 設計
3. `idcm-cmb-null-testing` — GRF null vs permutation null 之區別
4. `idcm-audit-verify` — 四件套／傳播／水位
5. `idcm-problem-register` — register 生命週期（跨單元只可 handoff）

## 背景（一句）

v10 P4 嘅 A₅ 節點圖像已被否證（11 異常 **0/11**、A₅ 軸對齊全失敗）。
M17 唔係「換個方式救 A₅」，而係**由觀測異常清單出發**，問 v11 現有結構有無事前可算嘅預言。

## 前置資料

```
專案: IDCM-v11.0--Reality-Architecture/main_problems/m17_cmb_large_scale_anomalies/
docs/     # 空，待承接（每發現 → docs/NN_slug.md，append-only）
scripts/  # 空，待承接
data/     # 空；地圖用 env 連結（SMICA_MAP），勿硬編碼 /tmp

SMICA 全解析度（v10 側，只讀）:
  ../IDCM-v10.0--It-From-Beat/open_problems_fresh/o26_node57_region/data/smica_iqu_2048.fits
  （舊記錄稱 2,013,312,960 B、nside 2048；承接時自行核 byte count 與 npix=50,331,648）

SMICA 已降解析度（v10 側，只讀）:
  ../IDCM-v10.0--It-From-Beat/main_problems/p14_global_routing_table/data/planck_pr3/smica2.fits
```

## 已立事實（不可重推）

```
1. Planck 11 異常存在（觀測事實）
2. A₅ 節點圖像 0/11；A₅ 軸對齊全失敗（P4/docs/33）
3. 冷斑 −134.19 μK、位置 (209°,−57°)、R8=4.17° 係觀測事實，與假說成敗無關
4. 任何「用觀測方向反推 preferred axis」都係循環 —— 禁止（CQ5 已判死）
5. A₅ orbit 有 2 dof（種子）⟹ 「A₅ 對稱軸」唔唯一（P4/docs/72）
6. 60 dof 硬牆：節點級統計上限 <2σ；全天空二點統計更差（P4/docs/73）
```

## 陷阱（已踩過，勿重蹈）

| # | 陷阱 | 後果 |
|:-:|:--|:--|
| 0 | **引用 v10 結論前唔查撤回狀態** | 本次重寫嘅唯一原因。任何 v10 數字先查 `P4/final_result.md` ＋ `docs/48/55/63/67/72/73` |
| 1 | 直接引用 v10 數值 | 違反版本隔離；須 v11 重算 |
| 2 | 全天空無 mask 直接比 Planck 官方值 | 偏差可達 10×（如 S_{1/2}）；CQ1 必須用官方 mask |
| 3 | **permutation null 當宇宙學 null** | 固定節點值 ⟹ 抹掉宇宙方差 ⟹ 假陽性（P4 由 2.51σ 降至 ~1σ）|
| 4 | 用「隨機標籤」作最優化統計量之 null | 假陽性；必須「打亂數值後**重跑最優化**」|
| 5 | 忽略代表元非唯一 | BFS 著色 60 根 → 30 相異集；單一代表元可造出 3σ（P4/docs/48）|
| 6 | 多重比較不校正 | 11 異常 × 檢定量 ≈ 20+ 測試；須 Bonferroni/FDR |
| 7 | 事後選子集 | 必須事前註冊預言與閾值 |
| 8 | NESTED vs RING ordering | NILC/Commander 需 `nest=True` + `reorder(n2r)`；SMICA 係 RING |
| 9 | 單位混淆 | Planck 地圖係 **K**（注入 μK 要 ×1e−6）；`w(θ)` 在 K²（顯示要 ×1e12）|

## 輸入門（開工前必答）

每一子問題開始前，先書面回答「用邊個 v11 常數、為何係佢」。

```
CQ1  唔需要（純觀測）          → 可直接開
CQ2  原初譜 IR 端屬邊個常數？   → 答唔出 ⟹ 🔴 無預言
CQ3  v11 手性／宇稱算子存在嗎？ → 答唔出 ⟹ 🔴 無預言
CQ4  M6 dof=1 能否給 σ²_16？    → 答唔出 ⟹ 🔴 無預言
CQ5  已判死（冇 preferred axis）
```

## 首步（建議順序）

1. 取 Planck 官方 mask（`COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits`），建 `data/maps/`
2. 跑 CQ1：自家 11 異常基線（mask 版），與 Planck 官方值對照表（偏差 <5%）
3. CQ1 通過 → 逐項過輸入門（CQ2/CQ3/CQ4）
4. 三項皆答唔出 ⟹ 按 `01_problem_record §6` 以 negative 閉合

## 邊界

- 本 M 只做觀測層檢驗；群選擇 → M16（已 closed）；DM/DE → M10；味 → M12–M14
- 不得在其他 M/O 目錄寫入任何檔案；跨單元只可用 register handoff
- 每發現 → 獨立 doc（`docs/NN_slug.md`），append-only

---

*建立: 2026-11-05｜重寫: 2026-11-09｜承接: 待*