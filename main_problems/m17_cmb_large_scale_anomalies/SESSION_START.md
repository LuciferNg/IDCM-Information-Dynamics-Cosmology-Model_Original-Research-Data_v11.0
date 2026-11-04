# M17 — Session 啟動文件（CMB 大尺度異常）

> **狀態**: ⚪ 已註冊、未開推。承接時請完整閱讀本檔與 `01_problem_record.md`。

## 所需 Skills

1. `idcm-cmb-analysis`（核心）— SMICA 下載/讀取/平滑/節點取樣、mask 處理、NESTED pitfall
2. `idcm-audit-verify` — 四件套／傳播／水位
3. `idcm-problem-register` — register 生命週期
4. `falsifiable-check-design`（含 `references/optimization-statistic-null.md`）— 最優化統計量之 null 設計

## 背景（一句）

v10 `P4` 嘅 A₅ 節點圖像已被否證（Planck 11 異常命中 0/11、A₅ 軸對齊全失敗）。
M17 唔係「換個方式救 A₅」，而係**由觀測異常清單出發**，問 v11 現有結構有無事前可算嘅預言。

## 前置資料

```
專案: IDCM-v11.0--Reality-Architecture/main_problems/m17_cmb_large_scale_anomalies/
docs/     # 空，待承接
scripts/  # 空，待承接
data/     # 空；SMICA 地圖用 env 連結，勿硬編碼 /tmp

SMICA map: ../IDCM-v10.0--It-From-Beat/open_problems_fresh/o26_node57_region/data/smica_iqu_2048.fits
           2,013,312,960 B（用完請自行驗證 byte count 與 nside=2048、npix=50,331,648）
```

## 已立事實（不可重推）

1. A₅ 節點假說 🔴 否證（v10 register `P4` verdict = negative）
2. A₅ 群自生成 null：圖拓撲隨種子方向變（30/60/90/120 邊），與冷斑同構僅 12/100
3. 冷斑本體 −134.19 μK、位置 (209°,−57°)、尺度 R8=4.17° 屬**觀測事實**，與假說成敗無關
4. 任何「用觀測方向反推 preferred axis」都係循環 —— 禁止

## 陷阱（預先提醒）

| # | 陷阱 | 後果 |
|:-:|:--|:--|
| 1 | 直接引用 v10 數值 | 違反版本隔離；須 v11 重算 |
| 2 | 全天空無 mask 直接比 Planck 官方值 | 偏差可達 10× (如 S_{1/2})；CQ1 必須用官方 mask |
| 3 | 用「隨機標籤」作最優化統計量之 null | 假陽性（v10 P4 曾得 100% vs 實際 5.3%） |
| 4 | 多重比較不校正 | 11 異常 × 多檢定量 ≈ 20+ 測試；須 Bonferroni/FDR |
| 5 | 事後選子集 | 必須事前註冊預言與閾值 |
| 6 | NESTED vs RING ordering | NILC/Commander 需 `nest=True` + `reorder(n2r)`；SMICA 係 RING |

## 首步（建議順序）

1. 讀 v10 `P4/docs/33`，抄下 11 異常之官方定義與顯著度
2. 取 Planck 官方 mask（`COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits`），建 `data/maps/`
3. 跑 CQ1：自家 11 異常基線（mask 版），與官方值對照表
4. CQ1 通過 → 才開 CQ2–CQ4；CQ5 標 🔴 阻塞待 M16

## 邊界

- 本 M 只做觀測層檢驗；群選擇 → M16；DM/DE 身份 → M10
- 不得在其他 M/O 目錄寫入任何檔案；跨單元只可用 register handoff
- 每發現 → 獨立 doc（`docs/NN_slug.md`），append-only

---

*建立: 2026-11-04｜承接: 待*