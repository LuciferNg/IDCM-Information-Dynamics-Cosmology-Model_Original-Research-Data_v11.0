# M29 — E₈ 分塊全枚舉：3,679,075,400 態（R8 工程落地）

<!-- registry
pid:     M29
title:   E₈（V^{9→7}）3,679,075,400 態之分塊全枚舉 —— 未剪枝 kernel
state:   registered
verdict: -
flags:   opened-only; no-derivation; engineering
owner:   M29
created: 2026-11-27
updated: 2026-11-27
prereq:  M14 §四十四（S4 分塊 bit-identical 機制驗證）＋ R8（資源判死改判可行）
source:  調律者 2026-11-27「各自分派到新M另開session執行」；缺口出自 M14 docs/57 §十三 路線 Ⅳ
-->

> **狀態**: ⚪ 已註冊、未開推（只開不推）
> **定位**: 純工程 ＋ 正確性；**唔改判準，唔剪枝**。

---

## §〇 歸屬審計（fence）

| 準則 | 事實 | 裁定 |
|:--|:--|:--|
| ① 機制已驗 | 分塊 memmap 累加 vs 稠密 **bit-identical**、`max|diff| = 0.0`；`2^24` 塊 = 0.125 GiB ⟹ 220 塊 | ✅ 機制唔需重證，只需在目標格落地 |
| ② 硬規 | 「**不剪枝**」：`V^{9→7}` 必用未剪枝 kernel | ✅ 寫入禁區 |
| ③ 環境 | GTX 1060 3 GB；54 GB RAM；12 核；float64 27.411 GiB（29.43 GB） | ⚠️ 須分塊外存 |

**裁定**：開 M29。

---

## §一 問題陳述

| # | 問題 | 為何必要 |
|:-:|:--|:--|
| Q1 | 3,679,075,400 態之未剪枝 kernel 可否在 float32 13.706 GiB／float64 27.411 GiB 下全枚舉？ | R8 兩格 blocked 之解 |
| Q2 | 分塊（220 × 0.125 GiB）與稠密是否 bit-identical（目標格重驗）？ | 正確性閘 |
| Q3 | 全枚舉結果是否改變 `V^{9→7}` 之既有結論？ | 否則只係工程收據 |

---

## §二 待證命題（🔴 未證；判準跑前寫死）

| # | 待證 | 失敗判準 |
|:-:|:--|:--|
| T1 | 目標格分塊 ≡ 稠密（bit-identical，`max|diff| = 0`） | 任何非零差 ⟹ 🔴 |
| T2 | 未剪枝（kernel 全集；不得以 0 權重跳過） | 剪枝 ⟹ 整輪作廢 |
| T3 | 兩跑 stdout byte-identical；耗時只落 stderr | 違 ⟹ 作廢 |
| T4 | 結果若改變舊結論，須以 append-only 附註記錄 | 改寫舊文 ⟹ 作廢 |

---

## §三 已知前置（附出處）

| 項 | 值 | 出處 |
|:--|:--|:--|
| 態數 | 3,679,075,400 | M14 §四十四 |
| 記憶體 | float32 13.706 GiB／float64 27.411 GiB（29.43 GB） | M14 §四十四 |
| 分塊 | `2^24` = 0.125 GiB ⟹ 220 塊 | M14 §四十四 |
| 機制 | 分塊 vs 稠密 bit-identical；`max|diff| 0.0` | M14 §四十四 |

---

## §四 開放項（P）

| # | 項 | 狀態 |
|:-:|:--|:-:|
| P1 | 目標格 bit-identical（T1） | 🔴 |
| P2 | 全枚舉（T2–T3） | 🔴 |
| P3 | 與舊結論比對（T4） | 🔴 |

---

## §五 攻擊路線（未執行）

| # | 路線 | 成本 |
|:-:|:--|:-:|
| R1 | 先在小態數（如 `V^{3↔5}` 58,442,310）重驗 bit-identical | 低 |
| R2 | 目標格分塊落地（220 塊、checkpoint 續跑） | 高 |
| R3 | 全枚舉 ＋ 兩跑 byte-identical | 高 |

---

## §六 禁區

```
✗ 剪枝（V^{9→7} 必未剪枝）
✗ float32 化正典（float64 正典；float32 只作對照）
✗ 改容差／降 dps
✗ 改寫舊結論（append-only 附註）
```

---

## §七 驗證標準

bit-identical／未剪枝／兩跑同碼／耗時落 stderr。

---

## §八 狀態

```
state = registered（只開不推）｜verdict = -｜scripts 0｜docs 0｜data 0
```

---

## §九 相關文件

- M14 `docs/49` §四十四（S4 分塊驗證）
- M14 `scripts/oneoff/probes_m14/probe_s4_chunk_validate.py`
