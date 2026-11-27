# M29 SESSION_START — E₈ 分塊全枚舉

<!-- registry
pid:     M29
title:   E₈（V^{9→7}）3,679,075,400 態分塊全枚舉（未剪枝）
state:   registered
verdict: -
flags:   opened-only; no-derivation; engineering
owner:   M29
created: 2026-11-27
prereq:  M14 §四十四
source:  調律者 2026-11-27「各自分派到新M另開session執行」
-->

> **狀態**: ⚪ 已註冊、未開推（只開不推）

## 一、機制唔需重證

分塊 vs 稠密已 bit-identical（`max|diff| 0.0`）。本 M 只做**落地**。

## 二、第一步

```
R1  小態數重驗（V^{3↔5} 58,442,310）
R2  目標格 220 塊、checkpoint 續跑
R3  全枚舉、兩跑同碼
```

## 三、硬規

```
✗ 不剪枝｜✗ float32 化正典｜✗ 改容差｜✗ 改寫舊結論
```

## 四、資源

```
12 核｜54 GB RAM｜GTX 1060 3 GB｜float64 27.411 GiB ⟹ 必須外存分塊
```
