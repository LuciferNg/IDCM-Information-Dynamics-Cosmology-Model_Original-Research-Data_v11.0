# M6 SESSION_START — 尺度自由度定理：複合投影 ⟹ dof=1

<!-- registry
pid:     M6
title:   尺度自由度定理：複合投影 ⟹ dof=1（絕對標度非獨立自由度）
state:   closed
verdict: closed
flags:   -
owner:   M6
created: 2026-10-25
updated: 2026-10-26
prereq:  O38, O212, M2, O213, O121
source:  調律者2026-10-25「開M6, 只開不推」；理據源於 2026-10-25 對話（絕對標度／複合投影推論）
-->

> **Session**: 2026-10-25 開檔 → **2026-10-26 推導完成（本 session）**
> **狀態**: 🟢 **已閉合** — T1 ✅ · T2 `dof_scale = 1` ✅ · P4 反例 0 ✅ · T3 → O38 G2（handoff）
> **結果摘要**: `final_result.md`｜**推導**: `docs/01`–`docs/04`（40/40 PASS）
> **來源**: 調律者 2026-10-25「開M6, 只開不推」；2026-10-26 指定本 session 為 M6 專用
> **前置**: O38（docs/15 §4 單位制測試；docs/28 穹頂判據；docs/38 G1/G2 裁決）· O212（錨定理 dof=1）· M2（標度盲引理）· O213/O121（宏觀端錨）
> **位置**: `main_problems/m6_scale_dof_theorem/`

---

## 一、接手起點（第一步讀乜）

1. `O38 docs/15` §4 — 單位制測試：同一個量三個數字 ⟹ 絕對數值 = unit artifact
2. `O38 docs/28` §2–§4 — 穹頂判據：兩端密集、中間 30.9／25.1 dex 空洞、樹階梯命中 1/4
3. `O38 docs/38` §1–§4 — G1 關閉（M₀ 數值 = 偽命題）／G2 保留（4.095e38 計數身份）／「每級對應物」假設應放棄
4. `O212 docs/01` §1 — dof = 1 齊次性論證（錨唯一）
5. `M2 docs/01` §7–§9 — 標度盲引理、`Σ√γ=cot(π/176)+1`、振幅讀法份額 0.430≠Ω_m
6. `01_problem_record.md` §二 — 待證 T1–T3（未證，勿當結果引用）

## 二、唔准做（禁區）

- 🔴 禁 `φ^exp` 族（O136 已判過度吻合否證）
- 🔴 禁重推 O38 G1（數值可推 = 偽命題，已關閉）與 G2（生成器 → 屬 O38，本 M 只可 register handoff）
- 🔴 禁把「絕對標度無意義」講成「錨可刪」——O212 dof=1 仍然強制
- 🔴 禁在本 M 寫入其他 M/O 目錄（跨 O 只准 register handoff）

## 三、五處分工

| 檔 | 角色 |
|:--|:--|
| `README.md` | 對外導引（非賬本） |
| `01_problem_record.md` | 推導鏈（本 M 現時只有問題陳述，無推導） |
| `SESSION_START.md` | 本檔 = 起手式 |
| `final_result.md` | 未建（未開始，無結果可摘要） |
| `problem_register/` | 生命週期（pid `M6`） |