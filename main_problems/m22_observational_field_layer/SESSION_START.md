# M22 — SESSION_START（接手文件）

> **建立**: 2026-11-14｜**狀態**: ⚪ 已註冊、未開推
> **本 session 性質**: 只開不推 —— 零 scripts、零 docs 內容、零數據。

---

## §一 接手起點

1. 本檔（章程 + 禁區 + 重現閘）
2. `01_problem_record.md`（§〇 fence、Q1–Q4、T1–T4、R1–R5、§六 禁區）
3. `README.md`
4. M17 `01_problem_record.md` §3（輸入門 CQ2–CQ5）、§4 判準、§6 收口條件
5. M6 `final_result.md`（dof = 1 定理）
6. M16 `final_result.md`（A₅ 唯一解）
7. v10 `P4 docs/72`、`docs/73`、`docs/48`（只作約束，唔作輸入）

---

## §二 重現閘（硬性）

```
① 核 M17 §3 五個輸入門現狀（CQ1 可做、CQ2–CQ4 落空、CQ5 no-go）
② 核 M6 dof=1 定理陳述（唔重推）
③ 核 v10 P4 docs/73 之 60 dof 硬牆數字（−0.01× 靈敏度）
④ 唔准重跑 M17／M6／M16 腳本；唔准跑 v10 P4 腳本
```

---

## §三 特別紀律：投影 ≠ 解釋

```
本 M 只答「結構常數喺天空上點顯示」，唔答「11 異常解唔解得通」。
判死／判生屬 M17。M22 交付「輸入門通過與否」。
```

---

## §四 禁區（不得重推 / 不得復用）

```
✗ 60 節點直讀作 CMB 統計來源      P4 docs/73：靈敏度 −0.01×
✗ A₅ orbit 身份作鑑別             P4 docs/72：2 dof
✗ 冷斑 2.6σ / p=0.0050            M17 §2 已撤回
✗ 正交統計量合併至 5σ             P4 docs/63/69：上限 < 2σ
✗ 由觀測反推 v11 常數             M17 §3 明禁
✗ v10 W-field 數值作結構輸入       v10 已凍結，CMB 結論多撤回
```

---

## §五 與相鄰 M 之分工

| M | 做咩 | M22 唔做咩 |
|:-:|:--|:--|
| M17 | 11 異常之自家測量、預言檢驗、判死 | 唔判 M17 五項子問題生死 |
| M6 | `dof=1` 定理 | 唔重推 dof |
| M16 | 群選擇 | 唔重推 A₅ 唯一性 |
| M3 | Ω 三分量 | 唔重算比例 |
| M12 | CY₃ 投影（唔係天空投影） | 唔混用「投影」一詞（本 M 指結構→場） |

---

## §六 交接／依賴

- **上游**：M17（輸入門需求）、M6、M16、M3、M12
- **下游**：M17 CQ2／CQ3／CQ4
- **跨 M 交接只限 register**；禁在他 M 寫檔。

---

## §七 本 session 交付記錄

```
2026-11-14  開檔：01_problem_record.md + SESSION_START.md + README.md
            資料夾 m22_observational_field_layer/{docs,scripts,data} 建空
            零推導、零 scripts、零 data
2026-11-20  第一輪執行完（partial）：重現閘 35/35；docs/01 判準凍結
            scripts/01–05（34/20/9/8/9 PASS，合計 115 PASS／0 FAIL）
            docs/02 投影算子（T1 ✅）；docs/03 統計量帶；docs/04 三輸入門答覆
            final_result.md；新增禁區 5 條；register note →
            判詞：Q1 ✅／Q2 ✅形式／Q3 🟡／Q4 🔴；T1 ✅／T2 🟡→✅形式／T3 🟡／T4 🔴
2026-11-20  R6 殘項攻堅（docs/05；scripts/07–10 = 14/15/6/5 PASS）
            MR0 收窄（{3,3′} 二值、惰性）｜MR1 加強（NA4 Ω 不可落格；BR1′/BR1″）
            MR2 閉（負）：反極點恆等式 ⟹ 手性不可測（無條件）
            MR3 閉：ΛCDM 真譜 —— 缺損 40.18%、f_crit 0.23160／0.11948
            合計 196 PASS／0 FAIL
2026-11-20  R7：MR0 終局（docs/06；scripts/11 = 14/14 PASS）
            形式證死：線性化 ⟹ {3,3′}；度量 1 維 ∝ I（dof_scale=1）；自由軌道 2000/2000
            存在不可導出：|Aut 𝒯| = 352 = 2⁵·11、5 ∤ 352、A₅ 單群 ⟹ 同態平凡
            ⟹ 立 bridge axiom B1；全部結果 conditional on B1；合計 213 PASS／0 FAIL
```
