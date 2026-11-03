# M1 docs/07 — 入站 handoff 審計與裁定（2026-10-29）

**指令**：檢查並處理推送至 M1 register handoff｜**方法**：`register-handoff-handling` skill
**範圍**：`problem_register/register.jsonl` 全 3565 行（append-only 正本）＋ `register.sqlite` event 表交叉核

---

## §1 掃描

| 掃法 | 命中 |
|:--|:--|
| `target_pid == "M1"`（大小寫不敏感） | **1** |
| `target_pid` 為 M1 子 pid（`M1/…`） | 0 |
| `ref` 或 `note` 含 `m1_mass_file_schema` 而 `pid ≠ M1` | 0 |
| sqlite `event` 表 `upper(target_pid)='M1'` | 1（與 jsonl 同一條） |
| M1 目錄內非本 M 產生之檔（外部投遞） | 0 |
| 其他 register 檔（v11 本地） | 無 |

> 註：`status --pid M1` 之 `incoming (1)` 即下表唯一一條；子 pid 隱形問題在本 M 不適用（無子 pid）。

---

## §2 入境清單（唯一一條）

```
ts      2026-10-20T13:10:06+08:00
kind    merged          from pid m1 → target M1
ref     （無）
note    pid 跟隨 O 系列大寫化 (2026-10-18): m1 → M1
```

---

## §3 裁定

| 項 | 判 |
|:--|:--|
| 性質 | **別名合併事件**（`merged`），非請求型 handoff：無 ask、無 ref、無輸入供給 |
| 裁定 | **ACK closed-moot** —— 無需接收方任何動作；不需回覆來源（來源即自身舊 pid） |
| 對 M1 現狀影響 | 無。`v11.0/M1` 單一 pid、state=partial、無重複／無孤兒 alias（`alias` 事件 2026-10-20T13:10:15 已解析） |
| 行動 | 不改任何文件；僅本審計記錄 |

---

## §4 出境未覆（附帶核，非入境）

| 發出 | 目標 | 內容 | 現況 |
|:--|:--|:--|:--|
| 2026-10-28T11:56:30 | M4 | 請求把 M4 `docs/13`、`docs/15` 之「機制 M1／機制 M2」改名 K1／K2（與問題編號撞名） | **未覆**：M4 `docs/15` 仍寫「機制 M1/M2」（17:49 版）；M4 全日忙於 audit 與技術債，未處理此項 |
| 2026-10-28T14:12:34 | M10 | 供給型：ρ 引擎＋coverage 17.040 %，支撐 W5「質量值 dof=0」 | **未覆**（供給型，無 ask；M10 已 partial 收尾） |

**裁定（兩條皆判「不追」）**

1. **M4 命名**：障礙只在引用時嘅歧義，可由 M1 **單邊消歧**解決 ⟹ 本 M 此後引用一律寫
   **「M4 `docs/15` §T3『機制 M1』」**（帶來源與節號），並在文內標明與問題編號 `v11.0/M1` 無關。
   請求由「要求」降級為「已記錄之建議」，不再重發（避免 register 噪音）。
2. **M10 供給**：無 ask ⟹ 無待答；M10 若日後回覆，於新 docs 段 append 對帳即可。

---

## §5 一致性與完整性（四項核對）

| # | 核對 | 結果 |
|:-:|:--|:--|
| ① | register 現狀 | `v11.0/M1` state=partial、verdict=partial、flags `-` |
| ② | jsonl 雙向 | 自己事件 9 條（registered×2／merged／alias／note×3／progress／handed_off×2 見 §4）；入境 1 條（§2） |
| ③ | 漂移 | `register.py check`：**HARD DRIFT 0**、**BLOCK-DRIFT 0** |
| ④ | 文檔狀態一致 | README「🟡 partial」／`01_problem_record.md` 塊 `state: partial`／`SESSION_START.md` 塊 `state: partial` ⟹ 三處一致 |
| ⑤ | 檔案完整性 | `scripts/` 7 支 ↔ `data/` 6 個 json（`05_doc_consistency` 只出 stdout，無 json，屬設計）；`docs/01–07` 與 README 表對齊；`final_result.md` **未建**（partial 階段，close 時建；`check` 未報缺） |

---

## §6 結論

```
入境待裁定 = 0（唯一一條 merged 判 closed-moot）
新增工作項 = 0
出境未覆 = 2（皆已裁定不追：M4 改單邊消歧、M10 屬供給型）
來源目錄零寫入（M4／M10 未動一字）
```

---

## §7 更正（2026-10-29，append-only）

**M4 已於 2026-10-29T01:49:09 回覆**（`docs/46` E3，register note）：**接受改名 機制M1/M2 → K1/K2**，並作**範圍修正** —— 只 `docs/15` 有 6 處，`docs/13` **零命中**（M1 原請求範圍過寬）；`docs/15` 已加 append-only 附註。

**M1 獨立核（唯讀）**

```
grep -c "機制 M1|機制 M2" m4/docs/15_exponent8.md = 7 行（5/24/26/43/87/88/113；113 係新附註本身 ⟹ 原文 6 處 ✓）
grep -c "機制 M1|機制 M2" m4/docs/13_mdmu_second_instance.md = 0 ✓（只有 M_EW 一詞）
```

**裁定更正**

| 原判（§4） | 更正後 |
|:--|:--|
| M1→M4 改名請求「未覆 ⟹ 降級為建議、不再追」 | **已被接受**（M4 `docs/46` E3、`docs/15` 附註）⟹ 由「建議」升回**已結**；M1 引用一律寫 **K1／K2** |
| 「M4 `docs/13`、`docs/15` 撞名」 | **M1 原請求事實錯誤**：`docs/13` 零命中 ⟹ 撤回該部分，只 `docs/15` 成立（M1 自查錯，append 記錄，原文不刪） |
| M1→M10 W5 供給「未覆」 | 仍屬供給型（無 ask）；M10 已於同日自行重推（`docs/19`），M1 亦已補發 W5 加强版 handoff（`docs/08`） |

**後續引用規則**：M1 全部文件引用 M4 該兩機制時用 **K1**（`γ₁` per-level）／**K2**（`dim su(N_gen)=8`）；舊稱只可出現在引述原文時。
