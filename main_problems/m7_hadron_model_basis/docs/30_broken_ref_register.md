# M7 docs/30 — 斷鏈引用登記表（歷史／跨記錄）

> 機器可讀：`data/30_broken_ref_register.json`（7 條）｜閘門：`scripts/07` **D1**
> 觸發：調律者第 3 項「檢查所有文檔數據及腳本引用數字無漂移及斷鏈」。

## 一、原則

append-only ⟹ **歷史引用不改寫**。凡文檔內指向不存在之 `scripts/`／`data/` 檔者，於本表登記其身份；`07` D1 只放行本表所列，其餘缺檔一律 🔴。
**新文檔規則**：引用現行檔一律用現名；不得再引入舊名。

## 二、登記（7 個名字，15 個出現位置）

| 舊引用 | 類 | 身份 | 出處 |
|:--|:--|:--|:--|
| `scripts/01_repro_gate.py` | 歷史離碟 | 舊名（12:42 前），git 從未 commit，scratch `d577a467dd7f…`；現名 `01_reproduction_gate.py` | `docs/13` §三 |
| `scripts/02_definition_freeze.py` | 歷史離碟 | 舊名腳本，scratch `2cd41c46f7dc…`；內容遷為 `docs/11` | `docs/13` §三 |
| `data/01_repro_gate.json` | 歷史離碟 | scratch `0ad2f06e3d89…`；現名 `01_reproduction_gate.json` | `docs/13` §三 |
| `data/02_definition_freeze.json` | 歷史離碟 | scratch `30369fd256a4…` | `docs/13` §三 |
| `scripts/02_t2_scale_dof_theorem.py` | 跨記錄 | **M6** 之未追蹤 WIP（語料定義漂移之源） | `docs/14` §六 |
| `data/02_t2_scale_dof_theorem.json` | 跨記錄 | 同上（M6） | `docs/14` §六 |
| `scripts/16_define_d6.py` | 跨記錄 | **M10** D6 定義生成腳本（存在於 `m10_*/scripts/`） | `docs/25` §三 |

## 三、判詞

- 現行推導鏈（`01–20` 腳本 → `10–29` 文檔 → JSON）**零斷鏈**：21 腳本全部有 log ＋ JSON（`10` 用 `watermark.json`），20 文檔全部有 README 落點。
- 上述 15 處皆**歷史記錄或跨記錄引用**，非現行鏈斷裂 ⟹ 登記而不改寫。
- 數字漂移：`04` 文檔↔輸出 **220/220**；`07` B 層獨立重算（B1–B13）全 PASS。
