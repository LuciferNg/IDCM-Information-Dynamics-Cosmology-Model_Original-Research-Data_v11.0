# M19 — 模式指數註冊（mode index registry）

## 狀態

🟢 **`closed`（scoped：9 模式集 × F1 frame）** — 2026-11-20 收尾審計｜水位 **97/97**（driver rc=0）｜零跨 M 寫入

## 一句話

`k_q` 係**混合量**：模式集 ＋ 輕子類別有結構出處（4.392 bit），逐粒子值全由質量 rank 餵入（14.077 bit）——`s` 三讀結構零判別力、`c/b/t` 純結構規則 0/21 命中、`S4` 阻塞性質更正為「可算但不可判」。

## 為何要開

M12 `final_result.md` L201：

> 阻塞：`c`、`b` 之模式指數在 v11 **未註冊**（`u→2`、`d→3` 已裁；`s` 三讀 4/5/7 未定）
> ⟹ 數值評估即後驗指數選擇，依禁區不可作判據。

M14 `final_result.md` L173：

> **S4（非 ℤ₂ 定向）阻塞於「c／b 模式指數未註冊」**、**S5 無可算載體**。

兩處獨立確認**同一阻塞點** ⟹ 唔係「某 M 未做完」，係**登記冊缺一類物件**。

## 現況表

| 夸克 | 指數 | 現況 |
|:--|:-:|:--|
| u | 2 | ✅ 已裁（v10 `O137→O233`） |
| d | 3 | ✅ 已裁（同上） |
| s | 4/5/7 | 🔴 三讀未定 |
| c | — | 🔴 未註冊 |
| b | — | 🔴 未註冊 |
| t | — | 🔴 未註冊 |

## 檔案

| 檔 | 內容 |
|:--|:--|
| `01_problem_record.md` | §〇 歸屬 fence／Q1–Q4／T1–T4／P1–P4／R1–R5／禁區／驗證標準 |
| `SESSION_START.md` | 接手起點、重現閘、禁區、分工、交接 |
| `README.md` | 本檔 |
| `docs/` `scripts/` `data/` | docs 8 份、scripts 9 支 ＋ l1–l3 ＋ run_chain.sh、data 7 JSON ＋ 13 log |

## 開放問題

權威清單只有一份：**`docs/08_open_problems_inventory.md`**（OA 2／OB 6／OC 9）。本檔唔重覆。

## 依賴

```
上游：M4（γ_n 譜）· M13（gen-2^k 權重）· M14（u/d 裁決）
下游：M12 S4 · M14 S5 · M20（若 T3 需新載體）
```

## 註冊

`register.py` pid = `M19`｜state = `registered`｜verdict = `-`｜source = 調律者 2026-11-14「先開這幾個 M，但不在本 session 推」
