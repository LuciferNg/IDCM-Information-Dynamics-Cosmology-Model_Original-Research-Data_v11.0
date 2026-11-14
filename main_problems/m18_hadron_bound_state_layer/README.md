# M18 — 強子束縛態層（hadron bound-state layer）

## 狀態

⚪ **registered — 未開推**（2026-11-14 開檔）｜本 session 零推導

## 一句話

束縛態譜 `E_h(n, J^PC)` **由邊個算子生成** —— 呢層全 IDCM 冇任何 M 承擔，而 M4 `Q2/Q3/Q8` 同 M5 質量差 em 閉合全部卡喺呢度。

## 為何要開

M4 `§11.5` 原文：

> ★真缺口：強子束縛態譜（**IDCM 現無任何 O 承擔**）
> ⟹ Q2/Q3/Q8 的真正上游係一個**未建立嘅結構層**，唔係任何現存 O 嘅待辦項。

M5 `final_result.md` L199：

> 剩餘唯一前線 = mode 佔用／J^PC 算子（blocked-by-O140）或非成對（state-dependent）項。

## 缺口證據（M5 三層掃描全敗）

| 掃描層 | 結果 |
|:--|:--|
| contact 算子標度維度 | 非缺口 |
| 支撐（support） | 不敏感 |
| 勢衰減律 ρ ∈ {1/4, 1/3, 3^−1/2} | 缺口 62–68× |

⟹ 缺口必然落喺**供給側自由度**（mode 佔用／`J^PC`／徑向激發）。

## 檔案

| 檔 | 內容 |
|:--|:--|
| `01_problem_record.md` | §〇 fence／Q1–Q4／T1–T4／P1–P5／R1–R5／禁區／4 Hard Tests |
| `SESSION_START.md` | 接手起點、重現閘、禁區、分工、交接 |
| `README.md` | 本檔 |
| `docs/` `scripts/` `data/` | 空（只開不推） |

## 依賴

```
上游：M4（譜代數、Q2/Q3/Q8）· M5（⟨1/r⟩、J^PC 阻塞）· M6（dof=1）· M7（模基退化，負結果）· M12（CY₃ 投影）
下游：M4 Q2/Q3/Q8 · M5 質量差 em 閉合
```

## 註冊

`register.py` pid = `M18`｜state = `registered`｜verdict = `-`｜source = 調律者 2026-11-14「先開這幾個 M，但不在本 session 推」
