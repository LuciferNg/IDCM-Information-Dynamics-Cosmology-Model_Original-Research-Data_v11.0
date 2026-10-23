# M4 — 質量 = 譜讀數差值：全費米子譜差閉合

## 狀態: 部分閉合（2026-10-24）

質量 = 譜讀數差值：全費米子譜差閉合（差值閉合力學 v2，R7 定案承接）。
已確立：差值閉合力學（Δm²+R 直出譜，唔經 M_comp）＋ K_chain 譜代數（γ_n = 2−2T_n(c)，全譜 ∈ ℤ[γ₁]）。
攻堅路線：全費米子單點譜比（已全滅，density 閘）、lattice 精度重判（span≈12 升部分閉合）、A 譜身份（開放）。

詳見 `01_problem_record.md`／`SESSION_START.md`／`docs/13_script13_mdmu_second_instance.md`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 已確立結構 + R1–R4 結果（§八）＋ append-only 更正（§九）＋ 逐腳本文檔登記（§十） | 部分閉合 |
| `SESSION_START.md` | 起手式 + 接手起點 + 2026-10-24 執行記錄（§七） | 部分閉合 |
| `docs/` | 13 檔（`01–13`）：逐腳本記錄＝前設問題＋事前判準＋執行＋結論報告（原輪次報告已併入，見 §編號規則） | 部分閉合 |
| `scripts/` | 計算腳本（13 檔：01–08 遷入、09–13 新增） | 部分閉合 |
| `data/` | stdout log ＋ JSON（14 log ＋ 4 JSON） | 部分閉合 |

## 編號規則

- `docs/NN_scriptNN_<slug>.md`：**一對一**對應 `scripts/NN_<slug>.py`（NN = 01…13），內含「前設問題（腳本 docstring）＋ 事前判準 ＋ 執行記錄 ＋ 結論報告」
- 輪次／閘報告**已併入**對應腳本文檔，唔再另立檔（2026-10-24 調律者指令：只保留 01–13）：
  - 閘重建（8 支期望 vs 實跑）＋ 01 [C] bug 取證 ＋ 05 MC protocol 不符 → `docs/01–08` §六／`docs/01` §七／`docs/05` §七
  - density 統一 protocol ＋ 逐命中全表 → `docs/09` §六
  - γ₁↔γ₆ 代數身份 ＋ O89 階梯核證 ＋ log 等距 → `docs/10` §六–§九
  - FLAG lattice 數據來源、誤差傳播、1σ 帶、裁定表 → `docs/11` §六–§七
  - A 歸一化、族結構、M_comp、下一步 → `docs/12` §六–§七
  - T4b 第二實例全表 ＋ 全 session derived-vs-not／gap G1–G6 → `docs/13` §五–§六
- 所有腳本一律「假設從未執行過」重跑：canonical run log 存 `data/run_NN_<slug>.log`（附 sha256）

## 依賴

O257 (R7 docs/07 定案) + O66 (M_comp 降級) + O89 (γ₁ 冪階梯，2026-10-24 核證：M_GUT/PS/R 三階梯作廢) + O57 (m_ν 盲測) + P6 (§10 差值閉合) + O73 (c=Goldstone) + O212 (dof=1)

## 註冊

pid `M4`（created 2026-10-21；state partial 2026-10-24；handoff → P6 A 身份、O140 selection rule）
