# M12 / 全鏈審計（audit-verify）：R1–R3 全部產物之機械、獨立、文檔三層核

**日期**：2026-10-29｜**腳本**：`scripts/15_full_audit.py`、`scripts/l2_independent_verify.py`、`scripts/l3_doc_crosscheck.py`
**判詞**：**全鏈可重現、數字可回溯、無後驗湊數**；機械層 `FAIL 0`／`WARN 1`，L2 `14/14 PASS`（76 量），L3 `55/55 PASS`＋溯源 `0` 缺口。
**審計三件 sha**：`15_full_audit.py dc3a75fae0fc5b36`／`l2_independent_verify.py 901bebc1a8494b03`／`l3_doc_crosscheck.py ca8ee18b53be11b7`；log 行數 `105 / 124 / 38`。
**本節性質**：審計（不出新物理結論、不出預言）；所有被審數字皆為 R1–R3 既有者。

---

## 一、範圍與方法

| 層 | 做什麼 | 產物 |
|:--|:--|:--|
| **A（機械）** | 檔案對稱、可重現性、產物綁定、路徑衛生、文檔衛生、README 覆蓋 | `data/15_full_audit.json`、`data/run_15_full_audit.log` |
| **L2（獨立重推）** | 由**原始輸入**獨立重算全部承載數字；不 import／call 任何 L1 腳本，不讀 `docs/` | `data/l2_independent_verify.json`、`data/run_l2.log` |
| **L3（文檔交叉核）** | 每條文檔宣稱 ↔ L2 key；已否證值須保留且有更正節；文檔數值須溯源至機器產物 | `data/l3_doc_crosscheck.json`、`data/run_l3.log` |

L2 之**獨立性刻意設計**（與 L1 走不同算法）：

- κ 譜：L1 用 `sympy.charpoly`；L2 用 `rank_ℚ + tr K^(1) + tr((K^(1))²)` 三量推論非零特徵值必為 `±√(tr K²/2) = ±√5`（**不需要** charpoly）。
- RH 零模：L1 用 `Fraction` 全消去；L2 用多模數 `GF(1000003)` 消去得 rank 54。
- γ 恆等式：L2 同式以 `4sin²` 與 `2−2cos` 兩路互核。
- `⟨O⟩` 之 σ：L2 由 `σ_em`、`σ_X` 與質量差**重新推導絕對傳播**，不取用 M11 之 σ 欄。

---

## 二、A1／A2：檔案清單與對稱

- `docs/`：**12** 份（`01`–`11` ＋本審計 `12`）
- `scripts/`：17 支 = L1 `01`–`14` ＋ `15_full_audit.py` ＋ `l2_independent_verify.py` ＋ `l3_doc_crosscheck.py`
- `data/`：17 json ＋ 17 log（全部 canonical）
- NN 編號 **1..14 連續無缺號** ✅
- 根目錄：`01_problem_record.md`、`SESSION_START.md`、`README.md` ＋ 各自 `.pre_patch_*` 備份（append-only 紀律之實證）

## 三、A3：可重現性（14/14 支 ×2 次）

| 腳本 | log 行數 | sha256（前 16） | 兩次跑 byte-identical |
|:--|--:|:--|:--:|
| `01_r1_domain_additivity.py` | 122 | `f0750a3c073c040b` | ✅ |
| `02_reverse_gate.py` | 90 | `0eda71d1b79ae672` | ✅ |
| `03_paircount_gate.py` | 37 | `8fde8027bc622ca0` | ✅ |
| `04_boundary_vertex.py` | 59 | `e6834449d490ecc9` | ✅ |
| `05_fault_tolerance.py` | 29 | `8877a566521fada5` | ✅ |
| `06_fault_tolerance_bits.py` | 16 | `10d5d9e1cf845de6` | ✅ |
| `07_two_body_exhaustion.py` | 44 | `640d53e9417ae725` | ✅ |
| `08_resolution_density.py` | 94 | `5308e6c334fa8d8e` | ✅ |
| `09_readout_sigma_budget.py` | 59 | `3a61bcc77a932292` | ✅ |
| `10_pairtype_counting.py` | **69** | `6763f0e864b01f83` | ✅ |
| `11_splitting_ratios.py` | 87 | `b69eaa1ee06ef13a` | ✅ |
| `12_kappa_contraction.py` | 78 | `db23ba8549d05c33` | ✅ |
| `13_frustration_operator.py` | 77 | `25cbaf1d3b331082` | ✅ |
| `14_canonical_incidence.py` | 50 | `a3b3d8bb2daf86bd` | ✅ |

**★ 方法修正（審計自身之漏洞，F8）**：L1 腳本會**自行覆寫** `data/run_*.log`，故「跑完再與 canonical 比對」係**自證**。本輪改為：跑之前先快照 canonical bytes，再比 (i) 每次 stdout ↔ 快照、(ii) 跑後 canonical ↔ 快照。改後 14/14 仍全綠，且兩次跑彼此一致 ⟹ 結果非自證。

## 四、A4／A5：產物綁定與路徑衛生

- 每個 `data/NN_*.json` 皆有產者腳本；腳本宣告之輸出全部存在；**0 孤兒** ✅
- `scripts/*.py` 內**無** `/tmp/` 寫入 ✅（審計自身以 shell redirect 寫 `/tmp` 作比對緩衝，屬審計暫存，不進記錄）

## 五、A6／A7／A8：文檔衛生

| 項 | 結果 |
|:--|:--|
| `01_problem_record.md`／`SESSION_START.md`／`README.md` 狀態行 | 三者一致（`partial`）✅ |
| `final_result.md` | 未建 ⚪ 接受：本記錄**未 close**（依 O238 教訓，內部窮盡 ≠ 有結論；close 觸發器 = 新數據或新機制）|
| 11 份 docs 之 ``` 圍欄偶數、`U+FFFD` = 0 | 11/11 ✅ |
| README 檔案表覆蓋全部 docs | 修後 ✅（原僅列 R1，見 F3）|

## 六、L2 獨立重推：14 檢查全 PASS、76 量

五條恆等式（兩路互核）與九條交叉核：

| key 群 | 獨立重推結果 | 對應文檔值 |
|:--|:--|:--|
| `sum_gamma_44` | `90.0000000000` | `Σγ=90` |
| `O_np / O_K / O_S / O_X` | `1.1262875 / 2.4213281 / 1.9089844 / 2.171875` | `1.1263/2.4213/1.9090/2.1719` |
| `z_np / z_S / z_X / z_K` | `0.6835 / 5.7722 / 11.2923 / 14.1183` | `0.68/5.77/11.29/14.11` |
| `sigO_np / _S / _X / _K` | `0.184766 / 0.157476 / 0.103777 / 0.100673` | `0.1848/0.1575/0.1038/0.1007` |
| `z_decay` | `−3.2240` | `−3.223` |
| `seg_rel` | `0.6836`（段寬/σ）| `0.68σ` |
| `r_K / r_S / r_X / r_SS` | `−3.041756 / 3.716756 / 5.296397 / 6.246656` | 同值 |
| `sig_r_K … sig_r_SS` | `0.015940 / 0.029705 / 0.163837 / 0.058885` | 同值 |
| `kappa_T3 / T1 / nz_full / K1_nonzero` | `−8 / −380 / 1341 / 10` | 同值 |
| `kappa_K1_tr / tr2 / rank / absEig` | `0 / 10 / 2 / 2.2360679775` | `±√5` |
| `kappa_support` | `1·3·21·34 = 2142` | `{1,3,21,34}` |
| `rh_V / E / frust / H / rank / nullity` | `60 / 120 / 32 / −56 / 54 / 6` | 同值 |
| `rh_hist 0/1/2` | `8 / 40 / 12` | `8/40/12` |
| `sum_matter_n / sum_gamma_9` | `55 / 2.30906440147` | 同值 |
| `slot_col 1/2/3` | `14 / 18 / 23` | `(14,18,23)` |
| `z_K1_slot2/3`、`z_K2_slot2/3` | `−81.84 / −22.30`、`−74.30 / −18.59` | `−81.8/−22.3`、`−74.3/−18.6` |
| `sigT_3ch` | `0.26402` | `0.2640` |
| `coeff_h` | `[+1.00001, −1, +1, +1]` | 見 F5 |
| `floorX_np / floorX_Ξ` | `0.041356 / 0.079749` | 見 F1 |

⟹ **四錨、四 z、四 σ、四分割比、κ 五量、RH 六量、O137 八量，全部由原始輸入獨立重現**；無一項需回頭讀 L1 腳本或文檔。

## 七、L3 文檔交叉核：55/55 PASS

- **55 條宣稱**逐條綁定至 L2 key（`num`／`last`／`multi`／`set` 四種型）：**PASS 55、FAIL 0、MISSING 0**
- **REFUTED 機制**：已否證之舊值（`docs/06` 之 `0.2076`）須**仍存在**（append-only）且該文檔須有更正節 ⟹ PASS
- **溯源掃描**：文檔內全部 ≥3 位小數之數值（共 631 個機器數值池），須能在 `data/*.json`／`data/*.log`／`scripts/*.py` 找到文檔精度下四捨五入相等者 ⟹ **未溯源 = 0**
- **已分類接受 5 個**（披露而非隱藏）：`docs/02` 之 `−2.0695`、`4.9538`（跨記錄引 M2 `docs/02 §三`）、`4.2257`（★ R1 inline 枚舉結果，未落地腳本 ⟹ 僅旁證）、`docs/05` 之 `0.0026`、`0.026`（文內算術）

---

## 八、審計發現（逐條）

### F1 🔴 實質錯誤：`σ_O(Ξ)` 之傳播式

- 本記錄 `docs/05 §二`／`docs/06 §二` 曾報 `σ_O(Ξ)` 傳播應得 `0.2076`，並據此向 M11 發「`σ_E=0.2318` 疑 typo／不自洽」之 handoff。
- **M11 §14.6（23:10）判為導式錯**：把 `σ_E` 當相對誤差 ⟹ 項 = `|⟨O⟩|·σ_E/(cX) = 0.19666`，較正確項 `σ_Δm/(cX)` 多乘一個 `⟨O⟩ = 2.172`。
- **本輪 L2 獨立確認**：`sigO_X = 0.1037768`（檢查 C7）⟹ `0.1038` 係正典，`0.2076` 錯。
- 處置：`docs/05 §六`、`docs/06 §七` append 更正（原文不動）；handoff 結案。
- **連帶結論換位**：route (A) 之瓶頸**不是** em 錨，而是 `σ_X = 0.094` 之硬地板（`n−p 0.0414`、`K 0.0889`、`Σ 0.0701`、`Ξ 0.0798`，全 > 目標 `0.017`）⟹ 須 `σ_em ≤ 0.0435` **且** `σ_X ≤ 0.0386`（2.4×）。

### F2 🟡 `log` 與 `stdout` 不一致（5 支腳本）

`scripts/02`–`06` 之 `say("[written] …")` 寫在**寫檔之後** ⟹ canonical log 缺末兩行，與 stdout 不一致（違「data 與 stdout 一致」紀律）。
處置：五支腳本改為先 `say` 後寫檔（各自 `.pre_patch_*` 備份），重生成 canonical；實測新 canonical 與**修前完整 stdout 逐 byte 相同** ⟹ 數字零變動。

### F3 🟡 README 檔案表過期

原只列 R1 兩項（`docs/01`、`docs/02`）。處置：重寫（備份後），現覆蓋全部 docs／scripts／data。

### F4 🟡 行數記錯

`docs/07`／`SESSION_START` 曾稱 `scripts/10` 之 log = 68 行，實測 **69** 行。處置：`docs/07` 首行附註更正（原文句保留）。

### F5 🟡 `⟨O⟩` 含通道相依符號約定（新發現）

讀出式對 `K⁺−K⁰` 取 `Δm − E_ch = −2.5600` ⟹ `coeff_K = −1`（其餘三通道 `+1`）。
故「四錨全部 > 1」係**符號約定**之結果，`⟨O⟩` 實為**帶約定之絕對比值**，非有號荷；有號資訊只存於純分割比 `r`（`r_K = −3.0418`）。
處置：`docs/06 §七.4` 記錄；`scripts/l2..` C5/C6 逐位確認。**此點使 R1-Q2「加性失效」之敘述須以 `r` 而非 `⟨O⟩` 之號為準**——R1-Q2 之四 z 由 |⟨O⟩| 與 σ 算得，結論（加性不成立）不受影響。

### F6 🟡 `docs/02` 之 `4.2257 MeV` 不可重現

「全枚舉 34³ 組最小 max-dev 4.2257」係 R1 期間之 inline 運算，未落腳本。處置：於此列為**旁證、不作判據**；若日後需引用，須重寫為腳本並落地 log。

### F7 ⚪ `final_result.md` 未建

記錄未 close（`partial`），依記錄紀律於 close-out 時以 `RESULT_TEMPLATE.md` 建置。非缺陷。

### F8 ⚪ 審計自身漏洞（已修）

原 A3 以「跑後比對 canonical」判定可重現，但腳本會自行覆寫 canonical ⟹ 自證。已改為跑前快照法（見 §三）。

---

## 九、判詞

| 判準 | 結果 |
|:--|:--|
| 全鏈可否重現 | ✅ 14/14 ×2，跑前快照下 byte-identical |
| 數字可否獨立重推 | ✅ L2 由原始輸入重推 76 量，14 檢查全 PASS |
| 文檔與機器是否一致 | ✅ L3 55/55 ＋ 溯源 0 缺口 |
| 是否含後驗湊數 | ✅ 無：全部候選形式已於事先判準下判死；`D9`／`4.2257` 等後驗量已明文標為旁證 |
| 未決 | `F6`（1 個 inline 數值不可重現）、`F7`（close-out 時建 final_result） |

**⟹ M12 之 R1–R3 全部產物通過三層審計；唯一實質錯誤（F1）已由 M11 與本輪 L2 兩路確認並附註更正；殘餘為 1 個旁證數值之可重現性與 1 份待 close-out 文件。**

---

## 十、補記（append 2026-10-29；上列原文不動）

**F2b 🟡 審計腳本自身犯 F2 同一缺陷。**

`scripts/15_full_audit.py` 原於寫檔**之後**才 `say("[written] …")` ⟹ canonical `data/run_15_full_audit.log` 缺末兩行（105 行），與 stdout（107 行）不一致 —— 即本審計用以判 `docs/02`–`06` 之同一條紀律，審計工具自己違反。

處置：`say` 移至寫檔之前（`.pre_patch_*` 備份），重跑確認 **stdout ↔ canonical 逐 byte 相同**、107 行、`FAIL 0 / WARN 1`。

**故本文 §首之 sha／行數以本節為準**：

| 項 | 原文（§首） | 更正後 |
|:--|:--|:--|
| `scripts/15_full_audit.py` | `dc3a75fae0fc5b36` | **`89e955e2a7110eef`** |
| `data/run_15_full_audit.log` | 105 行 | **107 行**（`4d2e9770c1b3ba71`）|
| `run_l2.log` / `run_l3.log` | 124 / 38 行 | 不變 |

教訓（入 skill）：**凡「腳本自行寫 canonical log」者，`[written]` 必在寫檔之前 say**；否則 log 與 stdout 必然差兩行，而此差異會被任何以 log 為快照之審計判為內容漂移。

---

## 十一、補記二（append 2026-10-29：結果整理輪）

1. **F7 解除**：`final_result.md` 已建（依 M11 `final_result.md` 之精簡格式：狀態／問題／已確立表／已否證／殘餘規格／驗證水位／交接），標明 `partial`。
   ⟹ A6 由 `WARN` 轉 **PASS**；機械層現為 **`FAIL 0 / WARN 0`**（log 107 行，stdout ↔ canonical 逐 byte 相同）。
2. **L3 擴範圍**：溯源掃描與 CLAIMS 由 `docs/*.md` 擴至**根目錄四件套**（`final_result.md`／`README.md`／`01_problem_record.md`／`SESSION_START.md`）。
   結果：`CLAIMS 55/55 PASS`、未溯源 **0**、已分類接受由 15 增至 **19**（新增者全為轉引 docs 或 M11 值）。
3. 故本文 §五 A6 行、§九「未決」欄之記載以本節為準：**A6 PASS、未決只剩 F6（`4.2257` 之 inline 旁證）**。
---

## 十二、補記三（R4 加入後之 A2 判準修正，F9）

**日期**: 2026-10-29
**觸發**: `scripts/16`、`scripts/17` 加入後，審計報 `[FAIL] A2 NN 編號連續（1..17）  |  缺 [15]`。

**病因**: A2 原以 `nn`（已排除 `SELF = 15_full_audit.py`）計連續性。1..14 時無感；加入 16/17 後 `max(nums)=17` 而 15 被排除 ⟹ 誤報缺號。

**修法**: NN 連續性屬「磁碟真相」，改由全部 NN 檔（**含 SELF**）計算：
`nums = sorted(int(s[:2]) for s in scripts if re.match(r"^\d\d_", s))`。
`SELF` 仍只在 **A3 重跑**時排除（避免遞歸自跑），語義分離。

**性質**: 判準瑕疵（審計假陽性），非數據瑕疵；無任何數值變動。
**修後**: `FAIL 0 / WARN 0`（docs 13／scripts 19／data 17 json＋19 log；A3 16/16 ×2 byte-identical）。

**附註（同輪，數字更正，原文不動）**：`scripts/16`／`17` 加入後 `A3` 由 14/14 升至 **16/16**，故 canonical 隨之增長：

| 檔 | §十 記載 | 現值 |
|:--|:--|:--|
| `scripts/15_full_audit.py` | `89e955e2a7110eef` | **`3115ea2a4bd5903e…`**（A2 判準修正） |
| `data/run_15_full_audit.log` | 107 行／`4d2e9770c1b3ba71` | **113 行／`7d83fbe52f48eb90…`** |
| `data/run_l3.log` | 38 行 | **42 行**（docs 13 份、根目錄 4 件） |

已驗：`stdout ↔ canonical` 逐 byte 相同（16／17 亦各自 `diff` 為空）；`A3 16/16 ×2 byte-identical`；L3 `FAIL 0 / MISSING 0 / 未溯源 0`。

**附註（R5／R6 加入後，數字更新）**：NN 腳本增至 `01`–`19` ⟹ `A3` 由 16/16 升至 **18/18**；canonical `data/run_15_full_audit.log` 現為 **121 行**（sha `e92f7947113e4011…`）；`data/run_l3.log` 現為 **42 行**。三層水位：A `FAIL 0 / WARN 0`、L3 `55/55`＋未溯源 **0**（821 值池）。`docs/14`、`docs/15` 之數值全部可溯源至 `data/run_18_*`、`data/run_19_*`。（此行之首版因 shell heredoc 未引號、backtick 被當命令替換而損毀，已於同日就地修復；內容不變。）
