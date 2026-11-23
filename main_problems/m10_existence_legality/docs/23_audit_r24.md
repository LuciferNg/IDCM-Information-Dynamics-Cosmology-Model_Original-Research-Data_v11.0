# docs/23 — R24：全鏈審計（audit-verify）＋ R23 定義更正之傳播

> 觸發：調律者八點指令（腳本齊全／audit-verify／數字漂移／handoff 檢查處理／交付後重算／四件套同步／殘留清單）。
> 工具：`scripts/l1_reproducibility.py`、`scripts/l2_independent_verify.py`、`scripts/l3_doc_crosscheck.py`、
> `scripts/audit_verify.py`（本輪新建三層驗證閘）；報告 `data/l1|l2|l3_doc_crosscheck|audit_verify.json`。

---

## 一、四數（idcm-audit-verify §2.4）

| 層 | 結果 | 內容 |
|:--|:--|:--|
| **L1 原樣可重現** | **22/22**（R25 增 21、22 後） | 20 支腳本各跑兩次：rc 全 0、stdout 剝不穩定行後 sha256 相同 |
| **L2 獨立重推** | **10/10** | 自建 AME／NUBASE 解析 ＋ 獨立 11 通道實作，**不讀任何 data/*.json、不 import L1** |
| **L3 文檔交叉核** | **32/32** | 期望值只取 L2 `_derived`；細整數帶 anchor |
| **AUDIT 結構閘** | **15/15** | 檔案齊全／斷鏈 0／四件套＋水位／L1-L3 結果／06 交叉核／register／README 覆蓋 |
| L1 之 check 計數 | **199：187 PASS ／ 12 FAIL** | 12 全為**誠實負結果**（事前預期落空或樣本更正旗標） |

### 1.2 L2 首輪即捉到真錯

L2 首跑：`tier3 = 146`、`radius = 2.8355e−4`；L1 當時：`131`、`5.211e−4`。**兩者不一致 ⟹ 停線審計、定位根因**。

---

## 二、★ R23：`min_open_rel` 之通道集錯配（本輪最重要發現）

| 項 | 內容 |
|:--|:--|
| **位置** | `scripts/11_falsification_nubase.py`、`scripts/13_primordial_test.py` 之 `min_open_rel()` |
| **病徵** | docstring 寫「**11 通道**中最小之正 Q」，實作只計 **8**（β⁻／β⁺EC／α／p／n／2β⁻／2β⁺ECEC ＋ SF_opt） |
| **漏** | `SF_sym`、`cluster`、`ternary-α` |
| **後果** | 最小正 Q 偏高 ⟹ Tier 分層與失效半徑**系統性偏斜** |

### 2.1 更正後（11 通道）

| 量 | 舊（8 通道） | **新（11 通道）** |
|:--|--:|--:|
| 失效半徑 `rel*` | 5.21e−4 | **2.84e−4** |
| Tier III（S\P） | 131（72.0%） | **146（80.2%）** |
| Tier II 屏障帶 | 51（全 `SF_opt`） | **36（全 `cluster`）** |
| 屏障帶 A | [98,164] | **[110,164]** |
| `Nat\P` Tier 分佈 | 147＋51 | **162＋36** |
| G1 反例 | 8/51 | **7/36** |
| 單一同位素 11 條之通道／Tier | 多為 `SF_opt` | ¹⁰³Rh／¹²⁷I／¹³³Cs／¹⁴¹Pr／¹⁵⁹Tb → **cluster**；¹⁶⁹Tm／²⁰⁹Bi → **ternary-α** |

**不變**：`P = 90`、`A_max = 92`、`S\P = 182`、三機制 126／51／5、D5／D6 封印（定義不含數值）、
G1 仍不成立、G2 判別仍乾淨（gap 24.7 MeV）、單向性結論。

### 2.2 為何 D5／D6 不受污染

D5.5 定義為 `rel* := 由 data 之實測 max rel 給定`（**不含數值**）、D6.5 同 ⟹ 定義層穩健，
只有 data 值改變。此即 2026-10-28「定義不嵌計數」教訓之直接回報。

---

## 三、傳播（§4）

| 類 | 檔數 | 內容 |
|:--|--:|:--|
| 腳本更正 | 4 | `11`／`13`（min_open_rel 補 11 通道）、`12`（快照值＋D5.7）、`14`（樣本旗標＋E1c／動態分母） |
| 文檔數字 | **9** | `docs/10`、`12`、`13`、`14`、`16`、`17`、`20`、`README`、`final_result` |
| 補漏 | 5 | `docs/21`、`README`、`final_result`、`SESSION_START`、`01_problem_record` 之 `5.212e−4` |
| 交叉核 | — | `06_doc_crosscheck`：**368/368 claim** 全命中（前後皆 100%） |
| 更正節 | 9 | 每份受影響 doc append「R23 更正（2026-11-04）」表 |

**零未申報**：`06` 之 claim 全部對得上；無任何數字被靜默放寬容差。

---

## 四、腳本齊全性（指令第 1 點）

| 判 | 值 |
|:--|:--|
| 腳本 | **20**（`01`–`20`；R1–R22 每一步皆有對應腳本） |
| 文檔 | **23**（本檔） |
| 數據 JSON | 20 ＋ 2 refs |
| 每支腳本 → 對應 JSON | ✅（前綴配對；`06` 為交叉核輸出自身） |
| 每份 doc 見於 README | ✅（本輪補 `docs/22`／`docs/23`） |

---

## 五、斷鏈掃描（指令第 3 點）

- 引用掃描：docs 21／scripts 35／data 30 條 → **斷鏈 0**（唯一一條 M4 之 `27_em_sign_theorem.py`（在 `m4_mass_spectral_readout/scripts/`） 屬 M4，
  已改寫為明確跨記錄路徑 `../m4_mass_spectral_readout/…`）。
- 數值漂移：`06` 368 claim 全命中；L3 32 claim 全命中（獨立路徑）。

---

## 六、handoff（指令第 4–6 點）

- **入向 8 ＋ 他人回覆 note 2**：全部已裁決（R18／R20／R22），回執 8/8。
- **出向 19（5 作廢）**：live 14；待回 2 條（皆 M5：`em` 符號更正、σ 傳播）。
- **R23 之下游影響**：兩封對 M7 之 handoff（`9fa84abb`、`090e749a`）引「51 條／A∈[98,164]／8/51」⟹ **已發更正 handoff**（同版本 v11.0），並列明不變項。
- 超出 M10 範圍者（M5 之 `q²` 自能符號重檢、M7 之校準集數字）**只 register 交付**，未改他人檔案。

---

## 七、殘留與新開放問題（§8.1 六分類）

| 類 | 項 | 狀態 |
|:--|:--|:--|
| ① 已解（本輪新增） | R23 定義錯配 | ✅ 已修＋已傳播 |
| ② 本 O 可動殘留 | **W1/P7**（約定閘：v1 不封／v2 封 0.184σ；C 路抑制密度陰性；D γ 兩支皆不足） | 🟡 唯一殘留，槓桿在 M5 weight 約定 |
| ③ 結構上不可補 | `𝒞_B` 屏障能（R13）；`T2c` 對判準零影響（R21） | ✅ 已封 |
| ④ 外部接口（不計入本 M） | `O216` 色單態、`B−L`（O168）；M5 weight 約定；M4 之 F1 池定義對齊 | 🟡 只記依賴 |
| ⑤ 真未做（低價值） | 非 α 三分裂、β 延遲粒子（`docs/09` §八 判「對永久判詞無新增殺傷」） | ⬜ 未做，價值已評估 |
| ⑥ 新開放（本輪） | **①** `min_open_rel` 之 Tier 軸與 `Tier-M` 並存時，對外引用應以何者為準（本 M 判：**機制軸為主、裕度軸只作輔助**）；**②** `docs/07`／`scripts/08` 之 Tier 數字係 **6 通道階段值**，與 𝒞_v2 不可混引 | 🟡 需後續文檔標注 |

---

## 八、可重跑

```bash
/usr/bin/python3 scripts/l1_reproducibility.py --dir .
/usr/bin/python3 scripts/l2_independent_verify.py
/usr/bin/python3 scripts/l3_doc_crosscheck.py
/usr/bin/python3 scripts/audit_verify.py
```