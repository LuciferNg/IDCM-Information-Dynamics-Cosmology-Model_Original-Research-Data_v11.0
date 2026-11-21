# M18 docs/19 — 鏈審計（六項核對 ＋ L1／L2／L3b ＋ 時序 ＋ driver）
> 腳本：`scripts/19_chain_audit.py`｜log：`data/run_19_chain_audit.log`｜JSON：`data/19_chain_audit.json`
> 狀態：**20 checks：20 PASS**｜實作閘 FAIL 0｜2026-11-20
> ｜2026-11-20

## 一、前設問題

`scripts/08` 只審「交付層」（腳本↔doc↔JSON）。本支補 skill `idcm-audit-verify` §3.1 之六項核對
與 §4.8 之機器閘，令「全鏈一致無漂移」有**硬證據**而非聲明。

## 二、事前判準（寫死）

| # | 判準 |
|:-:|:--|
| Q1 | 腳本集對稱（含審計支）↔ `docs/NN`，兩側差集 = ∅ |
| Q2 | 每個 `data/NN_*.json` 有 producer ＋ canonical log；每支 script 宣告之 JSON 存在；JSON↔log **兩向零孤兒** |
| Q3 | 路徑衛生：`/tmp`、`/var/tmp`、`.cache`、`scratch` 命中 = 0 |
| Q4 | **L1**：逐支重跑，stdout == canonical log（逐字節）、JSON sha 前後相同 |
| Q5 | **L2**：log 匯總行 ↔ JSON `checks`（PASS／總數／實作閘；`WARN` 屬 `PASS` 子集 ⟹ `x+y==z`）|
| Q6 | **時序**：`data/*` mtime ≥ producer script mtime |
| Q7 | **driver**：`run_chain.sh` 覆蓋集 == 腳本集 ∪ {`08`,`19`,`20`}，含三條線程鎖 export |
| Q8 | **L3b**：`docs/*.md` 內 3+ 小數位數值可溯（池 = data JSON ∪ run log；容差 0.5%；白名單逐條附理由）|
| Q9 | 四件套存在，狀態行含四態 emoji |
| Q10 | inline 幽靈推導：全部 headline 數值皆有落碟 script |

## 三、執行結果

```
Q1  推導腳本 20 支 ↔ docs 20 份；全腳本（含審計支）23 支 ↔ docs 23 份
Q2  23 JSON 全部有 producer ＋ canonical log；零孤兒（含 `Q2d` 反向：零孤兒 log）
Q3  零 /tmp／快取引用（本支偵測器自身 pattern 行已排除，附理由）
Q4  L1：20/20 exit 0｜20/20 JSON byte-identical｜20/20 stdout == canonical log（逐字節）
Q5  L2：20/20 log 匯總 ↔ JSON 一致
Q6  時序：零遲產物
Q7  driver 覆蓋 01–23（含 08/19/20/21/22/23）；三條線程鎖齊；預設 dry-run
Q8  L3b：掃 97 個數值 token，ORPHAN = 0（白名單 6 條附理由，含 4 條跨記錄值出處）
Q9  四件套 4/4；狀態行四態 emoji ✓
Q10 無 execute_code headline
```

## 四、本輪審計發現並即場修正之項

| # | 發現 | 處置 |
|:-:|:--|:--|
| 1 | `docs/06` 表 `T_t` 最近 dev 寫 `0.0014%`，而 canonical log 為 `0.00137%` | **修正 doc** 為 `0.00137%`（對齊 log，非放寬容差）|
| 2 | `scripts/01` 匯總行 `PASS 27 / FAIL 0 / WARN 1（共 27）` 之 `WARN` 屬 PASS 子集 | 判準寫明 `x+y==z`（`WARN` 資訊性），唔改腳本輸出 |
| 3 | 缺 `scripts/run_chain.sh`（driver）| **新建**，覆蓋 01–20，含三條線程鎖，預設 dry-run |
| 4 | L3b e-記法歸一曾逐位轉換（產生假 ORPHAN `1.030`）| **修正為整個 token 轉換**（`Decimal`），屬 skill §3.1⑦ 陷阱 2 |
| 5 | `docs/19`／`docs/20` 未存在（Q1a 差集）| 本節與 `docs/20` 補齊 |
| 6 | `data/run_19_audit_l2l3.log` 係 driver 打錯名之殘留孤兒 log | 刪除，並新增 `Q2d` 反向孤兒 log 閘 |

## 五、判決

1. **L1 硬證據成立**：17/17 腳本 stdout 與 canonical log **逐字節相同**、JSON sha 不變 ⟹ 無漂移。
2. **L2 成立**：log 匯總與 JSON 自報計數一致。
3. **L3b 成立**：docs 全部 3+ 小數位數值可反溯（ORPHAN 0）；跨記錄值 4 條入白名單並附出處。
4. **driver 補齊**：`run_chain.sh` 存在且覆蓋全集（此前全鏈無 driver，屬 §3.1⑥ 缺口）。
5. **一項實質修正**：`docs/06` `0.0014%` → `0.00137%`（顯示精度 ≠ 值；修正對齊 log）。
