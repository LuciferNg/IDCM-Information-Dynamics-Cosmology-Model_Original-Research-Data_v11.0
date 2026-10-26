# M7 docs/17 — 全鏈審計（audit-verify）

> 腳本 `scripts/07_audit_full_chain.py`（**23/23 PASS**）｜log `data/run_07_audit_full_chain.log`｜JSON `data/07_audit_full_chain.json`
> 依據 skill `idcm-audit-verify` §2（三層驗證）／§3（推導鏈一致性）／§4（傳播）；觸發：調律者 2026-10-26「先進行 audit-verify M7 session 全鏈」
> 本輪**只讀**取證，唯一寫入 = 本 M 目錄內（README／`run_chain.sh`／§水位／本檔）

## 一、審計範圍與方法

| 層 | 檢查 | 方法（**唔用原腳本同一算法**） |
|:-:|:--|:--|
| A 結構 | A1 檔案對稱／A2 產物綁定／A3 路徑衛生／A4 driver／A5 版本檔歸屬 | 掃 README ↔ 磁碟；JSON ↔ producer；`/tmp`／scratch 正則；driver 存在性 |
| B 重推導（L2） | B1–B9 headline 數字 | γ 用複數模方（非 4sin²）；Σγ 用結構計數；p 用 mpmath 上不完全 Γ；b̂/σ_b 用 **Fraction 精確正規方程**；MC 用解析常態尾互核；G3 導數用 **sympy 符號微分**；諧振用 **sympy eigenvals** |
| C 記錄 | C1 四件套＋register／C2 水位覆蓋／C3 append-only | 狀態行解析；register jsonl 最新 state；`git show --unified=0` 刪除行**逐字重現**判據 |

## 二、B 層：9 項獨立重推（全部 PASS）

| # | 量 | 本輪獨立值 | 記錄值 | 判 |
|:-:|:--|--:|--:|:--:|
| B1 | γ_n 兩法 max dev（n=0..44） | `1.07e-50` | — | ✅ |
| B1 | 配對恆等 `γ_n+γ_{44−n}=4` | `1.07e-50` | EXACT | ✅ |
| B2 | `Σγ`（n=0..44） | `90`（數值＝結構計數 22×4+2） | `90` | ✅ |
| B3 | `χ²(⟨O⟩≡1)` | `18.220726674` | `18.220726674` | ✅ |
| B3 | `p`(3 dof) | `3.96066e-4` | `3.9607e-4` | ✅ |
| B4 | `b̂` | `0.675330350` | `0.675330350` | ✅ |
| B4 | `σ_b` | `0.161347924` | `0.161347924` | ✅ |
| B5 | MC `p` vs 解析 | 比值 `1.0196` | — | ✅ |
| B6 | G3 導數 `∂ΔM/∂m_u /coeff` | `{N:1, Σ:1, Ξ:1}` | 讀法 A/B 同值 | ✅ |
| B7 | 三體相關（pairwise 諧核） | `−18.350342%` | `−18.35%` | ✅ |
| B8 | 每模容量 | `6 > N_q=3` | `6>3` | ✅ |
| B9 | `γ₄M₀`／所需 | `−2.2603662`／`+3.31008` MeV | 同 | ✅ |

### 2.1 ★ 自審：本 audit 首跑有兩處**自身算法錯**（已修，記錄在此）

| # | 首跑錯值 | 成因 | 修正後 |
|:-:|:--|:--|:--|
| B5 | `p_an = 5.689e-5`（比值 2.04） | 把 `Q(1/2, z²/2) = erfc(z/√2) = 2Φ(−z)` **再乘 2** ⟹ 四尾 | 比值 **1.0196**（MC `2.790e-5` ↔ 解析 `2.845e-5`） |
| B7 | `−59.18%` | `sympy.eigenvals()` 返回 `{0:1, 3:2}`，迭代 **keys** 未展開重數 ⟹ 少一個 √3 | **−18.350342%**（與 `docs/15` C5 一致） |

⟹ 兩處皆屬**審計器**錯，非記錄錯；若不修，會把「記錄正確」誤報成 MISMATCH（假陽性）。判據不變：跨方法一致才算 PASS。

## 三、C 層：記錄完整性

| 檢查 | 結果 |
|:--|:--|
| C1 register | `state=closed`、`verdict=negative`（`2026-10-26T12:59:45+08:00`） |
| C1 registry block | `01_problem_record`／`SESSION_START` 皆 `closed/negative` ⟹ **無 BLOCK-DRIFT** |
| C1 四件套狀態行 | README `🔴 closed / negative`；final_result `🔴 判死`；problem_record／SESSION_START 之首段為**開檔快照**（依 §0.1 不改）⟹ 以 registry block 為機器對齊點，一致 |
| C2 水位 | close 後 **28 檔**變更（docs/13–17 四輪 ＋ scripts/03–07 ＋ data 全重跑）；水位表**本輪補登**（`01_problem_record` §檔案閉合水位） |
| C3 append-only | 掃 10 commit；僅 2 處「刪 1 行」＝**尾行換行正規化**（`cat >>` 撞無 EOL），刪除行**逐字重現**於新增 ⟹ **零內容消失** |

## 四、A 層：5 項 finding 及處置

| # | finding | 處置 |
|:-:|:--|:--|
| A1 | README 缺 `scripts/07`、`scripts/02`（v1） | ✅ 已補兩行（07＝全鏈審計；02＝v1 歷史） |
| A1 | `docs/NN` 簡寫（`docs/12` 等）被判 missing | ✅ 審計器加簡寫解析（本 M／跨 M 分類） |
| A2 | `data/02_r5_mode_content_nogo.v1.json` 無 producer | ✅ 查明係 **byte-identical 重複副本**（sha `9f4912eab5572912`，與 `.json` 同）⟹ 白名單＋本檔 §A5 記錄 |
| A3 | scripts 內 `/tmp`／scratch 命中 1 | ✅ 係審計器**自身** pattern 字面 ⟹ 排除自身後 0 |
| A4 | 缺 `scripts/run_chain.sh`（§3.1 ⑥） | ✅ 已建：預設 dry-run，`--run` 逐 script 覆寫 log；鏈序 01→02→02b→03→04→05→06→07 |

### 4.1 A5 版本檔

`data/02_r5_mode_content_nogo.v1.json`／`.v1.log` 與現行輸出 **byte-identical**（`sha256 = 9f4912eab5572912…`／`c1576c42f9e3368e…`）。
屬早期命名遺留之重複副本。**建議清理**（刪 `.v1.*` 兩檔可減冗餘），但依政策需調律者逐項批准 ⟹ 現只記錄，不動手。

## 五、邏輯鏈檢查（§3.3）

| 項 | 判 |
|:--|:--|
| headline 數字有 producer | ✅ 全部：`docs/11`←03、`docs/12`←02b、`docs/15`←05、`docs/16`←06、本檔←07 |
| producer 用已否證之輸入 | ✅ 無（05 之 capacity 用 O216、06 用 M4 取向約定，皆在效） |
| script 硬編碼期望值 | 🟡 **1 處已記錄**：`02b` K3 之「sector ⟹ ≡2」row 用 `2*np.ones(3)` 硬寫；已於 `docs/12` §十 撤回、`docs/16` G3 用 sympy 重算證其為係數錯 |
| 同一量兩 script 不一致 | ✅ 無（18.35%／b̂／χ² 三處皆經 07 跨方法核對一致） |
| doc 引用數在 data/ 找不到 | ✅ 無（`04` 111 claim 全綠） |

## 六、總判

```
A 5／5（3 項即場修，2 項記錄待批）    B 9／9 獨立重推一致    C 3／3
合計 23/23 PASS；無阻塞級 finding；唯一待批＝刪 .v1 重複副本
```

**M7 全鏈可重現、可追溯、append-only 無破損。** 審計器自身兩處算法錯已修正並留痕（§2.1）。

## 七、機讀摘要


| 量 | 值 |
|:--|--:|
| 審計 npass | 23 |
| 審計 ntotal | 23 |
| 審計 b̂ | 0.675330350 |
| 審計 σ_b | 0.161347924 |
| 審計 χ² | 18.220727 |
| 審計 p | 0.000396066 |
| 審計 Σγ | 90 |
| 審計 B7 相關% | -18.350342 |
| 審計 MC/解析比值 | 1.0196 |
| 審計 C2 未覆蓋 | 0 |

## 八、C3 判據演進（append-only，交待兩次變更）

| 版 | 判據 | 觸發 |
|:-:|:--|:--|
| v1 | 刪除行是否**逐字重現**於新增（換行正規化可接受） | 首版 |
| v2 | ＋**水位 live-field 豁免**：`scripts/10_watermark.py` 每次重生 `### §檔案閉合水位` 全節（刪舊表 → 檔尾重貼），該節之刪除行不屬「內容消失」 | `64b65b7` 出現 43 行「未重現」，查係水位表重生（mtime 變 ⟹ 行文必變），非記錄改寫 |

豁免範圍嚴格限定：`^| N | \`path\`` 之表格行、節標題、`> live field` 註、表頭——其他任何刪除行仍須逐字重現。
現行判定：**23/23 PASS**（C3 之 `n_live` 欄另列，透明可審）。
