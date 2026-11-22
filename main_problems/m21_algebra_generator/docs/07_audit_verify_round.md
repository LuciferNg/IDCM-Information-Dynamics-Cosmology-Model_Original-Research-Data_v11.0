# M21 docs/07 — 全鏈 audit-verify（第三輪）：腳本落碟、全鏈重跑、文檔交叉核、斷鏈

**腳本**: `scripts/09_full_audit.py`（log `data/run_09.log`）
**性質**: 獨立審計（**唔 import** scripts/00–08；由第一原理重算）
**日期**: 2026-11-20（第三輪）

---

## §0 段設計（跑前寫死）

```
A 檔案／編號／四件套／registry block／script docstring
B 全鏈重跑 vs 入庫 log（byte-identical；rc 記錄；02/04 rc=1 屬預期）
C 獨立重算 19 項 headline（γ/Σγ/CA1/普查/A₅/orbit/Λ²／Sym²/維數池/等變配對）
D 文檔逐 claim 交叉核（7 份 docs；命中率須 100%）
E 斷鏈（docs→scripts、文檔→docs、final_result→檔案、data JSON 可解析）
F 表格 pipe（含已知既有白名單 1 條）
G 自身水位自查
```

## §1 段 A：檔案完整性

```
scripts 00–09（10 支，連續）｜docs 00–07（8 份，連續）
data/ 每支 script 有 run_NN.log（09 之首跑前例外，A3b 已容許）
四件套齊：01_problem_record.md ／ SESSION_START.md ／ README.md ／ final_result.md
registry block state/verdict = partial（與 register row 一致）
全部 script 有 shebang ＋ docstring（A7）
```

**腳本落碟對帳（checkpoint 1 之答）**：推導鏈每一步都有 script ＋ log ＋ docs 段：

| script | docs 段 | log | JSON |
|:-:|:--|:-:|:-:|
| 00 重現閘 | `docs/00` 判準源 ＋ `01_problem_record §二` | ✅ | ✅ |
| 01 A₅ 普查 | `docs/01` | ✅ | ✅ |
| 02 CY₃ 池閘 | `docs/02` | ✅ | ✅ |
| 03 兩因子 | `docs/03` | ✅ | ✅ |
| 04 互導圖 | `docs/04` | ✅ | ✅ |
| 05 登記查重 | `docs/05` | ✅ | ✅ |
| 06 水位自查 | 本檔 §3 ＋ `docs/05 §5` | ✅ | ✅ |
| 07 M22 驗收 | `docs/06 §1` | ✅ | ✅ |
| 08 BR1 裁定 | `docs/06 §2–§3` | ✅ | ✅ |
| 09 全鏈審計 | 本檔 | ✅ | ✅ |

## §2 段 B：全鏈重跑

```
00 rc=0 byte-identical ｜ 01 rc=0 ｜ 02 rc=1（2 FAIL = 結論）且 byte-identical
03 rc=0 ｜ 04 rc=1（2 FAIL = 結論）且 byte-identical ｜ 07 rc=0 ｜ 08 rc=0
05 = **非決定性**（輸出含 corpus 掃描總檔數，隨他人新增檔案變動：9,576 → 9,703 → 9,707）
   ⟹ 改為「PASS/FAIL 模式穩定 ＋ 快照 log 自我刷新」，並在 06 內標 NONDET 白名單
06 rc=0（17/17）
```

**判**：一切決定性腳本 byte-identical；唯一非決定性者已標明機制與處理（非降級 —— 判準仍核 PASS/FAIL 模式）。

## §3 段 C：獨立重算（19/19 PASS）

```
C1   γ_n + γ_{44−n} = 4（n=1..43，Decimal prec 60）違例 0
C2   Σγ = 90.000000000000000（容差 1e−45）
C3   度集唯一：{1,3,3,4,5}（combinations_with_replacement 獨立枚舉）
C4   Σd(d−1) = 44 = 4·Σ(d−1)
C5   Σ(d−1)(d−4) = 0（加權均值 = B）
C6   Σ(d−1)² = 33 ｜ Σ(d²−1) = 55
C7–C9 |A₅| = 60；#{g²=e} = 16、#非對合 = 44；共軛類 {20,15,12,12,1}
C10  鏈側 gcd 型別 {20,10,5,5,2,1,1}
C11  反射 (21 對, 1 不動點{22}, 1 無伴{44})
C12  Σ d(d−1)/2 = 22 ｜ Σ d(d+1)/2 = 38 ｜ 22+38 = 60
C13  Σ_(i<j) d_i d_j = 98 = ((Σd)²−Σd²)/2
C14  群 orbit {20,12,12}／{2×22}；模 orbit {20,10,5,5,2,1,1}／{21×2,1,1}
C15  等變配對 = 0
C16  維數池 36 命中 = 0（h11 真缺口）
C17  {20,24} cut 密度 = 2/43
C18–C19 44／98 於池內各唯一命中（Σd(d−1)／Σ_(i<j)）
```

## §4 段 D：文檔交叉核

7 份 docs 之關鍵數字逐條核（`{1,3,3,4,5}`、`|A₅| = 60`、`16 = 1+15`、`44 ∤ 60`、`249`、`q_obs 3.33%`、
`6 ≥ 5`、`53.53%`、`4.808e−03`、`7.583e−6`、`40.25571591`、`1.643452676`、`F=7 / M=3`、`13 邊`、
`9,576`、`21 對`、`22 對`、`2/43`、`23 ＋ ∅`、`131 PASS`）—— 命中率 **100%**。

## §5 段 E：斷鏈

```
docs→scripts 引用：0 條斷鏈
文檔→docs 引用：0 條斷鏈（已排除跨 O 路徑 `m22.../docs/07_...` 之偽命中）
final_result 引用（scripts/docs/data）：全部存在
data/*.json：全部可解析
```

## §6 段 F：表格

```
全目錄表格塊不一致 = 1（`01_problem_record.md` §三 已知前置表，2026-11-14 原文，未 escape `|G|`）
⟹ 列 KNOWN_BAD 白名單，不修改原文（append-only）
```

## §7 本輪發現並修正之漂移（4 項）

| # | 漂移 | 來源 | 處置 |
|:-:|:--|:--|:--|
| 1 | `01_problem_record` 三處（L27／L46／L77）缺口因子 `5.222686e35` 有誤 | M8 → M21 register note（11-22 00:38） | §十三 附註更正（原文不動）；值 = `5.223686e35` |
| 2 | 水位 `87 PASS`（第一輪）仍留在 §10.6／`final_result`／`SESSION_START` | 本 session 第二輪擴鏈未同步 | 就地更正為 `131 PASS`（附 §十三 記錄） |
| 3 | `06` 水位 `13/13` 已擴至 `17/17` | 第二輪加入 07/08 | 同上 |
| 4 | `05` 輸出含 corpus 總檔數 ⟹ 非 byte-identical | corpus 增長（他人新增檔） | 06 加 NONDET 白名單 ＋ log 快照刷新；判準改核 PASS/FAIL 模式 |
| 5 | `docs/05`／`§10.3` 查重快照數字（9,576／23…）已與活值不符 | corpus 9,576 → 9,707 | 標「快照＋活值」雙軌；09 加 `D-global5`（快照標籤 ＋ 各關鍵詞族仍 ≥1 命中） |

## §8 判

```
A 全綠 ｜ B 決定性全 byte-identical（05 已標非決定性並有機制）｜ C 19/19
D 8/8 文檔命中率 100%（含 docs/07 自身）｜ E 斷鏈 0 ｜ F 僅 1 已知既有 ｜ G 自查
⟹ 本檔 scripts/09 = **59 checks，59/59 PASS**
⟹ M21 全鏈一致、無漂移（除 §7 已登記並修正之 5 項）
```

## §9 Register：M21 相關 handoff 全帳（第四輪核查）

**入向（2，全部已處理）**

| # | ts | 來源 | 內容 | 處置 |
|:-:|:--|:--|:--|:--|
| 1 | 11-20 11:12 | **M22** `handed_off` | BR1 典範映射缺口（請求裁定） | ✅ note（11-22 01:16）＋ `handed_off` 回執（01:18）⟹ 待承接 **0** |
| 2 | 11-22 00:38 | **M8** `note` | `5.222686e35` 位數手誤（三處） | ✅ ACK note ＋ `§十三` 附註更正（原文不動） |

**出向（6）**

```
11-20  M15  note      Q3 供給：K_tree = Σ(d−1)² = 33
11-20  M6   note      殘餘措辭（乘積 1 object ＋ N_IDCM）＋ dof 分支 (a)
11-22  M22  note      回執：BR1 對象層判死（三重）
11-22  M22  handed_off 回執（正式，清「待承接」）
11-22  M8   note      ACK：更正已落 §十三；下游影響 0
11-22  M15/M6 note    第二輪後更新：BR1 判死 ⟹ 識別為 count-level 約定；dof 帳本 ＋1 約定
```

**下游影響檢查（§9 之答）**

```
① M8 更正（5.222686 → 5.223686）：只係因子位數，缺口 35.717977 dex 不變 ⟹ **下游重算 = 0**
② BR1 判死：M22（B2 仍唯一出處，已回執）；M17 前已由 M22 docs/12 得「需 B2 者不可作檢定」
   ⟹ 無新增影響；M15（錨判定改寫）／M6（dof 帳本）已發更新 note
③ 超越 M21 範圍者：全部只經 register note／handoff，本 O 零寫入他 O 目錄
   （`find` 檢查：本輪只改 m21 目錄 ＋ problem_register）
```

**register 狀態**

```
v11.0/M21  partial / partial ｜ incoming 1 待承接 0 ｜ STATUS-DRIFT 0 ｜ BLOCK-DRIFT 0
REF-MISSING 18（全部屬 M4/O124/O140/O142/O148 等既有項，非本 O）
```
