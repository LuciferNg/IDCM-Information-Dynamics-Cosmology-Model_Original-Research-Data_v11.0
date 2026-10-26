# M4 docs/32 — 「K_tree 味分支結構」偵察：**前提與已註冊結構衝突**，味屬 K_chain 而非 K_tree

> 方法：先查已註冊結構（RAG + v9 原文），再決定可否建設；**不得自行發明新結構**。
> 依據：`IDCM-v9.0--W-Field-Lattice/wfield_lattice/03_extended_derivation.md` §D.1、
> `.../paper/17_paper.md` L217–231、`.../fermion/11_dirac_lattice_derivation.md` L181、`expansion_plan.md` L43；
> 色荷歸屬見 O216（色 = 嵌入層標籤，非譜量）。
> 狀態：🔴 **任務 (a) 之前提不成立**；🟡 真缺口改述並發還（非 M4 可推）｜2026-10-25

## 一、任務 (a) 之前提

`docs/31` §十 認為 A2 之推導斷裂點可藉「`K_tree` 之味分支結構（`s` 與 `d` 佔不同分支類型）」修補。
本輪先查：**此結構是否已註冊存在？**

## 二、偵察結果（v9 原文，逐條引）

| 來源 | 內容 |
|:--|:--|
| `03_extended_derivation.md` §D.1 | 「Generation modes = chain eigenmodes `k ∈ {1, 4, 5}` with (s,d) grading：k=1 → Up sector；k=4 → Down sector；k=5 → Lepton sector」 |
| `paper/17_paper.md` L217–231 | 「Chain eigenmodes 按 Z₂ grading `d = {3, 2, 0}` 分級 ⟹ 三代，無 fine-tuning」 |
| `fermion/11_dirac_lattice_derivation.md` L181 | 「Z₂ grading `d = {3, 2, 0}` → 3 physical chiral modes」 |
| `expansion_plan.md` L43 | 「K_chain 的 Z₂ grading 自然避開 Nielsen-Ninomiya doubling；88 nodes 中只有一半 carry physical chiral modes」 |
| O216（已註冊）| 色荷 = **嵌入層標籤**、非譜量；色數 `3 = B − 1 = N_gen` |

**⟹ 味（flavour）之標籤住在 `K_chain`（模指標 `k` ＋ Z₂ grading `d`），唔係 `K_tree` 嘅分支。**
色荷則屬**嵌入層**標籤（O216）。故「`K_tree` 味分支結構」**在已註冊結構中不存在**；
`d` 亦唔係樹深度（係 Z₂ grading，取值 {3,2,0}），`docs/31` §十 之用語需更正。

## 三、此偵察對 `⟨O⟩` 問題之直接後果

若味 = 鏈模標籤，則 isovector 算符 `ūu − d̄d` 在**模空間為對角**：

```
ūu − d̄d  ⟶  P_{k=1} − P_{k=4}      （投影子之差）
⟹ ⟨h| ūu − d̄d |h⟩ = w_{k=1}(h) − w_{k=4}(h)
```

即 `⟨O⟩` 係**態之模佔有率差**，係「態組成」量，**唔係樹幾何量**。
此與 `docs/31` 之密度插入（樹幾何）被否證**互相印證**：機制唔喺樹上，喺鏈上。

**但**：`s` 與 `d` 同屬 down sector（`k = 4`），故「`s` 與多餘 `d` 之區別」**唔可能**來自 sector 標籤；
若要區別，只能來自同 sector 內之 grading／世代細分 —— 而該細分在 v9 條目中**未給出**（只列一 mode 對一 sector）。

## 四、結論報告

1. 🔴 **任務 (a) 之前提不成立**：無「`K_tree` 味分支結構」可建；味屬 `K_chain`（`k`、Z₂ `d`），色屬嵌入層（O216）。
   本輪**不**發明新結構（違反「結構必然」紀律）。
2. 🟡 **真缺口改述（第三次）**：`⟨O⟩_ch = w_{k=1} − w_{k=4}` 需要**強子在 `K_chain ⊗ K_tree` 基下之態**，
   而 M4／M5 現有對象只到「兩體徑向骨架波函數」（M5 41×41）與「味 sector 標籤」（v9），**兩者之間未有接駁**。
   ⟹ 缺嘅唔係「味分支規則」，而係「**強子之多夸克態在模基下之展開**」。
3. 🔴 **A2 不可由 M4 現有結構推導**（既無樹味分支、亦無模基強子態）；維持後驗。
4. ⚠️ **`docs/31` §十 用語更正**：`d` 係 Z₂ grading（{3,2,0}），唔係樹深度；「味分支」概念撤回。

## 五、發還（register；屬結構建設 O，非 M4）

| 項 | 內容 |
|:--|:--|
| 缺口 | 強子之多夸克態在 `K_chain ⊗ K_tree` 模基下之展開（`w_k(h)`）|
| 為何非 M4 | M4 之範圍係「質量譜讀出與 8 通道審計」，非「建立強子態基」；且此基一旦建立，`⟨O⟩` 只需一行投影，屬上游 |
| 驗收（若他 O 建立）| 須零參數給 `w_1 − w_4` 之三錨（0.891／1.498／2.293）或其等價，並解釋 `n_s` 依賴 |
| M4 保留 | `em(K⁰−K⁺) ∈ [−0.52, +0.55] MeV` 之型式無關預言；待外部 lattice QED 落刀 |

## 六、M4 現況（該線）

| 項 | 狀態 |
|:--|:--|
| 8 通道 em 歸因 | 🔴 已判死（`docs/18`；`docs/26` 窗撤銷後仍無可行機制）|
| `⟨O⟩ ≡ 1` | 🔴 否證（p = 2.0×10⁻⁵）|
| `⟨O⟩` 形式 | 🟡 不可辨識（7/12）；A2 後驗 |
| 密度插入 | 🔴 禁區（`docs/31`）|
| 樹味分支 | 🔴 前提不存在（本輪）|
| `em(K) ≈ 0` | ✅ 型式無關預言，待外部 |
| 模基強子態 | 🔴 未建，已發還 |

## 七、與其他檔案關係

- ← `docs/31`（斷裂點）、`docs/30`（包絡）、`docs/29`（形式）、v9 `wfield_lattice/03_extended_derivation.md` §D.1 / `paper/17_paper.md` / O216
- → `01_problem_record.md` §二十、register（note ＋ 發還）
- 更正：`docs/31` §十 之「味分支結構」用語（見本文 §四.4）
## 八、發還已正式執行（2026-10-26，append-only 附註）

§五 之發還已由 register 落地：

- `register.py handoff --pid M4 --target-pid M7`（target `M7` 已註冊；交接原文見 register event）
- 承接單位：`IDCM-v11.0--Reality-Architecture/main_problems/m7_hadron_model_basis/`（pid `M7`，2026-10-26 開檔，**只開不推**，推導留新 session）
- ⚠️ 命名註記：調律者原話「開新m6」，但 **M6 = `m6_scale_dof_theorem`（尺度自由度定理）已存在** ⟹ 本題實際開為 **M7**
- M4 側維持不變：A2 仍後驗；`em(K) ∈ [−0.52, +0.55] MeV` 型式無關預言仍待外部 lattice QED
