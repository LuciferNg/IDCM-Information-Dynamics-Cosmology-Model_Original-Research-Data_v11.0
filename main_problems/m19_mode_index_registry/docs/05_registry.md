# docs/05 — 登記冊（R4 / P4）：三態 ＋ frame 欄 ＋ 可執行引用閘

**日期**: 2026-11-14
**腳本**: `scripts/05_registry_build.py`｜**log**: `data/run_05_registry_build.log`（5/5 PASS）
**機器可讀**: `data/05_registry.json`

---

## 一、登記冊（正典 frame = F1，γ 升序）

| 粒子 | `n` | `γ_n` | 狀態 | 支撐 | 值之來源 |
|:--|--:|--:|:--|:--|:--|
| e | 1 | 0.005095771 | **STRUCT_CLASS** | {1,5,9,13} | T3 槽（類別）＋ 質量 rank（值） |
| u | 2 | 0.020357116 | INPUT_RANK | 全 9 模 | 質量 rank |
| d | 3 | 0.045706268 | INPUT_RANK | 全 9 模 | 質量 rank |
| s | 4 | 0.081014053 | INPUT_RANK | 全 9 模 | 質量 rank（4/5/7 結構皆可容） |
| μ | 5 | 0.126100550 | **STRUCT_CLASS** | {1,5,9,13} | 同上 |
| c | 7 | 0.244642021 | INPUT_RANK | 全 9 模 | 質量 rank |
| τ | 9 | 0.398917518 | **STRUCT_CLASS** | {1,5,9,13} | 同上 |
| b | 11 | 0.585786438 | INPUT_RANK | 全 9 模 | 質量 rank |
| t | 13 | 0.801444667 | INPUT_RANK | 全 9 模 | 質量 rank（13 非結構強制，見 `docs/04` §二 D6） |

**三態定義（可計算，非人手標）**

```
STRUCT_CLASS  容許支撐 ⊊ 9 模式集  ⟹ 值層受結構約束至「類別」
INPUT_RANK    容許支撐 = 全 9 模式 ⟹ 值層零結構約束，值由質量 rank 給
UNREGISTERED  冇任何 frame 可指（本輪 0 例）
```

## 二、★ 引用閘（F4，可執行；1 PASS / 2 BLOCK）

```
REQUIRED = (frame, state, provenance)

合法  {frame:"F1", state:"INPUT_RANK", provenance:"質量 rank"}   → PASS
缺欄  {state:…, provenance:…}                                    → BLOCK（缺 frame）
裸用  {k_q: 7}                                                   → BLOCK（缺 frame, state, provenance）
```

**規約（寫入本 M 及下游引用）**

```
R1  禁裸寫 "k_q = 7"；必寫 "F1 之 c = 7（state = INPUT_RANK）"
R2  任何含 c／b／t 之數值評估，須同時標 frame 與 state
R3  state = INPUT_RANK 之值，唔可以作為零參數預言之依據（M12 L201 禁區）
R4  state = STRUCT_CLASS 之值，只可用於「類別級」斷言（如 n ≡ 1 mod 4），
    唔可以宣稱「e 必須落 n = 1」
```

## 三、frame 欄（F3）：三 frame 之合法用法

| frame | 合法用途 | 禁止用途 |
|:--|:--|:--|
| **F1** γ 升序 | 正典登記；類內排序 | 宣稱值層有結構出處 |
| **F2** \|Δm\| 降序 | 對照／敏感度（已被 F1 之 S1 排除，`docs/03` B4b） | 當獨立理論 |
| **F3** CSSG-v3 | Layer B 支撐編碼（`m = M₀Σσγ`） | 當索引表（非雙射，見 `docs/02` §三） |

## 四、用法示例（三句可照抄）

```
① 查「c 喺邊個模式」
   → F1：c = 7｜state = INPUT_RANK｜provenance = 質量 rank
   → 唔可以寫「結構給出 c = 7」

② 查「輕子有冇槽」
   → STRUCT_CLASS：{e,μ,τ} ⊂ {1,5,9,13}（n ≡ 1 mod 4）
   → 只可講類別，唔可講「e = 1」

③ 要評 C(s,c) 數值
   → 先聲明 frame = F1；再標 state(s) = INPUT_RANK、state(c) = INPUT_RANK
   → 任何結論須寫「框架條件下」，不得宣稱零參數
```

## 五、唔可用嘅量（防誤用）

```
✗ 容許矩陣之權重比（720 vs 2880）當機率 —— 冇物理測度；純計數效應
✗ 「t 落 13 之剩餘 ≡1 槽」當結構 —— t 支撐 = 全 9 模（docs/04 §二 D6）
✗ 以 4.392 bit 宣稱「k_q 三分之二係結構」—— bit 數係全域，唔可逐粒子攤分
```

## 六、判

| 項 | 判 |
|:--|:--|
| 登記冊覆蓋 9 粒子 | ✅ 9/9 |
| 三態由判據決定 | ✅（3 結構類別 ＋ 6 輸入） |
| frame 欄強制 | ✅（禁裸 k_q） |
| 引用閘可執行 | ✅ 1 PASS / 2 BLOCK（真攔） |
| 跨記錄不一致 3 條已列（GX1–GX3） | ✅ 見 `docs/06` |
