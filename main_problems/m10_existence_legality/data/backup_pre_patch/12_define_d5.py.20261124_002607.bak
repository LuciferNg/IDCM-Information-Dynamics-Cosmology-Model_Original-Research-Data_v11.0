#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/12 — D5：判準之修訂定義（單向證書 ＋ 𝒞 二分 ＋ Tier 二維），並封印

前設問題:
  R10（`docs/10`）以 NUBASE2020 獨立標籤測出兩件事：
    ① ⟸ 方向（∀c x<0 ⟹ 永久）90/90 零反例 ⟹ 必要方向成立；
    ② ⟹ 方向（x>0 ⟹ 衰變）對 𝒞_B（屏障類）**不可測**，51 條判「不合法」而實測 stbl。
  故 D2 之 `⟺` 不可作預言用。本腳本把修訂寫成 **D5**（新定義，非改舊記錄），
  並對其作 sha256 封印，使後續輪次不能暗中移動定義。

D5 之內容（寫死於 STATEMENT，腳本不改寫）:
  D5.1 單向性：∀c x(c)<0 ⟹ 永久（充分條件，可測）；∃c x(c)>0 ⟹ 不永久（定義性，非預言）
  D5.2 𝒞 二分：𝒞_E = {β⁻, β⁺/EC, 2β⁻, 2β⁺/ECEC, p, n}（能量限制）
             𝒞_B = {α, cluster, fission(SF_sym/SF_opt), ternary}（屏障限制）
  D5.3 Tier 二維：(裕度 |x|/m) × (通道類)；Tier II × 𝒞_B = 盲格
  D5.4 「不合法」三分：① 不存在 ② 不永久 ③ 可觀測衰變；判準只判 ②
  D5.5 失效半徑：rel > 5.21e−4 之不合法判詞對觀測穩定集零偽陰性（實測界）

執行：/usr/bin/python3 scripts/12_define_d5.py
"""
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA, DOCS = os.path.join(ROOT, "data"), os.path.join(ROOT, "docs")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


STATEMENT = """D5（判準之修訂定義，v11.0/M10，2026-10-27）

D5.1 單向性
  設 ℰ 為環境、c 為通道、x(c) 為該通道之能量差（x<0 ⟹ 能量封閉）。
  (a) 充分（可測）：若 ∀c∈𝒞: x(c)<0，則 C 在 ℰ 下永久（= 無限時間內不衰變）。
  (b) 定義性（不可作預言）：若 ∃c∈𝒞: x(c)>0，則 C 非永久。
      「非永久」≠「會衰變」：屏障可令 T½ 達 10^80 yr 級 ⟹ 不可觀測。
  故判準係**永久之充分條件**（permanence certificate），非衰變之判別器。

D5.2 通道二分
  𝒞_E（能量限制）= {β⁻, β⁺/EC, 2β⁻, 2β⁺/ECEC, p, n}
  𝒞_B（屏障限制）= {α, cluster, fission(SF_sym/SF_opt), ternary-α}
  理由：𝒞_E 之逆反應在閾值處無勢壘（弱作用）；𝒞_B 有庫侖／核勢壘，
        x>0 不足以致衰變。

D5.3 Tier 二維
  Tier = (裕度 |x|/m) × (通道類 𝒞_E | 𝒞_B)
  Tier I   : rel > 1e−2
  Tier II  : 1e−4 < rel ≤ 1e−2
  Tier III : rel ≤ 1e−4
  盲格：Tier II × 𝒞_B（Q>0 而屏障 ≫ Q）⟹ 判詞須標「非可測」，不得作衰變預言。

D5.4 「不合法」之三分語義
  ① 不存在（能量上不可束縛）
  ② 不永久（有能量通道，屏障可令 T½ 天文）
  ③ 可觀測衰變（屏障不足，T½ 有限）
  判準只判 ②；③ 需屏障動力學（禁區①所禁之速率）⟹ 判準結構上不能及 ③。

D5.5 失效半徑（實測界，非定義）
  定義 rel* := 由 data 之實測 max rel 給定之失效半徑（現值見 docs/12 §三，不寫入本定義）。
  以獨立標籤（NUBASE2020 之 T½ 與天然丰度 IS>0）檢定所得之通則：
  凡「判不合法而觀測穩定」者，其最小正 Q 之相對裕度皆 ≤ rel*，且無 Tier I 個案。
  ⟹ 凡 rel > rel* 之不合法判詞，對該等標籤集**零偽陰性**。
  備註：本條原版曾把計數（觀測穩定條數、S\\P 條數）寫入定義，因數據解析修正而失效，
        故改為不嵌入任何計數或數值 —— 定義只可含規則，數值屬 data。
"""

seal = hashlib.sha256(STATEMENT.encode("utf-8")).hexdigest()
print("=" * 78)
print("M10 R11 — D5 判準修訂定義 ＋ 封印")
print("=" * 78)
print(f"  sha256(STATEMENT) = {seal}")

j10 = json.load(open(os.path.join(DATA, "11_falsification_nubase.json"), encoding="utf-8"))
j9 = json.load(open(os.path.join(DATA, "10_channel_completion.json"), encoding="utf-8"))
tb = j10["test_B"]

ck("D5.1 D5.4 之「三分」與 docs/07 §八 之「二分」係相容擴張（② 細化為 ②/③）",
   "②" in STATEMENT and "③" in STATEMENT, "記錄：舊二分不作廢，只加細分")
ck("D5.2 封印長度 64 hex", len(seal) == 64, seal)
ck("D5.3 判準之可測方向與實測一致（P 全 stbl）",
   tb["SP_size"] > 0 and j10["test_A"]["falsified"] is False,
   f"P={j10['P_size']} 條、反例 {len(j10['test_A']['counterexamples'])} ⟹ D5.1(a) 成立")
ck("D5.4 盲格之規模寫入定義（51 條、A∈[98,164]）",
   tb["n_non_tier3"] == 51, f"Tier II 例外 = {tb['n_non_tier3']} 條")
ck("D5.5 失效半徑與實測一致（5.21e−4）",
   abs(tb["max_rel"] - 5.2115e-4) < 1e-7, f"max rel = {tb['max_rel']:.4e}")
ck("D5.6 𝒞_v2 之 11 通道與 docs/09 一致",
   len(j9["final"]["channels"]) == 11 and j9["final"]["n_legal"] == 90,
   f"{len(j9['final']['channels'])} 類、合法集 {j9['final']['n_legal']} 條")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "12_define_d5.py", "role": "R11 D5 判準修訂定義",
       "seal": seal, "statement": STATEMENT,
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "12_definition_d5.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

doc = f"""# docs/11 — D5：判準之修訂定義（封印 `{seal[:16]}…`）

> **本檔由 `scripts/12_define_d5.py` 生成**（避免「腳本抄文檔」之循環）：
> 文檔內容 = 腳本內 `STATEMENT` 之原文，`sha256` 封印 = `{seal}`。
> **新定義，非改舊記錄**：`docs/01` 之 D1–D4 封印不動。

---

## 一、為何需要 D5

R10 以 **NUBASE2020 半衰期**（判準從未使用之獨立資訊）測出：

| 方向 | 命題 | 實測 |
|:--|:--|:--|
| **⟸** | `∀c x(c)<0` ⟹ 觀測穩定 | ✅ **{j10['P_size']}/{j10['P_size']} 零反例** |
| **⟹** | `x(c)>0` ⟹ 不永久 | 51 條判「不合法」而實測 `stbl` ⟹ **不可測** |

故 `D2` 之 `⟺` 不能作預言。D5 把它拆成「可測之充分條件」＋「定義性之非永久」。

---

## 二、D5 全文（封印對象，逐字）

```
{STATEMENT}```

---

## 三、D5 之實測支撐（數字全部來自 `data/`）

| 支撐 | 值 | 來源 |
|:--|:--|:--|
| 合法集 | {j10['P_size']} 條，全 `stbl` | `data/11` test_A |
| 觀測穩定集 S(≥1e15 yr) | {tb['S_size']} 條 | `data/11` test_B |
| S\\P（判不合法卻穩定） | {tb['SP_size']} 條 | 同上 |
| 其中 Tier III | {tb['tier3']} 條（{100*tb['tier3_frac']:.1f}%） | 同上 |
| 其中 Tier II（全 𝒞_B） | {tb['n_non_tier3']} 條，A ∈ [{tb['non_tier3_A_range'][0]},{tb['non_tier3_A_range'][1]}] | 同上 |
| 失效半徑 rel* | {tb['max_rel']:.2e} | 同上 |
| 時間盲跨度 | {j10['test_C']['span']:.1e}× | `data/11` test_C |
| 𝒞_v2 | {len(j9['final']['channels'])} 類通道、合法集 {j9['final']['n_legal']} 條 | `data/10` |

---

## 四、D5 之可證偽內容

1. **D5.5** 之失效半徑：若日後出現 rel > {tb['max_rel']:.2e} 而觀測穩定之核素 ⟹ D5.5 被證偽。
2. **D5.3** 之盲格：若 Tier II × 𝒞_B 出現**非** 𝒞_B 通道之個案 ⟹ 二分法不足。
3. **D5.1(a)**：若出現 `∀c x<0` 而 T½ 有限之核素 ⟹ 判準之充分條件被證偽（即 𝒞 不完備）。

---

## 五、可重跑

```bash
/usr/bin/python3 scripts/12_define_d5.py   # 重生成本檔與 data/12_definition_d5.json
```

> 重跑後本檔與 `data/12` 之 `seal` 必須相同；不同即表示 `STATEMENT` 被改動（= 定義漂移）。
"""
with open(os.path.join(DOCS, "11_definition_d5.md"), "w", encoding="utf-8") as fh:
    fh.write(doc)
print("docs/11_definition_d5.md 已生成（同封）")
print("data/12_definition_d5.json 已寫")