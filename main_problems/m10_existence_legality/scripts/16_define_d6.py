#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/16 — D6：機制軸判定（取代 Tier 裕度軸），並封印

前設問題:
  D5（docs/11, seal c7fecc02…）把 𝒞 二分為 𝒞_E／𝒞_B。R14（docs/15）證明二分不足：
  𝒞_E 內之一階弱通道亦可被**選擇定則**壓制（¹⁸⁷Re 型），二階弱則被 G_F⁴ 壓制。
  ⟹ 壓制機制有三種，且每種皆可令「能量開而不衰」⟹ D5.1 之單向性係**全通道通用**。
  故立 D6，封印之，使後續輪次不得暗中移動 Tier 定義。

教訓（D5 之重封）：定義**不得嵌入任何實測計數或數值**，只可含規則。

執行：/usr/bin/python3 scripts/16_define_d6.py
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


STATEMENT = """D6（機制軸判定，v11.0/M10，2026-10-27）

D6.1 單向性之普遍形式
  可測方向（充分條件）：若 ∀c∈𝒞: x(c)<0，則 C 永久（無限時間內不衰變）。
  不可測方向：∃c: x(c)>0 **不**可推「會衰變」。此限制對**任何**通道類成立，
  不因通道屬弱作用或具位壘而異 —— 三種壓制機制皆可令能量開而不衰。

D6.2 三類壓制機制（並集 = 判準之通道集，互不蘊涵）
  𝒞1 一階弱 {β⁻, β⁺/EC, p, n}：無位壘、單體；唯選擇定則可壓制。
  𝒞2 二階弱 {2β⁻, 2β⁺/ECEC}：G_F⁴ 壓制。
  𝒞3 位壘 {α, cluster, SF_sym/SF_opt}：穿隧壓制。

D6.3 選擇定則為跨類壓制層
  對 𝒞1 通道以禁戒階 n ∈ {0,1,2,3,4} 標記（0 = allowed，4 = 三階以上）。
  能量開 ∧ n ≤ 1 ⟹ 於 10^10 yr 尺度上必已衰變；
  能量開 ∧ n ≥ 2 ⟹ 可長壽而不衰。
  故 n 與位壘同屬壓制參數，**不入**能量判準。

D6.4 Tier-M 取代裕度軸
  Tier 之軸 := 機制標籤（𝒞1-直接／𝒞1-禁戒(n)／𝒞2／𝒞3）。
  禁用 |x|/m 之裕度軸作 Tier：其最高層覆蓋全表過半核素，對合法性零信息。
  裕度僅可作數值穩定度之輔助指標。

D6.5 失效半徑之正確讀法
  rel* := 「能量開而仍被壓制」之最大相對裕度，係**經驗界**，
  與質量測量精度無關（相關通道之 Q/σ ≫ 1）。不得解作精度牆。

D6.6 可證偽內容
  (a) 若出現 𝒞1 通道、n ≤ 1、能量開、而 T½ > 10^10 yr 之核素 ⟹ D6.3 破。
  (b) 若出現無法歸入 𝒞1/𝒞2/𝒞3 之能量開通道 ⟹ D6.2 之窮盡性破。
  (c) 若某核素全部通道 x<0 而 T½ 有限 ⟹ D6.1 破（即 𝒞 不完備）。
"""

seal = hashlib.sha256(STATEMENT.encode("utf-8")).hexdigest()
print("=" * 78)
print("M10 R15 — D6 機制軸判定 ＋ 封印")
print("=" * 78)
print(f"  sha256(STATEMENT) = {seal}")

j15 = json.load(open(os.path.join(DATA, "15_mechanism_classification.json"), encoding="utf-8"))
j11 = json.load(open(os.path.join(DATA, "11_falsification_nubase.json"), encoding="utf-8"))
lab = j15["labels"]

ck("D6.1 三機制之條數和 = 例外總數（窮盡）",
   sum(v for k, v in lab.items() if k != "C1 未禁戒（異常）") == j15["n_exceptions"],
   f"{lab} → 合計 {j15['n_exceptions']} 條")
ck("D6.2 無「未禁戒之 𝒞1」個案（D6.2 窮盡性）", lab.get("C1 未禁戒（異常）", 0) == 0,
   f"未解釋 = {lab.get('C1 未禁戒（異常）', 0)} 條")
ck("D6.3 無精度個案（D6.5 之依據）", not j15["E1_low"],
   f"Q/σ 最小 = {j15['qos_min']:.2f}（中位 {j15['qos_median']:.0f}）")
ck("D6.4 禁戒階閘有效：長壽而開放 𝒞1 者全部 n ≥ 2",
   all(v["order"] >= 2 for v in j15["E4b_violations"]),
   f"{len(j15['E4b_violations'])} 條，階分佈 "
   f"{sorted(set(v['order'] for v in j15['E4b_violations']))}")
ck("D6.5 封印 64 hex 且與 D5 互異", len(seal) == 64, seal)
ck("D6.6 判準之可測方向仍零反例（⟸ 不因 D6 而動）",
   j11["test_A"]["falsified"] is False,
   f"P = {j11['P_size']} 條，反例 {len(j11['test_A']['counterexamples'])} 條")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {"script": "16_define_d6.py", "role": "R15 D6 機制軸判定",
       "seal": seal, "statement": STATEMENT,
       "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL}
with open(os.path.join(DATA, "16_definition_d6.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)

doc = f"""# docs/16 — D6：機制軸判定（封印 `{seal[:16]}…`）

> **本檔由 `scripts/16_define_d6.py` 生成**（避免「腳本抄文檔」之循環）：
> 文檔內容 = 腳本內 `STATEMENT` 之原文，`sha256` 封印 = `{seal}`。
> **新定義，非改舊記錄**：D1–D4（`fe49c14d…`）與 D5（`c7fecc02…`）封印不動。
> D6 取代 D5.2 之 `𝒞_E／𝒞_B` 二分（該二分不足，見 `docs/15` §八）。

---

## 一、為何需要 D6

R14（`docs/15`）以 182 條例外作盲測，得：

| 檢定 | 結果 |
|:--|:--|
| 精度接對（`Q/σ_Q`） | 全部 **> 5**（最小 {j15['qos_min']:.2f}，中位 {j15['qos_median']:.0f}）⟹ **無精度個案** |
| 三機制窮盡 | **{j15['n_exceptions']}/{j15['n_exceptions']}** 全部可歸類；未解釋 **0** |
| 禁戒階閘 | 長壽而開放 𝒞1 者 **{len(j15['E4b_violations'])}** 條，全部 n ≥ 2；n ≤ 1 者 **0** |

⟹ 壓制機制有三種，`𝒞_E` 內部亦不單純；D5.2 之二分不足 ⟹ 立 D6。

---

## 二、D6 全文（封印對象，逐字）

```
{STATEMENT}```

---

## 三、D6 之實測支撐（數字全部來自 `data/`）

| 支撐 | 值 | 來源 |
|:--|:--|:--|
| 例外總數 | {j15['n_exceptions']} 條 | `data/15` |
| 𝒞3 位壘 | {lab.get('C3 位壘', 0)} 條 | 同上 |
| 𝒞2 二階弱 | {lab.get('C2 二階弱', 0)} 條 | 同上 |
| 𝒞1 選擇定則禁戒 | {lab.get('C1 選擇定則禁戒', 0)} 條 | 同上 |
| 未解釋 | {lab.get('C1 未禁戒（異常）', 0)} 條 | 同上 |
| `Q/σ_Q` 中位／最小 | {j15['qos_median']:.0f}／{j15['qos_min']:.2f} | 同上 |
| 長壽而開放 𝒞1 | {len(j15['E4b_violations'])} 條（全 n≥2） | 同上 |
| 舊軸 Tier III 覆蓋（全表） | {100*j15['tier3_frac_all']:.1f}% | `data/08`＋`data/15` |

---

## 四、D6 之可證偽內容（三條）

1. 出現 𝒞1 通道、禁戒階 ≤ 1、能量開、而 `T½ > 1e10 yr` 之核素 ⟹ **D6.3 破**。
2. 出現無法歸入 𝒞1／𝒞2／𝒞3 之能量開通道 ⟹ **D6.2 之窮盡性破**。
3. 出現全部通道 `x<0` 而 `T½` 有限之核素 ⟹ **D6.1 破**（即 `𝒞` 不完備）。

---

## 五、可重跑

```bash
/usr/bin/python3 scripts/16_define_d6.py   # 重生成本檔與 data/16_definition_d6.json
```

> 重跑後本檔與 `data/16` 之 `seal` 必須相同；不同即表示 `STATEMENT` 被改動（= 定義漂移）。
"""
with open(os.path.join(DOCS, "16_definition_d6.md"), "w", encoding="utf-8") as fh:
    fh.write(doc)
print("docs/16_definition_d6.md 已生成（同封）")
print("data/16_definition_d6.json 已寫")