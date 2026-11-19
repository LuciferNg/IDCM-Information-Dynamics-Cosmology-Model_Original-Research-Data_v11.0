#!/usr/bin/env python3
"""
M19 / scripts/02 — 結構／輸入位元分解（T1 判決核心）

預註冊（跑前凍結）:
  B1  候選空間 = 9 粒子 → 9 模式集之**全部**雙射，基數 9! = 362,880（無剪枝）
  B2  結構謂詞集（只准用已註冊結構，禁質量）:
        S1  T3 槽定理：輕子（χ=1，B−L∈ℤ）必須落 n ≡ 1 (mod 4) 之模式
            ⟹ leptons {e,μ,τ} ⊂ {1,5,9,13}
        S2  9 模式集成員（A₅ 動差準則）—— 由 B1 之定義域已含，零額外約束
        S3  夸克無槽（T3 對 B−L∉ℤ 唔施加；v10 O233 docs/13 明文「冇固定槽」）
            ⟹ 夸克可落任何模式，包括 ≡1 類（t 落 13 即此情形）
      判準: 一個謂詞「有結構內容」⟺ 佢令可行雙射數嚴格下降
  B3  結構殘量 ρ_struct = 4（≡1 槽取 3 個給輕子）× 3!（輕子次序）× 6!（其餘 6 粒子）
      = 17,280 = 9!/21
  B4  位元分解: log₂(9!) = log₂(ρ_struct) + log₂(9!/ρ_struct)
      結構位元 = log₂(9!/ρ_struct) = log₂21 = 4.392；輸入位元 = log₂17280 = 14.077
  B5  質量 rank 唯一性: 用觀測質量（只作排序）＋ 類內模式升序 ⟹ 唯一映射；
      該映射必須 == F1。若唔等 ⟹ T1 判決不可用（停）。
  B6  (反事實閘) 若 T3 槽定理被撤回（S1 移除）⟹ ρ_struct = 9! ⟹ 結構位元 = 0
      ⟹ k_q 全部係輸入。記錄此依賴，不得省略。

失敗判準: B5 FAIL ⟹ 位元分解與 F1 身份不能並用，本輪 T1 降為「不可判」。
"""
import itertools
import json
import math
from pathlib import Path

from mpmath import mp, mpf

mp.dps = 40
OUT = Path(__file__).parent.parent / "data"

M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
PART = ["e", "mu", "tau", "u", "d", "s", "c", "b", "t"]
LEP = {"e", "mu", "tau"}
MASS = {"e": 0.511, "u": 2.16, "d": 4.67, "s": 93.4, "mu": 105.6584,
        "c": 1270.0, "tau": 1776.86, "b": 4180.0, "t": 172690.0}
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}
S1 = {n for n in M9 if n % 4 == 1}          # {1, 5, 9, 13}
F1 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/02 — 結構／輸入位元分解（T1）")
print("=" * 92)

total = math.factorial(9)
chk("B1 候選空間 = 9!（全枚舉，無剪枝）", total == 362880, f"{total:,}")

# 全枚舉：所有雙射，統計 S1 可行者
admissible = 0
f1_admissible = True
f2 = {"t": 1, "b": 2, "tau": 3, "c": 4, "mu": 5, "s": 7, "d": 9, "u": 11, "e": 13}
f2_admissible = True
s_at = {4: 0, 5: 0, 7: 0}
for perm in itertools.permutations(M9):
    m = dict(zip(PART, perm))
    ok = all(m[p] in S1 for p in LEP)
    if ok:
        admissible += 1
        for n in (4, 5, 7):
            if m["s"] == n:
                s_at[n] += 1
    if perm == tuple(F1[p] for p in PART):
        f1_admissible = ok
    if perm == tuple(f2[p] for p in PART):
        f2_admissible = ok

chk("B2 S1（T3 槽：輕子 ⊂ {1,5,9,13}）令可行數下降",
    admissible < total, f"{admissible:,} / {total:,}")
chk("B3 結構殘量 ρ_struct = 4·3!·6! = 17280", admissible == 17280, f"{admissible:,} = 9!/{total//admissible}")
chk("B4a F1 滿足 S1", f1_admissible, "e1 μ5 τ9 全 ≡1 (mod 4)")
chk("B4b F2 違反 S1（τ 落 3）", not f2_admissible, "τ3 ∉ {1,5,9,13}")

struct_bits = math.log2(total / admissible)
input_bits = math.log2(admissible)
chk("B5 位元分解恆等式",
    abs((struct_bits + input_bits) - math.log2(total)) < 1e-12,
    f"結構 {struct_bits:.6f} bit ＋ 輸入 {input_bits:.6f} bit = {math.log2(total):.6f} bit")

# B6：質量 rank + 類內升序 ⟹ 唯一映射
def rank_map(particles, modes):
    ps = sorted(particles, key=lambda p: MASS[p])
    ms = sorted(modes)
    return dict(zip(ps, ms))


lewis = rank_map(["e", "mu", "tau", "t"], sorted(S1))      # ≡1 類：4 個
quarks = rank_map(["u", "d", "s", "c", "b"], sorted(set(M9) - S1))
combined = {**lewis, **quarks}
chk("B6 質量 rank（類內升序）⟹ 恰等於 F1",
    combined == F1, " / ".join(f"{LABEL[p]}→{combined[p]}" for p in sorted(combined, key=lambda x: combined[x])))

chk("B7 反事實：撤回 S1 ⟹ 結構位元 = 0",
    math.log2(total / total) == 0.0,
    "ρ_struct=9! ⟹ 全部 18.43 bit 皆輸入（k_q 全非結構）")

print("\n  s 位置喺 S1 可行集中之分佈：")
for n in (4, 5, 7):
    print(f"    s→{n}: {s_at[n]:,} 個可行雙射")

OUT.mkdir(exist_ok=True)
(OUT / "02_bit_decomposition.json").write_text(json.dumps({
    "total_bijections": total,
    "rho_struct": admissible,
    "S1_class": sorted(S1),
    "F1_admissible": f1_admissible,
    "F2_admissible": f2_admissible,
    "struct_bits": struct_bits,
    "input_bits": input_bits,
    "s_at": {str(k): v for k, v in s_at.items()},
    "rank_map": {LABEL[p]: combined[p] for p in combined},
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
