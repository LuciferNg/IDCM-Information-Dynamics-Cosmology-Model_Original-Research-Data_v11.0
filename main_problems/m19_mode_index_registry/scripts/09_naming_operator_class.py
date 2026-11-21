#!/usr/bin/env python3
"""
M19 / scripts/09 — OA1 閉合：逐粒子命名算子類之不存在定理

預註冊（跑前凍結）:
  N1  粒子側**已註冊**標籤集 = {χ = 3(B−L) mod 4, Q（電荷）, B−L}。
       `d ∈ {3,2,0}`（世代分級）已被 M14 `docs/03` 判**範疇錯置**（ℤ₂ 只給 2 類）
      ＋ AST 追溯為**純字面值** ⟹ 不得作標籤。世代（質量 rank）屬輸入（`docs/03` B6）。
  N2  命名算子 = `f: 粒子側標籤 tuple → M9`，要求 `f` 為**雙射**（9 粒子 → 9 模式）。
  N3  判據：`|image(f)| ≤ #distinguishable 粒子類`；若類數 < 9 ⟹ 雙射命名算子**不存在**。
  N4  逐標籤子集計算類數：{χ}､{Q}､{B−L}､{χ,Q}､{χ,B−L}､{Q,B−L}､{χ,Q,B−L}。
  N5  反事實：加入「世代（3 值）」標籤 ⟹ 類數 = ?；若 ≥ 9 ⟹ 阻塞精確定位於「世代標籤未註冊」。
  N6  判：若所有已註冊標籤子集類數 < 9 ⟹ OA1 判 🔴（不存在），並記下唯一救命條件。

失敗判準: N5 若加入世代後類數 < 9 ⟹ OA1 之阻塞不止世代標籤（結論須改寫，停）。
"""
import itertools
import json
from fractions import Fraction
from pathlib import Path

OUT = Path(__file__).parent.parent / "data"
M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
PART = ["e", "mu", "tau", "u", "d", "s", "c", "b", "t"]
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}
BL = {"e": -1, "mu": -1, "tau": -1, "u": Fraction(1, 3), "d": Fraction(-1, 3),
      "s": Fraction(-1, 3), "c": Fraction(1, 3), "b": Fraction(-1, 3), "t": Fraction(1, 3)}
Q = {"e": -1, "mu": -1, "tau": -1, "u": Fraction(2, 3), "d": Fraction(-1, 3),
     "s": Fraction(-1, 3), "c": Fraction(2, 3), "b": Fraction(-1, 3), "t": Fraction(2, 3)}
GEN = {"e": 1, "mu": 2, "tau": 3, "u": 1, "d": 1, "s": 2, "c": 2, "b": 3, "t": 3}

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def chi(p):
    return int(round(3 * BL[p])) % 4


def nclass(labels):
    keys = {tuple(f(p) for f in labels) for p in PART}
    return len(keys), keys


print("=" * 92)
print("M19 scripts/09 — OA1：逐粒子命名算子類之不存在定理")
print("=" * 92)

chk("N1 χ 與 O233 docs/13 表一致",
    [chi(p) for p in PART] == [1, 1, 1, 1, 3, 3, 1, 3, 1],
    " / ".join(f"{LABEL[p]}:{chi(p)}" for p in PART))

subsets = {
    "{χ}": [chi],
    "{Q}": [Q.__getitem__],
    "{B−L}": [BL.__getitem__],
    "{χ,Q}": [chi, Q.__getitem__],
    "{χ,B−L}": [chi, BL.__getitem__],
    "{Q,B−L}": [Q.__getitem__, BL.__getitem__],
    "{χ,Q,B−L}": [chi, Q.__getitem__, BL.__getitem__],
}
rows = {}
for nm, labs in subsets.items():
    n, keys = nclass(labs)
    rows[nm] = n
    print(f"    {nm:<12} 類數 = {n}  → {sorted(str(k) for k in keys)}")
chk("N3/N4 全部已註冊標籤子集類數 < 9",
    all(v < 9 for v in rows.values()), f"max = {max(rows.values())}")

best = max(rows.values())
chk("N4 最強已註冊標籤 = {χ,Q,B−L} 類數 3", rows["{χ,Q,B−L}"] == 3,
    "⟹ 最多 3 個可區分類 ⟹ 不可能雙射 9 模式")

n_gen, keys_gen = nclass([chi, Q.__getitem__, BL.__getitem__, GEN.__getitem__])
chk("N5 加入世代標籤 ⟹ 類數 = 9（恰好可行）", n_gen == 9, f"{n_gen} 類")

chk("N6 阻塞定位 = 世代標籤未註冊（非其他）",
    n_gen == 9 and best == 3,
    "3 類（已註冊）× 3 世代 = 9 = |M9| ⟹ 唯一救命條件就係世代標籤")

# 枚舉：以已註冊標籤定義之所有命名算子（函數空間 = 類 → 模式之任意指派）
# 類數 3 < 9 ⟹ 任何函數 image 至多 3 ⟹ 雙射 0（解析結論）；此處以計數確認
n_bijective = 0
for perm in itertools.permutations(M9):
    m = dict(zip(PART, perm))
    if len({m[p] for p in PART if chi(p) == 1}) <= 1 and \
       len({m[p] for p in PART if chi(p) == 3}) <= 1:
        n_bijective += 1
chk("N6b 以 χ 為唯一標籤之雙射命名數 = 0（枚舉核）", n_bijective == 0, f"{n_bijective}")

OUT.mkdir(exist_ok=True)
(OUT / "09_naming_operator_class.json").write_text(json.dumps({
    "class_counts": rows,
    "max_registered_labels_class": best,
    "with_generation_class": n_gen,
    "generation_label_status": "未註冊（M14 docs/03：範疇錯置＋純字面值）",
    "verdict": "OA1 🔴 不存在命名算子（已註冊標籤下）；唯一救命條件 = 世代標籤註冊",
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
