#!/usr/bin/env python3
"""
M19 / scripts/07 — 槽定理之作用域定理（T3 只能單向，係算術強制）

預註冊（跑前凍結）:
  H1  χ 由 T3 給：χ = 3·(B−L) mod 4（輕子 B−L = −1 ⟹ χ = 1；夸克 ±1/3 ⟹ χ ∈ {1,3}）
      全 9 粒子之 χ 值與 v10 O233 docs/13 表逐字對照
  H2  模式側 ℤ₄ 類大小（9 模式集內）：c1 = |{1,5,9,13}| = 4；c2 = |{2}| = 1；
      c3 = |{3,7,11}| = 3；c0 = |{4}| = 1
  H3  「全 T3」（χ ≡ n mod 4 對全部 9 粒子施加）之可行雙射數 = ?
  H4  「單向 T3」（只對輕子施加）之可行雙射數 = 17280（= scripts/02 之 ρ_struct）
  H5  判：若 H3 = 0 ⟹ 全 T3 唔可能，作用域受限係**算術強制**（非為遷就 F1 而設）
      若 H3 > 0 ⟹ 全 T3 可行，則「夸克無槽」係選擇而非定理（結論相反）
  H6  若 H3 = 0，須指出確切阻塞量：χ=1 之粒子數（6）> ≡1 模式數（4）

失敗判準: H1 與 O233 docs/13 表不符 ⟹ χ 定義讀錯（停）。
"""
import itertools
import json
from pathlib import Path

OUT = Path(__file__).parent.parent / "data"
M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
PART = ["e", "mu", "tau", "u", "d", "s", "c", "b", "t"]
LEP = {"e", "mu", "tau"}
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}
# B−L（夸克 ±1/3，輕子 −1）
BL = {"e": -1, "mu": -1, "tau": -1, "u": 1 / 3, "d": -1 / 3, "s": -1 / 3,
      "c": 1 / 3, "b": -1 / 3, "t": 1 / 3}
# v10 O233 docs/13 表之 χ 值（逐字）
CHI_REF = {"e": 1, "mu": 1, "tau": 1, "u": 1, "d": 3, "s": 3, "c": 1, "b": 3, "t": 1}

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/07 — 槽定理之作用域定理")
print("=" * 92)

chi = {}
for p in PART:
    v = 3 * BL[p]
    assert abs(v - round(v)) < 1e-12, p
    chi[p] = int(round(v)) % 4

chk("H1 χ = 3(B−L) mod 4 與 O233 docs/13 表逐字一致",
    chi == CHI_REF, " / ".join(f"{LABEL[p]}:{chi[p]}" for p in PART))

class_sizes = {c: sum(1 for n in M9 if n % 4 == c) for c in range(4)}
part_sizes = {c: sum(1 for p in PART if chi[p] == c) for c in range(4)}
chk("H2 模式類 = {c0:1, c1:4, c2:1, c3:3}", class_sizes == {0: 1, 1: 4, 2: 1, 3: 3}, str(class_sizes))
print(f"    粒子 χ 類：{part_sizes}")

n_full = 0
for perm in itertools.permutations(M9):
    m = dict(zip(PART, perm))
    if all(m[p] % 4 == chi[p] for p in PART):
        n_full += 1
chk("H3 全 T3（χ ≡ n mod 4 全施加）可行雙射 = 0", n_full == 0, f"{n_full}")

n_one = 0
for perm in itertools.permutations(M9):
    m = dict(zip(PART, perm))
    if all(m[p] in {1, 5, 9, 13} for p in LEP):
        n_one += 1
chk("H4 單向 T3（只施加於輕子）= 17280", n_one == 17280, f"{n_one:,}")

chk("H5 作用域受限 = 算術強制（H3 = 0）",
    n_full == 0 and n_one > 0,
    "全 T3 不可滿足 ⟹ 唔可以話「夸克無槽」係遷就 F1 之選擇")

chk("H6 阻塞量：χ=1 粒子 6 > ≡1 模式 4（鴿籠）",
    part_sizes[1] == 6 and class_sizes[1] == 4,
    f"{part_sizes[1]} 粒子 vs {class_sizes[1]} 槽 ⟹ 缺口 {part_sizes[1]-class_sizes[1]}")

OUT.mkdir(exist_ok=True)
(OUT / "07_slot_scope.json").write_text(json.dumps({
    "chi": {LABEL[p]: chi[p] for p in PART},
    "mode_class_sizes": {str(k): v for k, v in class_sizes.items()},
    "particle_class_sizes": {str(k): v for k, v in part_sizes.items()},
    "n_full_T3": n_full,
    "n_one_sided_T3": n_one,
    "pigeonhole_gap": part_sizes[1] - class_sizes[1],
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
