#!/usr/bin/env python3
"""
M19 / scripts/03 — 容許矩陣 ＋ 三讀裁決（T2）＋ c/b/t 判（T3）

預註冊（跑前凍結）:
  D1  容許矩陣 A[p][n] = #{S1 可行雙射 : p ↦ n}，全枚舉 362,880（無剪枝）
  D2  輕子列支撐 = {1,5,9,13}（類別強制）；夸克列支撐 = **全部 9 個模式**（無槽）
      ⟹ 結構內容止於**輕子類別**，夸克值層完全自由
  D3  T2（s 三讀 {4,5,7}）判準 = 本 record §二 T2 之凍結失敗判準：
      「若三讀皆可容 ⟹ 零判別力，🔴」。三讀之 S1 可行數全部 > 0 ⟹ 判 🔴。
  D4  三讀之取捨只能由**質量 rank**（輸入）給：u<d<s<c<b ⟹ 模式 2<3<4<7<11 ⟹ s=4。
      該 rank 步驟使用觀測質量 ⟹ 屬輸入，不得記為結構推導。
  D5  T3（c/b/t）：三者皆夸克 ⟹ T3 唔施加 ⟹ 值層零結構約束；F1 之 7/11/13 全屬 rank。
  D6  t→13 唔係「≡1 剩餘槽」：夸克可落任何模式（D2），故 t=13 亦係 rank 產物。
      此點推翻「t 落剩餘槽 ⟹ 半結構」之直覺讀法，必須明寫。

失敗判準: D1 矩陣總和 ≠ 17280 ⟹ 枚舉與 02 不一致（停）；D6 若 t 有結構強制 ⟹ T3 判決改寫。
"""
import itertools
import json
import math
from pathlib import Path

OUT = Path(__file__).parent.parent / "data"

M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
PART = ["e", "mu", "tau", "u", "d", "s", "c", "b", "t"]
LEP = {"e", "mu", "tau"}
QUARK = ["u", "d", "s", "c", "b", "t"]
MASS = {"e": 0.511, "u": 2.16, "d": 4.67, "s": 93.4, "mu": 105.6584,
        "c": 1270.0, "tau": 1776.86, "b": 4180.0, "t": 172690.0}
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}
S1 = {n for n in M9 if n % 4 == 1}
F1 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/03 — 容許矩陣 ＋ 三讀裁決 ＋ c/b/t 判")
print("=" * 92)

A = {p: {n: 0 for n in M9} for p in PART}
tot = 0
for perm in itertools.permutations(M9):
    m = dict(zip(PART, perm))
    if all(m[p] in S1 for p in LEP):
        tot += 1
        for p in PART:
            A[p][m[p]] += 1

chk("D1 容許矩陣總和 = ρ_struct = 17280", tot == 17280, f"{tot:,}")

print("\n  容許矩陣 A[p][n]（列 = 粒子，欄 = 模式）：")
print("  " + f"{'':>4}" + "".join(f"{n:>7}" for n in M9))
for p in PART:
    print("  " + f"{LABEL[p]:>4}" + "".join(f"{A[p][n]:>7}" for n in M9))

lep_support = {p: sorted(n for n in M9 if A[p][n] > 0) for p in LEP}
q_support = {p: sorted(n for n in M9 if A[p][n] > 0) for p in QUARK}
chk("D2a 輕子列支撐 = {1,5,9,13}", all(lep_support[p] == sorted(S1) for p in LEP), str(sorted(S1)))
chk("D2b 夸克列支撐基數 = 9（全模，無槽）", all(len(q_support[p]) == 9 for p in QUARK),
    f"每夸克支撐 {sorted(set(len(v) for v in q_support.values()))} 個模式（≡1 欄權重 720，非 ≡1 欄 2880）")
chk("D2c 輕子列於非 ≡1 模式全零", all(A[p][n] == 0 for p in LEP for n in M9 if n not in S1), "6/9 欄 = 0")

s_counts = {n: A["s"][n] for n in (4, 5, 7)}
chk("D3 T2 凍結判準：三讀皆可容 ⟹ 零判別力",
    all(c > 0 for c in s_counts.values()),
    "s→4/5/7 可行數 = " + ", ".join(f"{n}:{c:,}" for n, c in s_counts.items()) + " ⟹ 🔴（判準觸發）")

rank_ok = F1["s"] == 4
chk("D4 質量 rank 取 s=4（輸入級）", rank_ok,
    "u2<d3<s4<c7<b11 ⟸ 質量 u<d<s<c<b（用觀測值）")

chk("D5 c/b/t 值層零結構約束（夸克支撐 = 全 9 模）",
    all(len(q_support[p]) == 9 for p in ("c", "b", "t")),
    "c7/b11/t13 全屬 rank 輸入（結構對夸克值層零約束）")
chk("D6 t→13 唔係結構強制（支撐 9 模；13 之可行數 720 ≠ 非 ≡1 模之 2880）",
    A["t"][13] == 720 and A["t"][2] == 2880,
    "差異係「佔用 ≡1 槽 ⇒ 輕子少一個擺位」之計數效應，非結構必要性")

# 每粒子「結構剩餘自由度」（條件於其他粒子）——只作展示，位元以 02 之全域分解為準
resid = {}
for p in PART:
    sup = [n for n in M9 if A[p][n] > 0]
    resid[LABEL[p]] = {"support": sup, "log2_free": math.log2(len(sup))}

print("\n  每粒子值層剩餘自由度（log₂ 支撐大小，僅展示；全域位元見 02）：")
for p in PART:
    print(f"    {LABEL[p]:>2}: 支撐 {resid[LABEL[p]]['support']}  log₂={resid[LABEL[p]]['log2_free']:.4f}")

OUT.mkdir(exist_ok=True)
(OUT / "03_admissibility.json").write_text(json.dumps({
    "A": {LABEL[p]: {str(n): A[p][n] for n in M9} for p in PART},
    "rho_struct": tot,
    "lepton_support": sorted(S1),
    "quark_support": q_support,
    "s_counts": {str(k): v for k, v in s_counts.items()},
    "resid": resid,
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
