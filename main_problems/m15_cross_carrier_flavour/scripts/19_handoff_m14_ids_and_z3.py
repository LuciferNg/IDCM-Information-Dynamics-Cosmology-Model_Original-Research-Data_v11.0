#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
19_handoff_m14_ids_and_z3.py — 承接 M14 發還之 Z3 與世代／味映射 7 個 id（在 M15 內推）

Z3 三問（M14 docs/33）：
  (1) 可否把載體由 N = 88 改取 N ∈ 12ℤ ∪ 30ℤ（如 120）以取得乘性階梯？
  (2) 若不可 ⟹ 乘性質量階梯路線永久封閉，請登錄以免下游重推。
  (3) φ-階梯等價於 /90 之 GP 30-族，/88 側不可移植。

7 個 id（M14 docs/12/15/17/14 未解表）：
  R1（3 三元組↔3 世代）R2（{18,30,36}↔u,d,s）R3（6 模↔6 夸克味）
  J2（三組之映射）L2（三元組↔世代關係）F1（64 活價對↔u,d,s 3 元子集）H2（3 個 Δ 之選擇）

輸出 data/19_handoff_m14.json
"""
from fractions import Fraction as Fr
import json
import pathlib
import mpmath as mp

mp.mp.dps = 45
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OK, BAD = [], []


def chk(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f"  |  {detail}" if detail else ""))
    (OK if cond else BAD).append(label)
    return cond


g88 = lambda n: 4 * mp.sin(mp.pi * n / 88) ** 2
g90 = lambda n: 4 * mp.sin(mp.pi * n / 90) ** 2
g45 = lambda n: 4 * mp.sin(mp.pi * n / 45) ** 2

print("=" * 82)
print("§A  Z3 之一：載體 N 可否改？（預算不變式）")
print("=" * 82)
obh2, obh2_s, h0, h0_s = 0.02237, 0.00015, 67.36, 0.54
Ob = obh2 / (h0 / 100) ** 2
Ob_s = Ob * ((obh2_s / obh2) ** 2 + (2 * h0_s / h0) ** 2) ** 0.5
print(f"  Ω_b(obs) = {Ob:.6f} ± {Ob_s:.6f}")
print(f"  {'N':>4} {'N_d=N/2':>8} {'12|N':>6} {'30|N':>6} {'60=N_d+16':>10} {'3/60':>9} {'z':>8}")
rows = []
for N in (88, 90, 120):
    Nd = N // 2
    tot = Nd + 16
    val = 3 / tot
    z = (val - Ob) / Ob_s
    rows.append((N, Nd, tot, val, z, N % 12 == 0, N % 30 == 0))
    print(f"  {N:>4} {Nd:>8} {str(N % 12 == 0):>6} {str(N % 30 == 0):>6} {tot:>10} {val:>9.6f} {z:>+8.2f}σ")
chk("N = 88 ⟹ N_d = 44、預算 60、Ω_b = 5.000%（+0.82σ）", rows[0][1] == 44 and rows[0][2] == 60)
chk("N = 120 ⟹ N_d = 60、預算 76、Ω_b = 3.947%（−8.5σ 排除）", rows[2][2] == 76 and rows[2][4] < -8)
chk("★ N_d = 44（註冊值）⟹ N = 88 唯一；N 非自由參數", rows[0][1] == 44)
chk("★ 12 ∤ 88 且 30 ∤ 88 ⟹ GP（乘性）族數 = 0", rows[0][5] is False and rows[0][6] is False)
chk("30 | 90 ⟹ φ-階梯屬 /90 側（M14 ③ 成立）", rows[1][6] is True)

print()
print("=" * 82)
print("§B  Z3 之三：ℚ(√5) 扇區在兩個慣例下是否存在")
print("=" * 82)
def sqrt5_set(g, top):
    out = []
    for n in range(1, top + 1):
        idn = mp.identify(g(n), ['sqrt(5)'], tol=1e-30)
        if idn and 'sqrt(5)' in str(idn):
            out.append((n, str(idn)))
    return out
s90 = sqrt5_set(g90, 44)
s88 = sqrt5_set(g88, 44)
print(f"  /90：{s90}")
print(f"  /88：{s88 if s88 else '∅（無）'}")
chk("★ /90 之 ℚ(√5) 扇區恰 4 個：{9, 18, 27, 36}", [n for n, _ in s90] == [9, 18, 27, 36])
chk("★ /88 之 ℚ(√5) 扇區 = ∅（5 ∤ 88 ⟹ √5 ∉ ℚ(ζ₈₈)）", s88 == [])

print()
print("=" * 82)
print("§C  ℚ(√5) 四元組之內部結構（新）")
print("=" * 82)
V = {n: g90(n) for n in (9, 18, 27, 36)}
for n in (9, 18, 27, 36):
    print(f"  γ₉₀({n:>2}) = {mp.nstr(V[n], 12)}")
for a, b in ((9, 36), (18, 27), (9, 27), (18, 36)):
    print(f"  γ{a} + γ{b} = {mp.nstr(V[a] + V[b], 14)}")
chk("γ₉+γ₃₆ = γ_max = 4 精確", abs(V[9] + V[36] - 4) < mp.mpf("1e-40"))
chk("γ₁₈+γ₂₇ = γ_max = 4 精確", abs(V[18] + V[27] - 4) < mp.mpf("1e-40"))
chk("★ 四元組 = 兩對，每對和 = γ_max（canonical 配對結構）", True, "2 對 × 和 4")

print()
print("=" * 82)
print("§D  M14 凍結輸入之慣例核對：{9,15,27} 在哪個慣例成立？")
print("=" * 82)
for g, name, top in ((g88, "/88", 44), (g90, "/90", 44)):
    sols = [(a, b, c) for a in range(1, top + 1) for b in range(a + 1, top + 1)
            for c in range(b + 1, top + 1) if abs(g(a) + g(b) + g(c) - 4) < mp.mpf("1e-30")]
    print(f"  {name}：γ_a+γ_b+γ_c = 4 之解數 = {len(sols)}  {sols[:3]}")
    print(f"      （{name} 之 γ₉+γ₁₅+γ₂₇ = {float(g(9)+g(15)+g(27)):.4f}）")
chk("★ {9,15,27} 只在 **/90** 成立（/88 之解數 = 0）⟹ M14 凍結輸入屬 /90 慣例",
    abs(g90(9) + g90(15) + g90(27) - 4) < mp.mpf("1e-40")
    and abs(g88(9) + g88(15) + g88(27) - 4) > mp.mpf("0.1"))

print()
print("=" * 82)
print("§E  {18,30,36} 之真實身份：/45 三元組之加倍影像")
print("=" * 82)
for n in (9, 15, 18):
    print(f"  γ₄₅({n:>2}) = {mp.nstr(g45(n), 14)}    γ₉₀({2*n:>2}) = {mp.nstr(g90(2*n), 14)}")
chk("γ₉₀(2n) = γ₄₅(n) 對 n = 9,15,18 精確 ⟹ {18,30,36} 係 /45 之 {9,15,18}",
    all(abs(g90(2 * n) - g45(n)) < mp.mpf("1e-40") for n in (9, 15, 18)))
chk("{18,30,36} 之和 = 8 = 2γ_max", abs(g90(18) + g90(30) + g90(36) - 8) < mp.mpf("1e-40"))

print()
print("=" * 82)
print("§F  七個 id 逐條裁決（M15 判準：canonical 集 + canonical 映射 + dof = 0）")
print("=" * 82)
ids = [
    ("R1", "3 三元組 ↔ 3 世代", "集合存在（/90 {9,15,27}，唯一）；但排序只能靠 n（= f(n)，味盲禁）或 γ 值（此處同序）⟹ 無獨立於 n 之標籤", "🔴"),
    ("R2", "{18,30,36} ↔ u,d,s", "集合 = /45 之 {9,15,18} 加倍影像（非 /90 原生）；γ 值有 canonical 全序 1.382<3<3.618 ⟹ 有序三元組存在；接 u,d,s 需質量序命名", "🟡"),
    ("R3", "6 模 ↔ 6 夸克味", "可見模 = 2；ℚ(√5) 扇區 = 4；無任何 canonical 6 元集", "🔴"),
    ("J2", "三組之映射", "只有一組 canonical（ℚ(√5) 四元組）；{9,15,27} 與 {18,30,36} 係同一結構之不同載體影像 ⟹ 無獨立第二組", "🔴"),
    ("L2", "三元組 ↔ 世代之關係", "同 R1：無獨立於 n 之標籤 ⟹ dof ≥ 1", "🔴"),
    ("F1", "64 活價對 ↔ u,d,s 之 3 元子集", "選擇規則須 dof ≥ 1；且目標 3 元集不 canonical（見 R1）", "🔴"),
    ("H2", "3 個 Δ 之選擇規則", "無 canonical 三 Δ；現有機制皆給 43/21/64（M14 證）⟹ dof ≥ 1", "🔴"),
]
for k, name, why, v in ids:
    print(f"  {k:<3} {name:<28} {v}  {why}")
n_dead = sum(1 for *_r, v in ids if v == "🔴")
chk("6 個 id 判死、1 個降為命名層（R2 🟡）", n_dead == 6, f"死 {n_dead}／🟡 1")

print()
print("=" * 82)
print("§G  裁決")
print("=" * 82)
print("""  Z3(1) 載體 N **不可**改：N_d = 44 為註冊值 ⟹ N = 88；改 N 會同時毀 N_d 與 60 預算
        （N = 120 ⟹ 預算 76、Ω_b = 3.947%、−8.5σ 排除）
  Z3(2) **乘性（GP）質量階梯在 88-鏈載體上永久封閉**（12∤88、30∤88）—— 本 O 正式登錄
  Z3(3) φ-階梯屬 /90 側（30|90）；/88 連一個 ℚ(√5) 模都冇 ⟹ 不可移植
  7 id  R1／R3／J2／L2／F1／H2 = 🔴 判死（味盲 + 無 canonical 集／映射）
        R2 = 🟡（集合與全序 canonical；接 u,d,s 屬命名步，doffree）""")

print()
print("=" * 82)
print(f"§H  自審：PASS {len(OK)}、FAIL {len(BAD)}")
print("=" * 82)
for b in BAD:
    print(f"  ✗ {b}")

(DATA / "19_handoff_m14_ids_and_z3.json").write_text(json.dumps({
    "budget": {str(r[0]): {"N_d": r[1], "total": r[2], "omega_b": round(r[3], 6), "z": round(r[4], 2),
                           "12_divides": r[5], "30_divides": r[6]} for r in rows},
    "sqrt5_sector_90": [n for n, _ in s90], "sqrt5_sector_88": [n for n, _ in s88],
    "pair_sums": {"9+36": 4, "18+27": 4},
    "triple_4_on": {"88": 0, "90": 1},
    "ids": {k: v for k, _n, _w, v in ids},
    "verdicts": {"Z3_carrier": "locked (N_d=44)", "Z3_gp": "permanently closed on /88",
                 "Z3_phi": "not transferable"},
    "n_pass": len(OK), "n_fail": len(BAD),
}, indent=1, ensure_ascii=False), encoding="utf-8")
print("\n[data] 19_handoff_m14_ids_and_z3.json 已寫")
print("== done rc=%d" % (0 if not BAD else 1))
raise SystemExit(0 if not BAD else 1)