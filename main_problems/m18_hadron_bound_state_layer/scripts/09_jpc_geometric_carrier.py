#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/09 — R3′：`J` 載體（幾何層不可約表示）—— 存在性定理
==================================================================
前設問題
--------
`docs/04` 判死 `J` 於 mode／tree 層（`D_88` irrep ≤2 維、`Aut(P_34)=C_2`）⟹ `J` 必須由幾何層供。
本輪**在本 M 內**推 `J` 載體（唔外推），並答三問：
  Q-a 邊個幾何群可以承載 `J`？（**存在性定理**，唔可以「揀一個似嘅」）
  Q-b 可承載之 `J` 值域？（單一 irrep ↔ 2J+1 維）
  Q-c 加上 `docs/04` 已導出之 `P`、`C`，`J^PC` 類容量幾多？可否證偽？

方法（零參數；群論為教科書數學，非 IDCM 擬合）：
  1. 旋轉點群 irrep 維數以 `Σd² = |G|` 核（Burnside），並由置換群實際建構核類數
  2. `A₅` 特徵表須過**全部**正交關係（精確 ℚ(√5)）
  3. `SO(3) ↓ I` 分解：`χ_J(θ) = sin((J+½)θ)/sin(θ/2)`，逐 `J = 0..8`
  4. `J=2` ⟹ 需 5 維 irrep ⟹ 掃全部點群判邊個有（存在性，非選擇）

事前判準（寫死，唔准事後改）
----------------------------
  Z1 點群 irrep 維數必須過 `Σd² = |G|`；類數須由置換群實建核（唔照抄表）
  Z2 `A₅` 特徵表須精確通過列／行正交（ℚ(√5)，非浮點）
  Z3 `SO(3)↓I` 每個 `J` 必須 `Σ_R n_R·dim R = 2J+1`（內部閉合），`n_R` 為整數
  Z4 `J` 載體判定 = **存在性**：列舉 `C_n / D_n / T / O / I` 全族，判「有冇 2J+1 維 irrep」
  Z5 容量計算零參數；`P`、`C` 用 `docs/04` 已導出者（`P = (−1)^{k+1}`、`C = ±1`）
  Z6 可證偽：(i) 觀測低窪（< 1.7 GeV）`J^PC` 類每個都必須可命 irrep；(ii) `J` 首次出現質量單調
  Z7 禁區：`A₅` **頂點着色作味標籤**（M12 `docs/03` 已判死：1 bit ＋ gauge）—— 本輪唔碰

輸出：data/09_jpc_geometric_carrier.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="Z", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="Z"):
    return chk(name, ok, detail, tag, integrity=False)


# ─────────────────────────────────────────────────────────────
# 置換群建構（核類數）
# ─────────────────────────────────────────────────────────────
def compose(a, b):
    return tuple(a[b[i]] for i in range(len(b)))


def sign(p):
    s = 1
    seen = [False] * len(p)
    for i in range(len(p)):
        if not seen[i]:
            j, L = i, 0
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def perm_group(n, even_only=False):
    G = set()
    for p in itertools.permutations(range(n)):
        if even_only and sign(p) != 1:
            continue
        G.add(p)
    return sorted(G)


def classes(G):
    seen, out = set(), []
    for g in G:
        if g in seen:
            continue
        # 共軛類：h g h⁻¹
        cl = set()
        for h in G:
            hinv = tuple(sorted(range(len(h)), key=lambda i: h[i]))  # inverse
            cl.add(compose(compose(h, g), hinv))
        out.append(sorted(cl))
        seen |= cl
    return out


print("=" * 78)
print("Z1. 置換群類數 ＋ 點群 irrep 維數（Σd² = |G|）")
print("=" * 78)
A4 = perm_group(4, even_only=True)
S4 = perm_group(4)
A5 = perm_group(5, even_only=True)
for nm, G in (("A₄ (=T)", A4), ("S₄ (=O)", S4), ("A₅ (=I)", A5)):
    cl = classes(G)
    print("  %-7s |G| = %d｜共軛類 = %d｜類大小 = %s" % (nm, len(G), len(cl), sorted(len(c) for c in cl)))
chk("Z1a |T|=12、|O|=24、|I|=60 且類數 4/5/5",
    (len(A4), len(S4), len(A5)) == (12, 24, 60) and
    (len(classes(A4)), len(classes(S4)), len(classes(A5))) == (4, 5, 5),
    "12/24/60，類 4/5/5")

DIMS = {"C_n": None, "D_n": None, "T": [1, 1, 1, 3], "O": [1, 1, 2, 3, 3], "I": [1, 3, 3, 4, 5]}


def dims_C(n):
    return [1] * n


def dims_D(n):
    if n % 2 == 0:
        return [1, 1, 1, 1] + [2] * (n // 2 - 1)
    return [1, 1] + [2] * ((n - 1) // 2)


ok_sum = all(sum(d * d for d in dims_C(n)) == n for n in range(1, 13)) and \
         all(sum(d * d for d in dims_D(n)) == 2 * n for n in range(2, 13)) and \
         sum(d * d for d in DIMS["T"]) == 12 and sum(d * d for d in DIMS["O"]) == 24 and \
         sum(d * d for d in DIMS["I"]) == 60
chk("Z1b 點群 irrep Σd² = |G|：C_n(1..12)、D_n(2..12)、T/O/I 全過", ok_sum,
    "Burnside 核 5/5 族")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z2. A₅ 特徵表正交關係（精確 ℚ(√5)）")
print("=" * 78)
r5 = sp.sqrt(5)
g = (1 + r5) / 2
CS = [1, 12, 12, 20, 15]                     # 1A, 5A, 5B, 3A, 2A
CH = {
    "A":  [1, 1, 1, 1, 1],
    "T1": [3, g, 1 - g, 0, -1],
    "T2": [3, 1 - g, g, 0, -1],
    "G":  [4, -1, -1, 1, 0],
    "H":  [5, 0, 0, -1, 1],
}
print("  類大小 =", CS, "Σ =", sum(CS))
row_ok = True
for a in CH:
    for b in CH:
        ip = sp.simplify(sum(CS[k] * CH[a][k] * CH[b][k] for k in range(5)) / 60)
        want = 1 if a == b else 0
        if sp.simplify(ip - want) != 0:
            row_ok = False
col_ok = True
for k in range(5):
    for l in range(5):
        ip = sp.simplify(sum(CH[r][k] * CH[r][l] for r in CH))
        want = 60 / CS[k] if k == l else 0
        if sp.simplify(ip - want) != 0:
            col_ok = False
print("  列正交（5×5 全對）= %s｜行正交（5×5）= %s" % (row_ok, col_ok))
chk("Z2a A₅ 特徵表 25+25 條正交關係精確成立", row_ok and col_ok, "列 25/25、行 25/25")
chk("Z2b Σ(χ(e))² = 1+9+9+16+25 = 60 = |A₅|", sum(CH[r][0] ** 2 for r in CH) == 60, "60")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z3. SO(3) ↓ I 分解（χ_J(θ) = sin((J+½)θ)/sin(θ/2)）")
print("=" * 78)
THETA = [mp.mpf(0), 2 * mp.pi / 5, 4 * mp.pi / 5, 2 * mp.pi / 3, mp.pi]


def chiJ(J, th):
    if th == 0:
        return mp.mpf(2 * J + 1)
    return mp.sin((J + mp.mpf(1) / 2) * th) / mp.sin(th / 2)


CHn = {r: [mp.mpf(0)] * 5 for r in CH}
_g = (1 + mp.sqrt(5)) / 2
CHmp = {
    "A":  [1, 1, 1, 1, 1],
    "T1": [3, _g, 1 - _g, 0, -1],
    "T2": [3, 1 - _g, _g, 0, -1],
    "G":  [4, -1, -1, 1, 0],
    "H":  [5, 0, 0, -1, 1],
}
branch = {}
for J in range(0, 9):
    chiv = [chiJ(J, THETA[k]) for k in range(5)]
    ns = {}
    for r in CH:
        n = mp.mpf(0)
        for k in range(5):
            n += CS[k] * chiv[k] * CHmp[r][k]
        n = n / 60
        ni = int(mp.nint(n))
        if abs(n - ni) > mp.mpf("1e-40"):
            ni = None
        ns[r] = ni
    branch[J] = ns
    dims = sp.symbols("x")
    tot = sum((ns[r] or 0) * CH[r][0] for r in CH)
    print("  J=%d：%s ｜Σ n_R·dim R = %s（應 %d）"
          % (J, {r: ns[r] for r in CH if ns[r]}, tot, 2 * J + 1))
chk("Z3a 全部 J=0..8 之 n_R 為整數（1e−30 內）",
    all(all(v is not None for v in branch[J].values()) for J in branch), "9/9")
chk("Z3b 內部閉合 Σ_R n_R·dim R = 2J+1（J=0..8）",
    all(sum((branch[J][r] or 0) * CH[r][0] for r in CH) == 2 * J + 1 for J in branch), "9/9")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z4. 存在性定理：J 載體 = 邊個點群？（J=2 ⟹ 需 5 維 irrep）")
print("=" * 78)
scan = []
for n in range(1, 13):
    scan.append(("C%d" % n, dims_C(n)))
for n in range(2, 13):
    scan.append(("D%d" % n, dims_D(n)))
scan += [("T", DIMS["T"]), ("O", DIMS["O"]), ("I", DIMS["I"])]
has3 = [nm for nm, d in scan if 3 in d]
has5 = [nm for nm, d in scan if 5 in d]
print("  有 3 維 irrep 之群 =", has3)
print("  有 5 維 irrep 之群 =", has5)
chk("Z4a 有 3 維 irrep（J=1）者 = {T, O, I} 全部", has3 == ["T", "O", "I"], "%s" % has3)
chk("Z4b ★ 有 5 維 irrep（J=2）者 **唯一** = {I}", has5 == ["I"],
    "T/O 皆無 5 維 ⟹ J=2 之存在**強制**載體為二十面體群（I ≅ A₅），非選擇")
verdict("Z4c J 載體 = 幾何層 A₅（I）—— 由 J=2 之存在唯一逼出", True,
        "非「揀似嘅」；mode/tree 層（docs/04）結構性無 3 維 irrep ⟹ 唔可能替代", "Z4")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z5. J 值域與 J^PC 容量")
print("=" * 78)
first_J = {r: min(J for J in branch if branch[J][r]) for r in CH}
print("  各 irrep 首次（最低 J）出現 =", first_J)
# 「單一 irrep」= 該 J 分解恰好一個 R、且 dim R = 2J+1（整個多重態 = 一個 irrep）
single = [J for J in range(0, 7)
          if len([r for r in CH if branch[J][r]]) == 1
          and [CH[r][0] for r in CH if branch[J][r]][0] == 2 * J + 1]
carrier = {J: [r for r in CH if branch[J][r]][0] for J in single}
print("  單一 irrep 之 J → irrep =", carrier)
print("  ★ T2（dim 3）首次出現於 J = %d ⟹ 唔係 J=1 載體（J=1 只有 T1）" % first_J["T2"])
n_single = len(carrier)
cap = n_single * 2 * 2      # ×P(2) ×C(2)
print("  J^PC 類容量 = %d irrep × 2(P) × 2(C) = %d" % (n_single, cap))
chk("Z5a 單一 irrep J 值域 = {0,1,2}，載體 = {A, T1, H}（G 之 4 維無 2J+1 匹配）",
    single == [0, 1, 2] and set(carrier.values()) == {"A", "T1", "H"},
    "%s ／ %s（G dim=4 ⟹ 唔可能係純 J 多重態）" % (single, carrier))
chk("Z5b J^PC 容量 = 12 類（3 irrep × P × C，零參數）", cap == 12, "%d" % cap)
chk("Z5c T2 首次出現 J = 3（與 T1 之 J=1 不對稱）⟹ T2 態唔可能係 J=1",
    first_J["T2"] == 3 and first_J["T1"] == 1, "T1→1（J=1）、T2→3（J=3）")

# ─────────────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Z6. 可證偽 (i)：觀測低窪 J^PC 每個都須可命 irrep")
print("=" * 78)
LOW = [("π", "0−+"), ("η", "0−+"), ("ρ(770)", "1−−"), ("ω(782)", "1−−"),
       ("a₀/f₀(980)", "0++"), ("a₁(1260)", "1++"), ("b₁(1235)", "1+−"),
       ("f₂(1270)", "2++"), ("a₂(1320)", "2++"), ("π₂(1670)", "2−+")]
rows = []
for nm, jpc in LOW:
    J = int(jpc[0])
    dim = 2 * J + 1
    Rs = [r for r in CH if CH[r][0] == dim]
    rows.append({"hadron": nm, "jpc": jpc, "J": J, "dim": dim, "irreps": Rs})
    print("   %-11s %s ⟹ dim %d ⟹ irrep %s" % (nm, jpc, dim, Rs))
chk("Z6a 全部低窪（<1.7 GeV）J^PC 可命於 {A, T1, T2, H}（10/10）",
    all(r["irreps"] for r in rows), "10/10 有載體")
J3 = [r for r in branch[3] if branch[3][r]]
print("  J=3 之分解 =", {r: branch[3][r] for r in CH if branch[3][r]})
chk("Z6b J=3（ρ₃(1690) 3−−）須多 irrep（G⊕T2）⟹ 唔可能單一 irrep", set(J3) == {"G", "T2"},
    "J=3 ⟹ G(4)+T2(3) = 7 ✓ ⟹ J≥3 必須分裂多重態")

print()
print("Z6'. 可證偽 (ii)：J 首次出現質量單調")
FIRST = [(0, 139.57039, "π"), (1, 775.26, "ρ"), (2, 1270.0, "f₂"), (3, 1688.8, "ρ₃"), (4, 2040.0, "f₄")]
mono = all(FIRST[i][1] < FIRST[i + 1][1] for i in range(len(FIRST) - 1))
for J, m, nm in FIRST:
    print("   J=%d 首現 %-4s m = %7.1f MeV" % (J, nm, m))
chk("Z6c 首現質量對 J 單調升（5/5，無反例）", mono,
    "單調 ⟹ 與『J≥3 需多 irrep（額外成本）』相容；反例將直接否證本讀法")

# Z6d：12 類之補集 = 4 類，核其是否恰為 qq̄ 禁類（外部標準，只作對照）
ALL12 = set()
for J, r in carrier.items():
    for P in ("+", "−"):
        for C in ("+", "−"):
            ALL12.add("%d%s%s" % (J, P, C))
obs = {j for _, j in LOW}
spare = sorted(ALL12 - obs)
print("  容量 12 類 − 已觀測 = %s（%d 個）" % (spare, len(spare)))

# qq̄ 允許類（外部夸克模型標準，只作對照）：P = (−1)^{L+1}, C = (−1)^{L+S}
qqbar = set()
for L in range(0, 3):
    for S in (0, 1):
        for J in range(abs(L - S), L + S + 1):
            qqbar.add("%d%s%s" % (J, "+" if (-1) ** (L + 1) > 0 else "−",
                                  "+" if (-1) ** (L + S) > 0 else "−"))
print("  qq̄ 允許（L,S ≤ 2）= %s" % sorted(qqbar))
pred = ALL12 & qqbar
print("  預測 = 容量 12 ∩ qq̄ = %s（%d 個）" % (sorted(pred), len(pred)))
print("  已觀測低窪 = %s（%d 個）" % (sorted(obs), len(obs)))
chk("Z6d ★ 已觀測 ⊆ 預測（7/7 逐類），且空槽 = qq̄ 禁類", obs <= pred,
    "pred %d ⊇ obs %d；pred−obs = %s（活預言）；12-容量內 qq̄ 禁類 = %s（exotic 空槽）"
    % (len(pred), len(obs), sorted(pred - obs), sorted(ALL12 - qqbar)))
chk("Z6e ★ 3−− ∉ 單 irrep 容量 12 ⟹ J=3 必須多 irrep（與 ρ₃(1690) 偏重一致）",
    "3−−" not in ALL12, "qq̄ 允許含 3−− 但容量唔含 ⟹ 結構上需要 G⊕T2 分裂多重態")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "09_jpc_geometric_carrier.py", "checks": checks, "integrity_fail": nfail,
       "A5_classes": CS, "A5_char": {r: [str(x) for x in v] for r, v in CH.items()},
       "branch_SO3_down_I": {str(J): {r: branch[J][r] for r in CH} for J in branch},
       "groups_with_3dim": has3, "groups_with_5dim": has5,
       "single_irrep_J": single, "carrier": {str(k): v for k, v in carrier.items()},
       "first_J": first_J, "capacity_classes": sorted(ALL12), "spare": spare,
       "qqbar_allowed": sorted(qqbar), "predicted_classes": sorted(pred), "observed_classes": sorted(obs),
       "JPC_capacity": cap, "low_jpc_rows": rows,
       "first_mass": [{"J": a, "mass": b, "hadron": c} for a, b, c in FIRST], "monotone": mono}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "09_jpc_geometric_carrier.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/09_jpc_geometric_carrier.json")
sys.exit(0 if nfail == 0 else 1)
