#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/01 — P1／P2：歸約定理之形式地位 ＋ 跨載體第二標籤之不可能性
================================================================================
判準來源：docs/01_frozen_criteria.md §三（C1–C6）、§四（盲測）。跑前凍結。
慣例：/88 = 環（44 模，Σγ=90）；/90 = Dirichlet（44 站，Σγ=88=2N_d）。不得混用。

§1  C1  恆真審計：O137 §一 清單 ①–⑤ 含 `n` 本身（①）且含對 n 單射之量（④ γ_n）
         ⟹ 「任何幾何供給之量皆 f(n)」對「由 n 算出之量」恆真 ⟹ 零鑑別力。
         量化：γ_n 對 n=1..44 之最小間隙 > 0 ⟹ ④ 單獨已經決定 n ⟹ 清單冗餘。
§2  C4  A₅ 完美性精確核（[A₅,A₅] = A₅，階 60）⟹ 任何與 H 交換之 A₅-作用必然平凡。
§3  C5  簡併逃逸之精確否證：γ_{n+1} − γ_n = 4 sin((2n+1)π/90) sin(π/90) > 0（n=1..43），
        且 946 對全部 γ 互異（兩載體同做）。
§4  C6  ℤ_m 分級之 dof：
        (a) 反射等變分級之自由度（精確計數）；
        (b) 算術分級 n ↦ n mod m 之族：類大小、是否含「同族子鏈」類、反射等變性；
        (c) 判詞：m 自由 ⟹ dof > 0 ⟹ 「3 類 mod 3」屬循環（把答案當前提）。
§5  盲測 B1 之前置（A₅-rep 於模式空間之存在性）：A₅ 之 44 維表示分解計數（存在 ≠ 可用）。

性質：確定性、零觀測量、N_trials = 0（§4(c) 之 m 掃描為窮舉，非統計試驗）。
"""
import itertools
import json
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 60
LOG = []
RES = {}


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail="", warn=False):
    ok = bool(cond)
    st = "PASS" if ok else ("WARN" if warn else "FAIL")
    line = f"[{st}] {tag}" + (f"  |  {detail}" if detail else "")
    say(line)
    RES[tag] = {"status": st, "detail": detail}
    return ok


# ─────────────────────────────────────────────────────────────────────────────
# 慣例對象
# ─────────────────────────────────────────────────────────────────────────────
ND = 44          # 站／模數
X90 = 90         # Dirichlet 分母
X88 = 88         # 環分母


def gamma90(n):
    """Dirichlet /90：γ_n = 4 sin²(nπ/90)，n = 1..44，Σγ = 88"""
    return 4 * mp.sin(mp.pi * n / X90) ** 2


def gamma88(n):
    """環 /88：γ_n = 4 sin²(nπ/88)，n = 1..44，Σγ = 90"""
    return 4 * mp.sin(mp.pi * n / X88) ** 2


say("=" * 78)
say("M15 scripts/01 — P1／P2 歸約定理形式地位 ＋ 第二標籤不可能性")
say(f"慣例宣告：/90 Dirichlet（44 站，Σγ=88）｜/88 環（44 模，Σγ=90）｜dps={mp.mp.dps}")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  C1 恆真審計 — O137 §一 清單之形式地位")
say("-" * 78)

# O137 §一 清單（逐字重列，跨版本只讀）
LIST = {
    "①": "n 本身（整數標籤）",
    "②": "parity(n) = n mod 2",
    "③": "block(n) = ceil(rank/3)",
    "④": "γ_n = 4 sin²(nπ/88)（K_chain 特徵值）",
    "⑤": "#{l : h(l) ≠ 0, l ≤ n, l ≡ n (mod 2)}",
}
say("O137 §一 清單（v10，只讀）：")
for k, v in LIST.items():
    say(f"  {k} {v}")

# γ_n（/88）之最小間隙 ⟹ ④ 單獨對 n 單射
g88 = [gamma88(n) for n in range(1, ND + 1)]
gaps = [g88[n] - g88[n - 1] for n in range(1, ND)]
min_gap = min(gaps)
max_gap = max(gaps)
say("")
say(f"  γ（/88）最小相鄰間隙 = {mp.nstr(min_gap, 12)}（n = {int(gaps.index(min_gap)) + 2}→{int(gaps.index(min_gap)) + 1}）")
say(f"  γ（/88）最大相鄰間隙 = {mp.nstr(max_gap, 12)}")
chk("C1-a 清單含『n 本身』(①)",
    True,
    "⟹ 任何由 n 定義之量皆 f(n) ⟹ 定理對『由 n 算出之量』恆真")
chk("C1-b ④ γ_n 單獨對 n 單射（最小間隙 > 0）",
    min_gap > 0,
    f"min gap = {mp.nstr(min_gap, 8)} ⟹ 清單中 ②／③／⑤ 對『決定 n』皆冗餘")
say("")
say("  ★ 結論 C1：定理之邏輯形式 = ∀D(由 n 定義) : D = f_D(n) —— **恆真**。")
say("    其鑑別力 = 0（與 M14 docs/22 E3「恆真判 ⟹ 零檢定力」同型）。")
say("    O137 §一 之真正內容唔係『幾何冇第二自由度』，而係：")
say("      (i) 幾何之模式指標 = 鏈之模指標（同一性），(ii) m = f(n) 之 f 低複雜度不存在（§二 屏障）。")
say("    (ii) 係維數論證（9 點／9 參數 = 飽和），唔係幾何必然 ⟹ 見 docs/02。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  C4 A₅ 完美性 — 非阿貝爾群與 H 交換之作用必然平凡")
say("-" * 78)

def perm_mul(a, b):
    """(a·b)(i) = a(b(i))"""
    return tuple(a[b[i]] for i in range(5))


def perm_inv(a):
    inv = [0] * 5
    for i, v in enumerate(a):
        inv[v] = i
    return tuple(inv)


def _sign(p):
    s = 1
    for i in range(5):
        for j in range(i + 1, 5):
            if p[i] > p[j]:
                s = -s
    return s


A5 = [p for p in itertools.permutations(range(5)) if _sign(p) == 1]
assert len(A5) == 60, len(A5)
A5set = set(A5)

# 交換子子群生成集：閉包下逐層擴張（有限群，直接 BFS 到不動點）
gens = []
for a in A5:
    for b in A5:
        gens.append(perm_mul(perm_mul(a, b), perm_mul(perm_inv(a), perm_inv(b))))
gens = list(set(gens))
comm = {tuple(range(5))}
frontier = {tuple(range(5))}
while frontier:
    new = set()
    for x in frontier:
        for g in gens:
            y = perm_mul(x, g)
            if y not in comm:
                comm.add(y)
                new.add(y)
    frontier = new
say(f"  |A₅| = {len(A5)}；交換子生成元 {len(gens)} 個；⟨[A₅,A₅]⟩ 之階 = {len(comm)}")
chk("C4-a [A₅,A₅] = A₅（完美群）", len(comm) == 60, f"階 {len(comm)}/60")

# 非阿貝爾核（存在 a,b: ab ≠ ba）
nonab = None
for a in A5:
    for b in A5:
        if perm_mul(a, b) != perm_mul(b, a):
            nonab = (a, b)
            break
    if nonab:
        break
chk("C4-b A₅ 非阿貝爾", nonab is not None, "存在 ab ≠ ba 之對")

say("")
say("  定理（本輪）：設 G 以表示 ρ 作用於模空間 V（dim 44）且 [ρ(G), H] = 0。")
say("    則 ρ(G) ⊆ 𝒞(H)；𝒞(H) 對角（F3）⟹ 阿貝爾 ⟹ ρ(G) 阿貝爾。")
say("    A₅ 完美且非阿貝爾 ⟹ 任何同態 ρ: A₅ → Aut(V) 滿足上條件者必平凡。")
say("  ⟹ **A₅ 幾何不可能在模式空間上供給任何非平凡標籤。**")
chk("C4-c 推論：A₅ 作用平凡 ⟹ 無 A₅-標籤", len(comm) == 60 and nonab is not None,
    "完美 + 非阿貝爾 + 交換子群 = 自身")

# A₅ 之 44 維表示存在性（存在 ≠ 可用；反面核）
# 44 = 5a+4b+3c+3d+e over dims {5,4,3,3,1}
sols = []
for a in range(0, 9):
    for b in range(0, 12):
        for c in range(0, 15):
            for d in range(0, 15):
                r = 44 - 5 * a - 4 * b - 3 * c - 3 * d
                if r >= 0:
                    sols.append((a, b, c, d, r))
say(f"  （核）A₅ 之 44 維表示分解數 = {len(sols)}（例如 5a+4b+3c+3d+1e = 44）⟹ 表示存在，")
say("       但由上述定理，任何此類表示皆不可與 H 交換 ⟹ 『存在』與『可用』分離。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  C5 簡併逃逸之否證")
say("-" * 78)

for name, g, x, xlabel in (("/90", gamma90, X90, "Σγ=88"), ("/88", gamma88, X88, "Σγ=90")):
    vals = [g(n) for n in range(1, ND + 1)]
    # 精確公式核：γ_{n+1} − γ_n = 4 sin((2n+1)π/x) sin(π/x)
    maxdev = mp.mpf(0)
    for n in range(1, ND):
        lhs = g(n + 1) - g(n)
        rhs = 4 * mp.sin((2 * n + 1) * mp.pi / x) * mp.sin(mp.pi / x)
        maxdev = max(maxdev, abs(lhs - rhs))
    mind = min(vals[n] - vals[n - 1] for n in range(1, ND))
    # 946 對互異
    npair = ND * (ND - 1) // 2
    coll = 0
    for i in range(ND):
        for j in range(i + 1, ND):
            if abs(vals[i] - vals[j]) < mp.mpf("1e-50"):
                coll += 1
    say(f"  {name}（{xlabel}）：Δγ 閉式 max|dev| = {mp.nstr(maxdev, 6)}，最小相鄰差 = {mp.nstr(mind, 12)}，"
        f"{npair} 對中簡併 = {coll}")
    chk(f"C5 γ 嚴格遞增且無簡併 {name}", maxdev < mp.mpf("1e-50") and coll == 0 and mind > 0,
        f"閉式 dev {mp.nstr(maxdev, 4)}；簡併 0/{npair}")

say("")
say("  ⟹ 無簡併 ⟹ 冇退化子空間可容非平凡標籤 ⟹ P2 之『非 𝒜 載體』路線之唯一逃逸口關閉。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  C6 ℤ_m 分級之 dof（(a) 等變分級計數；(b) 算術分級族；(c) 循環性判詞）")
say("-" * 78)

# (a) 反射 R: n ↦ 45−n（/90 之補償模，flip 宇稱）；分級 g: {1..44} → ℤ_m
# 等變要求：g(45−n) = −g(n)（mod m）。
# 配對 {n, 45−n} 把 {1..44} 分成 22 對（無固定點：45−n = n ⟹ n=22.5 非整數）
pairs = []
seen = set()
for n in range(1, ND + 1):
    if n in seen:
        continue
    m_ = 45 - n
    assert 1 <= m_ <= ND, (n, m_)
    pairs.append((min(n, m_), max(n, m_)))
    seen.add(n)
    seen.add(m_)
say(f"  反射配對數 = {len(pairs)}（每對 2 元素；45−n = n 無整數解 ⟹ 無固定點）")


def count_equiv_gradings(m, sign):
    """每對 (a, 45−a)：選 g(a) ∈ ℤ_m 自由，g(45−a) = sign·g(a)。自由位 = 22 ⟹ m^22"""
    return m ** len(pairs)


for m in (2, 3, 4, 5):
    n_neg = count_equiv_gradings(m, -1)
    n_pos = count_equiv_gradings(m, +1)
    say(f"    m = {m}：反射反等變分級 = {n_neg} 個（= {m}^22）；等變 = {n_pos} 個")
chk("C6-a 反射等變 ℤ_3 分級之自由度 > 0",
    count_equiv_gradings(3, -1) > 1,
    f"3^22 = {3 ** 22} 個 ⟹ dof = 22 個三元選擇 ⟹ 違 A3（除非另有結構強制）")

# (b) 算術分級 g(n) = n mod m：類大小／同族子鏈／反射等變
say("")
say("  (b) 算術分級 g(n) = n mod m（canonical，dof = 0）之窮舉 m = 2..11：")
say("      m | 類大小 | 類 0 是否為同族子鏈（γ_{mk} ≡ γ_k^{(x/m)}） | 反射 g(45−n) = −g(n)?")
hdr = f"      {'m':>2} | {'類大小':<16} | {'類0同族?':<8} | {'反射反等變?':<10}"
say(hdr)
arith = {}
for m in range(2, 12):
    sizes = [0] * m
    for n in range(1, ND + 1):
        sizes[n % m] += 1
    # 類 0 之同族性：γ_{mk} vs γ_k^{(90/m)}；需 90/m 為整（m | 90）且 x′ 偶（Dirichlet 型）
    same_fam = "-"
    if X90 % m == 0:
        xp = X90 // m
        ok = True
        kmax = ND // m
        maxdev = mp.mpf(0)
        for k in range(1, kmax + 1):
            maxdev = max(maxdev, abs(gamma90(m * k) - 4 * mp.sin(mp.pi * k / xp) ** 2))
        same_fam = "是" if (maxdev < mp.mpf("1e-50") and xp % 2 == 0) else f"否(型x′={xp})"
    equiv = all((45 - n) % m == (-n) % m for n in range(1, ND + 1))
    arith[m] = {"sizes": sizes, "same_fam": same_fam, "equiv": equiv}
    say(f"      {m:>2} | {str(sizes):<16} | {same_fam:<8} | {'是' if equiv else '否':<10}")

# C10 前置：同型子鏈之 m 集合（m | 90、x′ 偶、站數 ≥ 4）
admissible = []
for m in range(2, 46):
    if X90 % m:
        continue
    xp = X90 // m
    if xp % 2:
        continue
    stations = xp // 2 - 1
    if stations >= 4:
        admissible.append((m, xp, stations))
say("")
say("  (c) 同型子鏈（m | 90 且 x′=90/m 偶 且 站數 = x′/2−1 ≥ 4）：")
for m, xp, st in admissible:
    say(f"      m = {m:>2}：x′ = {xp:>2}，站數 = {st}")
chk("C10 前置：m = 3 非唯一（≠ 盲測 B2 之否定）",
    len(admissible) >= 2,
    f"合格 m 集合 = {[a[0] for a in admissible]}（{len(admissible)} 個 ≥ 2）⟹ NON-UNIQUE")

say("")
say("  ★ 判詞 C6／C10：")
say("    (i) 任意 ℤ_3 分級（含反射等變者）= 3^22 個 ⟹ dof > 0 ⟹ 違 A3。")
say("    (ii) 算術分級 n ↦ n mod m 係 canonical（dof=0），但 m 自由：")
say("         要求『恰 3 類』即等於**假設** m = 3 ⟹ 循環論證（把答案當前題）。")
say("    (iii) 非循環之路：由同族子鏈唯一選出 m。C10 前置顯示合格 m ≥ 2 個 ⟹ 不唯一 ⟹ 路線死。")
say("    ⟹ **『3 = 模 3 類數』係循環；『3 = 同族子鏈指數』係非唯一。兩者皆不可作結構身份。**")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§5  盲測 B1 之前置：A₅-rep 於模式空間之 44 維分解（存在性計數，不含物理）")
say("-" * 78)
say(f"  分解數（dims 1,3,3,4,5 之非負整係數解）= {len(sols)} ⟹ 44 維 A₅-表示存在。")
say("  依 §2 定理，其皆不可與 H 交換 ⟹ 不構成標籤。B1 之真正測試（重數 m_ρ(n) ≥ 3）在 scripts/02。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("=" * 78)
nfail = sum(1 for v in RES.values() if v["status"] == "FAIL")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
say("  C1 歸約定理 = 恆真（零鑑別力）；C4 A₅ 作用必然平凡；C5 無簡併逃逸；")
say("  C6 ℤ_m 分級 dof > 0 或非唯一；C10 前置：m = 3 不唯一（B2 成立）。")
say("=" * 78)

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "arith_gradings": {str(k): v for k, v in arith.items()},
               "admissible_m": admissible}, f, ensure_ascii=False, indent=1)

sys.exit(1 if nfail else 0)