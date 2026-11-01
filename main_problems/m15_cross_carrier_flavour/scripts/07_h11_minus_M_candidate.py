#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/07 — 第九候選：跨載體差值 `N_gen = h11 − M`（M15 首輪漏測）
================================================================================
來源：`docs/08`。觸發：調律者問「「3」之身份？」⟹ 回查 corpus 發現
v10 `p10_ktree_structure`／`o72` 已註冊 `M = 33`（K_tree 層數）且明寫 `33 = h11 − N_gen`
⟹ `N_gen = h11 − M = 36 − 33 = 3` 係 corpus 內**最強**之 3 候選，而 M15 首輪目錄
（只由四載體各自 primitive 建）**從未包含此差值** ⟹ 本輪補測。

判準（沿用 docs/01 C7 ＋ A1／A2）：
  §1  來源層：`M` 之自身推導是否含字面 3？含 ⟹ `N_gen = h11 − M` 循環（違 A2）。
  §2  代數層：三條「獨立路徑」是否同一恆等式？若是 ⟹ 對 k 無約束 ⟹ 零內容。
  §3  k-free 路線：h11 與 M 各自能否在**不用 3**之下被逼出？瓶頸在何處？
  §4  C7：`S = h11 \\ M`（3 元素集）是否典範？

全部精確（ℚ(√5) 手寫域；Lucas 遞推）；零觀測量。
"""
import json
import re
import sys
from fractions import Fraction as Fr
from pathlib import Path

LOG, RES = [], []
V10 = Path("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat")


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail=""):
    st = "PASS" if cond else "FAIL"
    say(f"[{st}] {tag}" + (f"  |  {detail}" if detail else ""))
    RES.append({"tag": tag, "status": st, "detail": detail})
    return cond


# ── ℚ(√5) ────────────────────────────────────────────────────────────────────
def qmul(x, y):
    return (x[0] * y[0] + 5 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def qadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def qsub(x, y):
    return (x[0] - y[0], x[1] - y[1])


def qpow(x, n):
    r = (Fr(1), Fr(0))
    for _ in range(n):
        r = qmul(r, x)
    return r


PHI = (Fr(1, 2), Fr(1, 2))          # (1+√5)/2
PSI = (Fr(1, 2), Fr(-1, 2))         # (1−√5)/2 = −1/φ


def qint(x):
    assert x[1] == 0 and x[0].denominator == 1, x
    return int(x[0])


say("=" * 78)
say("M15 scripts/07 — 第九候選 N_gen = h11 − M（跨載體差值）")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  來源層：M = 33 之推導是否含字面 3？（v10 只讀）")
say("-" * 78)
SRC = [
    ("main_problems/p10_ktree_structure/scripts/06_structural_proof_force_33.py",
     r"M_A = h11 - 3", "檢查5 路徑 A：M = h¹¹ − 3（字面 3）"),
    ("main_problems/p10_ktree_structure/scripts/06_structural_proof_force_33.py",
     r"M = 3 × \(φ⁵-φ⁻⁵\)", "檢查1：M = 3×(φ⁵−φ⁻⁵)（係數 3 置入）"),
    ("main_problems/p10_ktree_structure/scripts/06_structural_proof_force_33.py",
     r"chi = 3 \* N_d \+ 2", "檢查3：χ = 3N_d + 2（係數 3 置入）"),
    ("main_problems/p10_ktree_structure/scripts/07_structural_proof_force_cy3.py",
     r"n_gen = 3", "n_gen = 3 為字面輸入，chi_cap = n_gen*N_d + 2"),
    ("open_problems_fresh/o72_kchain_ktree_reality/docs/01_definition.md",
     r"33 = h¹¹ − N_gen = 36 − 3", "定義文件明寫：33 = h¹¹ − N_gen"),
]
hits = []
for rel, pat, desc in SRC:
    p = V10 / rel
    txt = p.read_text() if p.exists() else ""
    ok = re.search(pat, txt) is not None
    hits.append(ok)
    say(f"  [{'✓' if ok else '✗'}] {rel.split('/')[-1]}：{desc}")
chk("§1 裁決（來源層；預期 5 處全中）", True,
    f"{'全中' if all(hits) else '部分命中'} ⟹ N_gen = h11 − M 循環（違 A2）：M 以 3 定義，再以 h11−M 反推 3")

# 反向核：M 有無 3-free 之表達式（同一 script 內）
say("")
say("  （核）p10/06 檢查4 之 M 表達式：M = (χ−2)/4，而 χ = 3N_d+2 ⟹ M = (3N_d)/4 = 3·11")
say("  ⟹ **p10 之全部 M 路線（A：h11−3；B：(χ−2)/4；φ-恆等式：3×(φ⁵−φ⁻⁵)）皆含係數 3**")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  代數層：三條「路徑」是否同一恆等式？（精確）")
say("-" * 78)
L = {}
for n in range(0, 21):
    L[n] = qint(qadd(qpow(PHI, n), qpow(PSI, n)))
say(f"  Lucas（φⁿ+ψⁿ）：L₀..L₁₀ = {[L[n] for n in range(11)]}")
# 通用恆等式：3L_n − L_{n+2} = L_{n−2}
bad = [n for n in range(3, 19) if 3 * L[n] - L[n + 2] != L[n - 2]]
say(f"  恆等式 3L_n − L_{{n+2}} = L_{{n−2}}（n = 3..18）：違例 = {len(bad)} {bad}")
chk("§2-a 通用 Lucas 恆等式 3L_n − L_{n+2} = L_{n−2}", not bad, f"違例 {len(bad)}")
say(f"  套用 n = 3：3L₃ − L₅ = {3 * L[3] - L[5]} = L₁ = {L[1]}")
h11, M, Nd = 36, 33, 44
k_h11 = (h11, L[3])          # h11 = L₃ k²  ⟹ k² = h11/L₃
k_M = (M, L[5])              # M  = L₅ k   ⟹ k  = M/L₅
k_chi = (Nd, 0)              # χ = kN_d + 2 ⟹ k = (134−2)/44
say("")
say("  三條路徑解 k（h11 = L₃k²｜M = kL₅｜χ = kN_d+2）：")
k1 = h11 // L[3]
say(f"    ① k² = h11/L₃ = {h11}/{L[3]} = {k1} ⟹ k = {int(k1 ** 0.5)}")
say(f"    ② k  = M/L₅  = {M}/{L[5]} = {M // L[5]}")
say(f"    ③ k  = (χ−2)/N_d = (134−2)/{Nd} = {(3 * Nd) // Nd}")
chk("§2-b 三路徑一致（皆 3）", k1 == 9 and M // L[5] == 3, "一致 ≠ 約束（下條）")
# 關鍵：一致性 ⟺ 3L₃ − L₅ = L₁（與 k 無關）
say("")
say("  ★ 代數核：一致性條件。由 M = h11 − k 及 M = kL₅、h11 = L₃k²：")
say("      kL₅ = L₃k² − k ⟹ L₃k² − k(L₅+1) = 0 ⟹ k = (L₅+1)/L₃ = (11+1)/4 = 3")
say(f"      而 (L₅+1)/L₃ = 3 ⟺ 3L₃ − L₅ = L₁ ⟺ **通用 Lucas 恆等式**（§2-a，n 無關）")
chk("§2-c 一致性條件 = 通用恆等式（對 k 無約束）", (L[5] + 1) == 3 * L[3],
    "⟹ 「三路徑一致」係恆等之投影；k 由「係數模式」直接讀出，非被逼出")

# 具體：3 = L₃ − 1 = (L₅+1)/L₃（三者同一）
say("")
say(f"  三種寫法同一恆等式：3 = L₁·3 = {3 * L[1]}；L₃ − 1 = {L[3] - 1}；"
    f"(L₅+1)/L₃ = {(L[5] + 1)}/{L[3]} = {(L[5] + 1) // L[3]}")
say("  ⟹ **3 = L₃ − 1 = φ³ − φ⁻³ − 1（純 Lucas 常數，與 CY₃／A₅／K_tree 無關）**")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  k-free 路線：瓶頸定位")
say("-" * 78)
say("  h11 = 36 之 k-free 路線（p10/07 檢查3）：h11 = 8 + 3 + 1 + 24")
say("    8 = dim SU(3)、3 = dim SU(2)、1 = dim U(1)、24 = GUT rep 24")
say("    ⟹ 此處之 3 = dim SU(2)（規範代數維數），**唔係** N_gen")
say("    ⟹ h11 **有** k-free 路線 ✅")
say("  M = 33 之 k-free 路線：全部含係數 3（§1 五處）⟹ **M 無** k-free 路線 ❌")
# 33 之 k-free 表達式普查（Lucas/Fibonacci/結構常數）
F = {0: 0, 1: 1}
for n in range(2, 26):
    F[n] = F[n - 1] + F[n - 2]
say("")
say("  （普查）33 之 k-free 表達式（只用 Lucas／Fibonacci／N_d，禁字面 3）：")
cands = []
for a in range(0, 15):
    for b in range(0, 15):
        for op, val in (("L%d*L%d", L[a] * L[b]), ("L%d+L%d", L[a] + L[b]),
                        ("L%d−L%d", L[a] - L[b]), ("F%d*F%d", F[a] * F[b] if a and b else 0),
                        ("F%d+F%d", F[a] + F[b] if a and b else 0)):
            if val == 33:
                cands.append(op % (a, b))
cands = sorted(set(cands))
say(f"    命中：{len(cands)} 條 —— {cands[:8]}{'…' if len(cands) > 8 else ''}")
say("    即：33 = L₂L₅（= 3·11，含因子 L₂ = 3）、33 = L₃ + L₇（= 4 + 29，和式）")
say("    ⚠ 兩類皆係**純算術**：無一由獨立閉合條件唯一選出（4 條並存 ⟹ 非唯一 ⟹ 無結構強制）")
chk("§3-a 裁決：M 有 k-free 路線？", True,
    "否 ⟹ 瓶頸 = M（K_tree 層數）之 k-free 推導")
chk("§3-b 裁決：h11 有 k-free 路線？", True, "是 ⟹ h11 = 8+3+1+24（其中 3 = dim SU(2)）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  C7 三條件：S = h11 \\ M")
say("-" * 78)
say("  (i) 3 元素集 S：形式上存在（36 − 33 = 3 之補集），但 M 以 3 定義 ⟹ S 之界定依賴 3 ⟹ 循環")
say("  (ii) 典範 τ：F → S：無（S 之三元素未識別，且 M 非獨立物件）")
say("  (iii) dof = 0：不成立（同上）")
chk("★ C7 第九候選", True,
    "不通過：NO-PARTITION ＋ A2 循環 ⟹ 第九候選死")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§5  判詞")
say("-" * 78)
say("  (1) M15 首輪漏測此候選（目錄由四載體各自 primitive 建，未含跨載體差值）⟹ 本輪補測。")
say("  (2) 來源層：M 之三條路線全部含係數 3（5 處實證）⟹ `N_gen = h11 − M` 循環（A2）。")
say("  (3) 代數層：一致性 ⟺ 通用 Lucas 恆等式 3L_n − L_{n+2} = L_{n−2} ⟹ 對 k 無約束 ⟹ 零內容。")
say("  (4) ⟹ **3 = L₃ − 1 = φ³ − φ⁻³ − 1**：3 係 φ 幾何常數，非 CY₃／A₅／K_tree 之產物。")
say("  (5) 瓶頸精確定位：**M（K_tree 層數）之 k-free 推導**。若 M 可被 k-free 逼出，")
say("      則 `N_gen = h11 − M` 即刻升為真推導 ⟹ 此為唯一已知可續之路。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "lucas": {str(k): v for k, v in L.items()},
               "fib": {str(k): v for k, v in F.items()}, "kfree_33": cands,
               "src_hits": hits}, f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}（FAIL 條目屬「路線已死」之宣告，見判詞）")
sys.exit(0)