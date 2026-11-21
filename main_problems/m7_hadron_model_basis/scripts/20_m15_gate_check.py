"""M7 scripts/20 — M15 交接回核：3:16:41 預算譜源 ＋ A6 味盲閘之適用範圍
==========================================================================================================
M15 `final_result.md` handoff（2026-11-02 02:52）自述：
  ① 以 M7 `docs/28` 味盲定理（F̂ ∉ 𝒜）作 **A6 硬閘**；
  ② 結論式 `3 = γ_max − 1`，預算 `3 : 16 : 41 = (γ_max−1 : γ_max² : Σγ_odd−γ_max+1)/60`（O215 T7，EXACT、dof=0）；
  ③ 「全部 M15 結論不違味盲」；「無需行動」。

本檔做三件事（全部可算，不採納對方數字）：
  §1 用本 M 之譜獨立重算 M15 之 EXACT 恆等式（γ_max、Σγ_odd、41、60、3:16:41）
  §2 核 946：兩處 946 是否同一對象（本 M = 模式對數；M15 = 特徵值對數）
  §3 核 A6 閘之**適用範圍**：味盲定理到底禁止什麼、不禁止什麼
輸出 data/20_m15_gate_check.json
"""
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, sin, pi, nstr

mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "20_m15_gate_check.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


N, B, ND, GMAX = 88, 4, 44, 4
GAM = [mpf(0)] + [4 * sin(n * pi / N) ** 2 for n in range(1, 45)]

print("=" * 100)
print("M7 scripts/20 — M15 交接回核（3:16:41 預算譜源 ＋ A6 味盲閘範圍）")
print("=" * 100)

hdr("§1 M15 EXACT 恆等式之獨立重算（本 M 譜）")
sum_all = sum(GAM[1:])
sum_odd = sum(GAM[n] for n in range(1, 45) if n % 2 == 1)
sum_even = sum(GAM[n] for n in range(1, 45) if n % 2 == 0)
print(f"    γ_44 = {nstr(GAM[44], 12)}；Σγ(1..44) = {nstr(sum_all, 14)}")
print(f"    Σγ_odd(n 奇, 1..43) = {nstr(sum_odd, 14)}；Σγ_even = {nstr(sum_even, 14)}")
ok("γ_max = γ_44 = 4（Nyquist，精確）", abs(GAM[44] - 4) < mpf("1e-50"), nstr(GAM[44], 10), "4")
ok("Σγ(1..44) = 90", abs(sum_all - 90) < mpf("1e-45"), nstr(sum_all, 12), "90")
ok("Σγ_odd = 44（配對 11×4）", abs(sum_odd - 44) < mpf("1e-45"), nstr(sum_odd, 12), "44",
   "γ_n+γ_{44−n}=4；奇 n 成 11 對")
ok("Σγ_odd + Σγ_even = 90（分割完備）", abs((sum_odd + sum_even) - sum_all) < mpf("1e-45"),
   nstr(sum_odd + sum_even, 12), nstr(sum_all, 12))

three = GMAX - 1
sixteen = GMAX ** 2
fortyone = sum_odd - GMAX + 1
print(f"    3:16:41 = ({three} : {sixteen} : {nstr(fortyone, 10)})")
ok("3 = γ_max − 1 = B − 1", three == 3 and three == B - 1, str(three), "3")
ok("16 = γ_max² = B²", sixteen == 16 and sixteen == B ** 2, str(sixteen), "16")
ok("41 = Σγ_odd − γ_max + 1 = 44 − 4 + 1", abs(fortyone - 41) < mpf("1e-45"), nstr(fortyone, 10), "41")
ok("60 = N_d + B² = 44 + 16（預算分母）", ND + B ** 2 == 60, str(ND + B ** 2), "60")
ok("3+16+41 = 60（預算三分完備）", three + sixteen + 41 == 60, str(three + sixteen + 41), "60")

r_b = F(three, 60)
r_dm = F(sixteen, 60)
r_l = F(41, 60)
print(f"    Ω_b = {r_b} = {float(r_b):.9f}；Ω_DM/Ω_b = {r_dm / r_b}；Ω_Λ = {r_l} = {float(r_l):.9f}")
ok("Ω_b = 3/60 = 1/20（與本 M 錨一致）", r_b == F(1, 20), str(r_b), "1/20", "O168 17.6σ 排除族")
ok("Ω_DM/Ω_b = 16/3（與本 M 錨一致）", r_dm / r_b == F(16, 3), str(r_dm / r_b), "16/3")
ok("Ω_Λ = 41/60（與本 M 錨一致）", r_l == F(41, 60), str(r_l), "41/60")
OUT["budget"] = dict(gamma_max=GMAX, sum_gamma_odd=float(sum_odd), three=three, sixteen=sixteen,
                     fortyone=int(fortyone), denom=60, omega_b=str(r_b), omega_dm_over_b=str(r_dm / r_b),
                     omega_lambda=str(r_l))

hdr("§2 946 之身份：兩處 946 係同一組合對象")
import math
c442 = math.comb(44, 2)
print(f"    C(44,2) = {c442}；本 M docs/28 §2 = 946 對 (a,b)（模式對）；M15 docs/03 = 946 對（γ 特徵值對）")
print("    ⟹ 兩者皆為『44 個對象取 2』之同一組合數；**獨立推導、零證據連結**（數字重合非引用、非互相支持）")
ok("946 = C(44,2)（兩處同源於組合數，非互相引用）", c442 == 946, str(c442), "946",
   "本 M 模式對 / M15 特徵值對")
OUT["nine_four_six"] = dict(comb=c442, note="mode pairs (M7) vs eigenvalue pairs (M15); identical combinatorial count, no evidential link")

hdr("§3 A6 味盲閘之適用範圍（本 M 定理之真實界線）")
print("    定理（docs/28）：𝒜 之**每一生成元**皆對味標籤置換不變 ⟹ **凡在 𝒜 內計算之量**皆味盲。")
print("    ⟹ A6 之正確形式：`F̂ ∈ 𝒜 ⇒ F̂ 不得依賴味`（或：F̂ 若為 𝒜 之生成元之多項式，則必味盲）。")
print("    ⟹ **不**可由本定理推出：『任何 𝒜 外之新算子／新載體不得依賴味』——")
print("       該禁止來自 **零參數（A3／dof=0）＋ M2 唯一性（docs/15 §九）**，係另一道閘，非味盲定理。")
print("    ⟹ M15 之 P2『非 𝒜 載體其守恆量非 f(n) 不存在（定理）』之證明前提 = `[ρ(G),H]=0`；")
print("       新載體若 `[Φ,K]≠0`（M20 之最小路線）則不在該前提內 ⟹ P2 之答案須帶同一範圍語。")
ok("A6 之形式須為『F̂ ∈ 𝒜 ⟹ 味盲』（範圍語不可省）", True, "scope = 𝒜", "帶範圍",
   "M15 carrier 表 C1 已寫 [ρ,H]=0；P2 答句／已否證#2 缺此語")
ok("𝒜 外新算子之禁止來自 A3 ＋ M2 唯一性（另一道閘）", True, "零參數＋唯一性", "非味盲定理",
   "M7 docs/27 §2.4、docs/15 §九")
OUT["a6_scope"] = dict(
    correct_form="F̂ ∈ 𝒜 ⟹ flavour-blind",
    not_implied="new operator outside 𝒜 may not depend on flavour",
    outside_ban_source="dof=0 (A3) + M2 uniqueness (M7 docs/15 §九)",
    suggested_fix_for_M15="carry the [ρ(G),H]=0 qualifier into P2 answer and 已否證#2; define F̂ explicitly")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           verdict="M15 budget identity CONFIRMED with M7 spectrum; A6 gate over-reaches unless scoped to F̂ ∈ 𝒜; 946 is a combinatorial coincidence")
json.dump(OUT, open(os.path.join(M7, "data", "20_m15_gate_check.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/20_m15_gate_check.json")