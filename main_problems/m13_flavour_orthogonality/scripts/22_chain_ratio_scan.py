#!/usr/bin/python3
"""22 — R1 續推：鏈譜比候選掃描（α_s/α_u 之零參數來源）＋ R2 BESIII 預測

R1：R1 之缺口係「每代因子 2」／「α_s/α_u = 4」之來源。本腳本把來源限制在
    **已註冊之鏈譜** γ_n = 4sin²(nπ/88)（n=1..44），窮舉全部比值 γ_a/γ_b，
    問：① 幾多個落觀測 1σ 窗？② 「4cos²(π/88)=γ₂/γ₁」係唯一抑或其中一個？
    ③ 窗內密度 vs 空（隨機）密度 ⟹ 有無鑑別力。
R2：給定 σ(Ξ⁰) 序列（現值、BESIII 系統、BESIII 統計），預測 B=4 vs B=5／15/4 之 z。
零 import 其他腳本；mpmath dps=40；浮點只作顯示。
"""
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data/run_22_chain_ratio.json"
out = []


def say(s=""):
    print(s)
    out.append(s)


def hdr(s):
    say()
    say("=" * 78)
    say(s)
    say("=" * 78)


N = 88
GAM = {n: 4 * mp.sin(mp.pi * n / N) ** 2 for n in range(1, N // 2 + 1)}

hdr("§0 鏈譜 γ_n = 4sin²(nπ/88)（n=1..44）")
for n in (1, 2, 3, 4, 44):
    say(f"  γ_{n:<2d} = {mp.nstr(GAM[n],12)}")
say(f"  γ₂/γ₁ = {mp.nstr(GAM[2]/GAM[1],12)}   4cos²(π/88) = {mp.nstr(4*mp.cos(mp.pi/N)**2,12)}")
say(f"  差 = {mp.nstr(GAM[2]/GAM[1] - 4*mp.cos(mp.pi/N)**2, 3)}（恆等 ✓）")

# ---------------------------------------------------------------- 觀測窗
A_D = mp.mpf("3.717438495")      # ρ_u=0 分支之 α_s/α_u
S_A = mp.mpf("0.412601")
LO, HI = A_D - S_A, A_D + S_A
hdr("§1 觀測窗（ρ_u=0 分支，1σ）")
say(f"  α_s/α_u = {mp.nstr(A_D,10)} ± {mp.nstr(S_A,6)} ⟹ [{mp.nstr(LO,8)}, {mp.nstr(HI,8)}]")
say("  註：精確反解窗 [3.34606, 4.18155]；本節用對稱窗作掃描基準（兩者窗內整數皆 4）")

# ---------------------------------------------------------------- 掃描 γ_a/γ_b
hdr("§2 鏈譜比掃描：γ_a/γ_b（a,b = 1..44，共 1936 個有序對）")
hits = []
for a in range(1, 45):
    for b in range(1, 45):
        if a == b:
            continue
        r = GAM[a] / GAM[b]
        if LO <= r <= HI:
            hits.append((a, b, r))
say(f"  落窗：{len(hits)} / {44*43} = {mp.nstr(mp.mpf(len(hits))/1892*100,5)}%")
say("  最簡者（|a−b| 小）：")
hits.sort(key=lambda t: (abs(t[0] - t[1]), t[0] + t[1]))
for a, b, r in hits[:10]:
    say(f"    γ_{a}/γ_{b} = {mp.nstr(r,10)}   （a/b = {a}/{b}）")

hdr("§3 密度對照：同寬窗下『隨機』應有幾多")
# 鏈譜比之分佈：均勻對數分佈不適用；改以實際排序作經驗分佈
allr = sorted(GAM[a] / GAM[b] for a in range(1, 45) for b in range(1, 45) if a != b)
say(f"  全 1892 個比值之範圍：[{mp.nstr(allr[0],6)}, {mp.nstr(allr[-1],6)}]")
say(f"  中位數 = {mp.nstr(allr[len(allr)//2],8)}；於 3–5 之間者 = "
    f"{sum(1 for r in allr if 3 <= r <= 5)} 個")
say(f"  窗寬 = {mp.nstr(HI-LO,6)}（佔 3–5 區間之 {mp.nstr((HI-LO)/2*100,5)}%）")
say(f"  ★ 判：窗內 {len(hits)} 個命中（{mp.nstr(mp.mpf(len(hits))/1892*100,4)}%）"
    "係『翻倍結構之聚集』（見 §4），並非隨機；")
say(f"     {len(hits)} 個命中已足以說明『落窗』本身**無鑑別力** —— 落窗只值 "
    f"−log₂({mp.nstr(mp.mpf(len(hits))/1892,4)}) = {mp.nstr(-mp.log(mp.mpf(len(hits))/1892,2),4)} bit")

hdr("§4 為何聚集：γ_{2n}/γ_n = 4cos²(nπ/88) 之近 4 結構")
for n in (1, 2, 3, 5, 11, 22, 44 // 2):
    if 2 * n <= 44:
        say(f"  γ_{2*n}/γ_{n} = {mp.nstr(GAM[2*n]/GAM[n],10)}   （4cos²({n}π/88)）")
say("  ⟹ 任何『翻倍指標』之比都近 4（n 細時偏差 <1.6%），故落窗者**未必**有結構意義")

hdr("§5 候選分級")
cands = [
    ("B = 4（整數，bare）", mp.mpf(4), "無機制：純標籤", "docs/12"),
    ("γ₂/γ₁ = 4cos²(π/88)", GAM[2] / GAM[1], "機制：翻倍映射 n→2n（已註冊：γ_n(88)=γ̃_2n(176)）", "本節"),
    ("4/cos²(π/88)", 4 / mp.cos(mp.pi / N) ** 2, "無機制", "—"),
    ("2/sin²(π/88)·…", mp.mpf(0), "構造無意義（略）", "—"),
]
say(f"  {'候選':26s} {'值':>12s} {'z':>8s}  機制")
for nm, v, why, ref in cands:
    if v == 0:
        continue
    say(f"  {nm:26s} {mp.nstr(v,10):>12s} {mp.nstr((v-A_D)/S_A,6):>8s}  {why}")

hdr("§6 兩候選之可分辨性（B=4 vs γ₂/γ₁）")
r4 = mp.mpf(3) / (2 * mp.mpf(4))
rg = mp.mpf(3) / (2 * (GAM[2] / GAM[1]))
say(f"  ρ_d(B=4)      = {mp.nstr(r4,12)}")
say(f"  ρ_d(γ₂/γ₁)    = {mp.nstr(rg,12)}")
say(f"  差 = {mp.nstr(abs(r4-rg),8)}；現 σ_ρd = 0.04478517819 ⟹ 差 = "
    f"{mp.nstr(abs(r4-rg)/mp.mpf('0.04478517819'),6)} σ")
say(f"  要 3σ 分辨需 σ_ρd ≤ {mp.nstr(abs(r4-rg)/3,8)}（現值之 "
    f"{mp.nstr(mp.mpf('0.04478517819')/(abs(r4-rg)/3),5)}×）⟹ **現精度永不可分**")
say("  ⟹ 兩候選在觀測上永久等價；分流只能靠結構（機制有無），不能靠數據。")

# ---------------------------------------------------------------- R2 預測
hdr("§7 R2 — BESIII 精度情境下之鑑別力預測")
dX = mp.mpf("6.85")
dN = mp.mpf("1.29333234")
s_Xm = mp.mpf("0.07")
say(f"  {'σ(Ξ⁰)':>10s} {'σ(trace)':>12s} {'z(B=5)':>9s} {'z(15/4)':>9s} {'判':>12s}")
res = {}
for sx0 in ("0.20", "0.10", "0.05", "0.013", "0.0015"):
    s0 = mp.mpf(sx0)
    st = mp.sqrt((2 * mp.mpf("5.0e-7") / dX) ** 2 + (2 * mp.mpf("2.9e-7") / dX) ** 2
                 + (2 * dN * s_Xm / dX ** 2) ** 2 + (2 * dN * s0 / dX ** 2) ** 2)
    tr = 2 * dN / dX
    z5 = (tr - mp.mpf(4) / 10) / st
    z154 = (tr - mp.mpf(3) / 8) / st
    verdict = "可分 5" if abs((mp.mpf(4)/10 - mp.mpf("0.4"))/st) > 5 and abs(z5) > 5 else (
        "可分 3" if abs(z5) > 3 else "不可分")
    res[sx0] = (float(st), float(z5))
    say(f"  {sx0:>10s} {mp.nstr(st,8):>12s} {mp.nstr(z5,6):>9s} {'—':>9s} {verdict:>12s}")
zmax = abs((2 * dN / dX - mp.mpf(4) / 10)
           / mp.sqrt((2 * mp.mpf("5.0e-7") / dX) ** 2 + (2 * mp.mpf("2.9e-7") / dX) ** 2
                     + (2 * dN * s_Xm / dX ** 2) ** 2))
say(f"  ★ 判：σ(Ξ⁰) ≤ 0.013（BESIII 系統位階）已給 |z(B=5)| = 5.70 > 5σ；")
say(f"     天花板（σ(Ξ⁰)→0，受 Ξ⁻ 0.07 限制）= {mp.nstr(zmax,5)}σ ⟹ "
    "瓶頸係**發表**而非量度能力，且 5.8σ 係此通道之極限")

hdr("§匯總")
say(f"  ① 鏈譜比落窗 {len(hits)} 個 ⟹ 「落窗」無鑑別力（密度太高）")
say("  ② 但 γ₂/γ₁ = 4cos²(π/88) 係唯一『有已註冊機制』者（翻倍映射）")
say("  ③ B=4 與 γ₂/γ₁ 之差 0.00048（0.011σ），現精度**永久不可分**")
say("  ④ ⟹ R1 由『無來源』升為『兩個零參數候選，一個有機制』；"
    "分流靠結構不靠數據")
say(f"  ⑤ R2：σ(Ξ⁰) ≤ 0.013 已給 |z(B=5)| = 5.70；天花板 {mp.nstr(zmax,4)}σ（受 Ξ⁻ 限制）"
    "⟹ 瓶頸係發表非能力")

RESULT = {
    "window": [float(LO), float(HI)],
    "n_chain_hits": len(hits),
    "frac_hits": len(hits) / 1892,
    "gamma2_over_gamma1": str(GAM[2] / GAM[1]),
    "rho_d_B4": str(r4), "rho_d_gamma": str(rg),
    "diff_rho_d_sigma": float(abs(r4 - rg) / mp.mpf("0.04478517819")),
    "besiii_forecast": {k: v for k, v in res.items()},
}
OUT.write_text(json.dumps(RESULT, ensure_ascii=False, indent=2), encoding="utf-8")
say()
say("JSON → data/run_22_chain_ratio.json")
sys.exit(0)