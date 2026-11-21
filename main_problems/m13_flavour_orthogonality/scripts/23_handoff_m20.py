#!/usr/bin/python3
"""23 — Handoff 處理（M20 請求 ＋ M15 回執）之算術核

M20 請求：α_s/α_u = 4 可否給零參數 provenance？可 ⟹ α 作第二權重來源。
本腳本提供三件事：
  ① 權威數值重算（ρ_u、ρ_d、其 σ）＋ 兩個 z（用正典 σ vs M20 引用之 stale σ）
  ② α=(1,2,4) 之兩候選（B=4 / γ₂/γ₁）狀態與可分辨性
  ③ 獨立性審計：α 係由同一批強子質量反解 ⟹ 循環風險之量化陳述
零 import 其他腳本；只讀 vendored PDG；mpmath dps=40。
"""
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40
ROOT = Path(__file__).resolve().parent.parent
PDG = ROOT.parent / "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt"
OUT = ROOT / "data/run_23_handoff_m20.json"
out = []


def say(s=""):
    print(s)
    out.append(s)


def hdr(s):
    say()
    say("=" * 78)
    say(s)
    say("=" * 78)


def row(pid):
    for line in PDG.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("*") or len(line) < 70:
            continue
        cols = [line[0:8], line[8:16], line[16:24], line[24:32]]
        if any(c.strip() == pid for c in cols):
            return line
    raise SystemExit(f"PDG 缺 {pid}")


def m(pid):
    return mp.mpf(row(pid)[33:51].strip()) * 1000


def s(pid):
    v = row(pid)[51:60].strip()
    return (mp.mpf(v) if v else mp.mpf(0)) * 1000


n, p = m("2112"), m("2212")
Sm, S0, Sp = m("3112"), m("3212"), m("3222")
Xm, X0 = m("3312"), m("3322")
dN, dS, dX = n - p, Sm - S0, Xm - X0
sN = mp.sqrt(s("2112") ** 2 + s("2212") ** 2)
sS = mp.sqrt(s("3112") ** 2 + s("3212") ** 2)
sX = mp.sqrt(s("3312") ** 2 + s("3322") ** 2)
d_s, d_d = dX / 2, dS - dX / 2
d_u = dN - d_d
rho_u, rho_d = d_u / d_s, d_d / d_s
# σ 傳播（與 l2 同式）
sd_s = sX / 2
sd_d = mp.sqrt(sS ** 2 + sd_s ** 2)
sd_u = mp.sqrt(sN ** 2 + sd_d ** 2)
sr_u = mp.sqrt((sd_u / d_s) ** 2 + (d_u * sd_s / d_s ** 2) ** 2)
sr_d = mp.sqrt((sd_d / d_s) ** 2 + (d_d * sd_s / d_s ** 2) ** 2)

hdr("§0 權威重算（ρ 值精確；σ 兩種約定）")
say(f"  ρ_u = {mp.nstr(rho_u,12)}   ρ_d = {mp.nstr(rho_d,12)}")
say(f"  σ（本檔之**獨立相加**，會低估）：σ_ρu = {mp.nstr(sr_u,8)}、σ_ρd = {mp.nstr(sr_d,8)}")
L2J = ROOT / "data/l2_independent_verify.json"
if L2J.exists():
    D = json.loads(L2J.read_text(encoding="utf-8"))["_derived"]
    sr_u = mp.mpf(D["sigma_rho_u"]["str"]); sr_d = mp.mpf(D["sigma_rho_d"]["str"])
    say(f"  σ（**正典**，含協方差；scripts/01／09／l2）：σ_ρu = {mp.nstr(sr_u,8)}、"
        f"σ_ρd = {mp.nstr(sr_d,8)}")
else:
    say("  ⚠ L2 JSON 缺席，只能用獨立相加值")
say(f"  迹  = {mp.nstr(rho_u+rho_d,12)}   Δ_N/Δ_Ξ = {mp.nstr(dN/dX,12)}")

hdr("§1 α=(1,2,4) 之預言與 z（兩種 σ 約定）")
tgt_u, tgt_d = mp.mpf(0), mp.mpf(3) / 8
say(f"  預言：ρ_u = 0（A₂ Cartan 零方向，精確）、ρ_d = 3/8")
say(f"  z（正典 σ）：ρ_u = {mp.nstr((tgt_u-rho_u)/sr_u,6)}；ρ_d = {mp.nstr((tgt_d-rho_d)/sr_d,6)}")
stale = mp.mpf("0.03510")      # M12 docs/38 §E 舊值，R27 後已 SUPERSEDED
say(f"  z（若用**獨立相加** σ_ρd = {mp.nstr(stale,5)}）：ρ_d = {mp.nstr((tgt_d-rho_d)/stale,6)}")
say(f"  ★ M20 來件報 z ≈ +0.79 / −0.81 ⟹ 佢哋用嘅係獨立相加（忽略 ρ_u–ρ_d 相關 corr=−0.99635）；")
say(f"    正典 z = (+{mp.nstr((tgt_u-rho_u)/sr_u,4)}, {mp.nstr((tgt_d-rho_d)/sr_d,4)})；"
    f"σ_ρd 被低估 {mp.nstr(sr_d/stale,5)}× ⟹ 判詞仍 1σ 內，但**幅度要改**")

hdr("§2 α 之零參數 provenance（兩候選狀態）")
G = lambda k: 4 * mp.sin(mp.pi * k / 88) ** 2
g21 = G(2) / G(1)
say(f"  (a) α_d = 2α_u：**精確結構**（A₂ Cartan 零方向；ρ_u=0 ⟺ α_d=2α_u）")
say(f"  (b) α_s/α_u：候選一 B=4（整數，無機制）；候選二 γ₂/γ₁ = {mp.nstr(g21,12)}（機制=翻倍，需世代序）")
say(f"      觀測 r = {mp.nstr(3/(2*rho_d),10)} ± {mp.nstr(3*sr_d/(2*rho_d**2),6)}")
say(f"      ρ_d 差 = {mp.nstr(abs(mp.mpf(3)/8 - mp.mpf(3)/(2*g21)),8)}"
    f"（{mp.nstr(abs(mp.mpf(3)/8 - mp.mpf(3)/(2*g21))/sr_d,5)}σ）⟹ 永久不可分")

hdr("§3 獨立性審計（M20 之『第二權重來源』可行性）")
say("  α 係由 **同一批強子質量**（Δ_N、Δ_Σ、Δ_Ξ）反解：")
say("    δ_s = Δ_Ξ/2、δ_d = Δ_Σ − δ_s、δ_u = Δ_N − δ_d ⟹ ρ ⟹ α")
say("  ⟹ 任何以強子質量為輸入之檢定，用 α 作『獨立第二權重』= **循環**（同一數據兩次入帳）")
say("  ⟹ α 只可在**不同扇區**（CKM／PMNS／宇宙學）作固定比 (1,2,4) 使用；")
say("     且其值之『固定性』受 M14 Q2（世代序）限制 ⟹ 現階段屬**條件性固定**")

hdr("§匯總")
say(f"  ① 正典：ρ_u = {mp.nstr(rho_u,8)} ± {mp.nstr(sr_u,6)}；ρ_d = {mp.nstr(rho_d,10)} ± {mp.nstr(sr_d,6)}")
say(f"  ② α=(1,2,4) ⟹ z = (+{mp.nstr((tgt_u-rho_u)/sr_u,4)}, {mp.nstr((tgt_d-rho_d)/sr_d,4)})")
say(f"  ③ M20 之 z=−0.81 來自獨立相加（忽略 corr=−0.99635），σ_ρd 低估 1.276× ⟹ 需更正")
say(f"  ④ provenance：α_d=2α_u 可（結構）；α_s/α_u 條件性（兩候選永久等價，一有機制）")
say(f"  ⑤ 作獨立權重：只在**不同扇區**可行；同批強子質量則循環")

RESULT = {
    "rho_u": str(rho_u), "sigma_rho_u": str(sr_u),
    "rho_d": str(rho_d), "sigma_rho_d": str(sr_d),
    "z_u_canonical": float((tgt_u - rho_u) / sr_u),
    "z_d_canonical": float((tgt_d - rho_d) / sr_d),
    "z_d_with_stale_sigma": float((tgt_d - rho_d) / stale),
    "alpha_ratio_candidate_gamma": str(g21),
    "degeneracy_sigma": float(abs(mp.mpf(3)/8 - mp.mpf(3)/(2*g21)) / sr_d),
}
OUT.write_text(json.dumps(RESULT, ensure_ascii=False, indent=2), encoding="utf-8")
say()
say("JSON → data/run_23_handoff_m20.json")
sys.exit(0)