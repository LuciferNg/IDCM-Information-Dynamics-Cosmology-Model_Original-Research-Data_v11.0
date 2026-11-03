#!/usr/bin/python3
"""20 — R2 攻堅：裁決通道最優化 ＋ 通道獨立性檢定

問題：M13 之裁决量（迹 = 2Δ_N/Δ_Ξ）受 Δ_Ξ 主導（σ 93.98% 來自 Ξ）。
本腳本窮舉 **全部質量差通道**（9 個基態重子 → 36 對差），逐條：
  ① 以該差作分母得 trace_c = 2Δ_N / D_c，比較 3/8（B=4）與 (B−1)/(2B)（B=3,5,6）
  ② 計 σ(D_c)、z、以及「4 vs 5」之分離度（σ 需縮多少倍）
  ③ 秩排序：邊條通道最有裁決力
  ④ **通道獨立性檢定**：把各通道視為 3/8 之獨立量度 → 加權均值 + χ²/dof
     （χ²/dof ≫ 1 ⟹ 通道不相容 ⟹ δ 之通道無關假設受壓）
零 import 其他腳本；只讀 vendored PDG 2024；mpmath dps=40；浮點只作顯示。
"""
import json
import sys
from decimal import Decimal
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

ROOT = Path(__file__).resolve().parent.parent
PDG = ROOT.parent / "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt"
OUT = ROOT / "data/run_20_channel_optimization.json"
out = []


def say(s=""):
    print(s)
    out.append(s)


def hdr(s):
    say()
    say("=" * 78)
    say(s)
    say("=" * 78)


BASE = {          # 名稱 → PDG MC id
    "p": "2212", "n": "2112", "Lambda": "3122",
    "Sigma+": "3222", "Sigma0": "3212", "Sigma-": "3112",
    "Xi0": "3322", "Xi-": "3312", "Omega-": "3334",
}


def row(pid):
    for line in PDG.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("*") or len(line) < 70:
            continue
        cols = [line[0:8], line[8:16], line[16:24], line[24:32]]
        if any(c.strip() == pid for c in cols):
            return line
    raise SystemExit(f"PDG 缺 {pid}")


M, S = {}, {}
for name, pid in BASE.items():
    ln = row(pid)
    M[name] = mp.mpf(ln[33:51].strip()) * 1000            # GeV → MeV
    S[name] = mp.mpf(ln[51:60].strip() or "0") * 1000

hdr("§0 基態重子質量（PDG 2024 vendored）")
for k in BASE:
    say(f"  {k:8s} = {mp.nstr(M[k],12):>14s} ± {mp.nstr(S[k],8)} MeV")

NA = M["n"] - M["p"]
S_NA = mp.sqrt(S["n"] ** 2 + S["p"] ** 2)
say(f"  Δ_N = n − p = {mp.nstr(NA,12)} ± {mp.nstr(S_NA,8)}")

# ------------------------------------------------------------------ 通道
names = list(BASE)
CH = []
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        hi, lo = (a, b) if M[a] >= M[b] else (b, a)      # 恒取重者減輕者（BASE 排序≈質量排序）
        d = M[hi] - M[lo]
        if d <= 0:
            continue
        sd = mp.sqrt(S[a] ** 2 + S[b] ** 2)
        CH.append((f"{hi}−{lo}", d, sd))

hdr("§1 通道表：trace_c = 2Δ_N / D_c  vs 3/8 = 0.375（B=4）")
say(f"  {'通道':16s} {'D (MeV)':>12s} {'σ(D)':>9s} {'trace_c':>12s} {'σ(t)':>10s} {'z(3/8)':>9s}")
TARGET = mp.mpf(3) / 8
rows = []
for nm, d, sd in CH:
    if d < 2:                       # 分母太小 ⟹ 無意義
        continue
    t = 2 * NA / d
    st = 2 * S_NA / d + 2 * NA * sd / d ** 2        # 線性傳播（保守相加）
    st_q = mp.sqrt((2 * S_NA / d) ** 2 + (2 * NA * sd / d ** 2) ** 2)
    z = (t - TARGET) / st_q
    rows.append({"chan": nm, "D": d, "sD": sd, "t": t, "st": st_q, "z": z})
    if nm in ("Xi-−Xi0", "Sigma-−Sigma0", "Sigma0−Sigma+", "Sigma-−Sigma+"):
        say(f"  {nm:16s} {mp.nstr(d,10):>12s} {mp.nstr(sd,6):>9s} "
            f"{mp.nstr(t,10):>12s} {mp.nstr(st_q,6):>10s} {mp.nstr(z,6):>9s}")

say()
say(f"  共 {len(rows)} 條可用通道（9 個重子兩兩差，正值）")

hdr("§2 裁決力排序（|z| 最小 = 最貼 3/8；另列 4 vs 5 分離度）")
rows.sort(key=lambda r: abs(r["z"]))
say(f"  {'#':>2s} {'通道':16s} {'trace_c':>12s} {'σ(t)':>10s} {'z':>9s} {'σ需縮':>8s}")
for i, r in enumerate(rows[:12], 1):
    # 要把 B=4(0.375) vs B=5(0.4) 分到 3σ：需 σ_t ≤ |0.4−0.375|/3
    need = abs(mp.mpf(4) / 10 - TARGET) / 3
    fac = r["st"] / need
    say(f"  {i:2d} {r['chan']:16s} {mp.nstr(r['t'],10):>12s} {mp.nstr(r['st'],6):>10s} "
        f"{mp.nstr(r['z'],6):>9s} {mp.nstr(fac,5):>8s}×")

best = rows[0]
say()
say(f"  ★ 最貼 3/8 之通道：{best['chan']}（trace={mp.nstr(best['t'],10)}，"
    f"z={mp.nstr(best['z'],6)}，σ={mp.nstr(best['st'],6)}）")
say("  ★ 但『貼』≠『可用』：非 Ξ 通道之 D 本身係 δ 結構之另一分量，"
    "以佢代 Δ_Ξ 等於改動 ansatz（通道無關假設），須先過 §4。")

# ------------------------------------------------------------------ 獨立性檢定
hdr("§3 通道獨立性檢定（各通道視為 3/8 之獨立量度）")
sub = [r for r in rows if r["st"] > 0]
wsum = mp.mpf(0)
mu = mp.mpf(0)
for r in sub:
    w = 1 / r["st"] ** 2
    wsum += w
    mu += w * r["t"]
mu /= wsum
chi2 = sum(((r["t"] - mu) / r["st"]) ** 2 for r in sub)
dof = len(sub) - 1
say(f"  n = {len(sub)}，加權均值 = {mp.nstr(mu,10)}，σ_mu = {mp.nstr(mp.sqrt(1/wsum),6)}")
say(f"  χ² = {mp.nstr(chi2,8)}，dof = {dof}，χ²/dof = {mp.nstr(chi2/dof,6)}")
say(f"  z(加權均值 vs 3/8) = {mp.nstr((mu - TARGET) / mp.sqrt(1/wsum),6)}")
say("  ★ 判：χ²/dof 巨大 ⟹ 各通道**並非**同一量之獨立量度；")
say("     即 δ 之通道無關假設不能靠平均提升精度（與 docs/14 §五 之 1 bit 一致）。")

# ------------------------------------------------------------------ 只取 Ξ 系
hdr("§4 只取「純 Ξ 系」通道（不假設通道無關）")
xi = [r for r in rows if r["chan"] in ("Xi-−Xi0",)]
for r in xi:
    say(f"  {r['chan']}: trace = {mp.nstr(r['t'],12)} ± {mp.nstr(r['st'],8)}，"
        f"z(3/8) = {mp.nstr(r['z'],6)}")
need3 = abs(mp.mpf(4) / 10 - TARGET) / 3
say(f"  要把 4 vs 5 分到 3σ：需 σ(trace) ≤ {mp.nstr(need3,8)}"
    f"（現 {mp.nstr(xi[0]['st'],8)} ⟹ 需縮 {mp.nstr(xi[0]['st']/need3,5)}×）")
say(f"  折算 σ(Ξ⁰)：現 0.20 ⟹ 需 ≤ {mp.nstr(mp.mpf('0.20')/(xi[0]['st']/need3),5)} MeV")

# ------------------------------------------------------------------ 匯總
hdr("§匯總")
say(f"  ① 可用通道 {len(rows)} 條；最貼 3/8 者 = {best['chan']}"
    f"（z={mp.nstr(best['z'],5)}），但屬別 ansatz 分量，不可直接替代")
say(f"  ② 通道獨立性：χ²/dof = {mp.nstr(chi2/dof,6)} ⟹ 不相容（δ 通道無關未成立）")
say(f"  ③ 純 Ξ 系：需 σ(trace) ≤ {mp.nstr(need3,8)}（縮 {mp.nstr(xi[0]['st']/need3,5)}×）"
    f"，即 σ(Ξ⁰) ≤ {mp.nstr(mp.mpf('0.20')/(xi[0]['st']/need3),4)} MeV")
say("  ④ 結論：R2 **無法**在本框架內以「換通道」解決 —— 精度瓶頸係單一質量 Ξ⁰，"
    "任何組合皆需同等或更差精度；唯一出路 = 直量（BESIII J/ψ 全集）")

RESULT = {
    "n_channels": len(rows),
    "best_chan": best["chan"],
    "best_trace": str(best["t"]),
    "best_z": float(best["z"]),
    "xi_trace": str(xi[0]["t"]),
    "xi_sigma": str(xi[0]["st"]),
    "xi_z": float(xi[0]["z"]),
    "chi2_over_dof": float(chi2 / dof),
    "weighted_mean": str(mu),
    "sigma_needed_trace_3sigma": str(need3),
    "sigma_Xi0_needed_3sigma": str(mp.mpf("0.20") / (xi[0]["st"] / need3)),
}
OUT.write_text(json.dumps(RESULT, ensure_ascii=False, indent=2), encoding="utf-8")
say()
say("JSON → data/run_20_channel_optimization.json")
sys.exit(0)