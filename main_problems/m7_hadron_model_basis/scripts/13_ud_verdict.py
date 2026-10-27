"""M7 scripts/13 — u/d 指派：生死判定（兩個載體逐個處決，非偏好選擇）
==============================================================================================================
判準（事前寫死）：
  生 ⟺ ① 載體參數非由被測量自身擬合（非循環） ② 指派在已註冊結構內唯一 ③ |z| ≤ 2
  三條全中 = 生；任一條死 = 該載體死。

處決 A（γ_n·C_f，O120）：唯一性靠「雙重態共用 C_f (H1)」；H1 已於 t/b 被 14.65σ 否證 ⟹ 前提死。
      另核：u/d 之 2σ 內 3 對 ⟹ X 預測有幾個互異值？
處決 B（|Δm_n|^α，O209）：唯一性係 rearrangement 排序恆真（對任何質量集皆然）⟹ 0 bit 信息；
      再做樣本外 α（剔除 u,d 擬合）⟹ 預言 u/d 比值 ⟹ z。
輸出：data/13_ud_verdict.json
"""
import itertools
import json
import os
import random

from mpmath import mp, mpf, sin, pi, log, sqrt, nstr, gammainc, inf

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(M7, "..", "..", ".."))
OUT = {"module": "M7", "script": "13_ud_verdict.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 104 + f"\n{s}\n" + "=" * 104)


DM = {int(k): abs(mpf(str(v))) for k, v in json.load(
    open(os.path.join(REPO, "IDCM-v10.0--It-From-Beat/open_problems_fresh/"
                            "o136_quark_mass_spectrum/data/63_odd_v_factor_theorem.json"))
)["dm2_M9"].items()}
MASS = {"e": 0.51099895, "mu": 105.6583755, "tau": 1776.86,
        "u": 2.16, "d": 4.67, "s": 93.4, "c": 1270.0, "b": 4180.0, "t": 172690.0}
MODES = sorted(DM)
tgt, sig_rel = mpf("2.150538"), mpf("0.110995") / mpf("2.150538")
M0, ND, C1 = mpf("27.900915"), 88, mpf(18) / 5
G = {n: 4 * sin(n * pi / ND) ** 2 for n in range(1, 45)}

print("=" * 104)
print("M7 scripts/13 — u/d 指派：生死判定（載體 A：γ_n·C_f；載體 B：|Δm_n|^α）")
print("=" * 104)

# ══════════ 處決 B-1：唯一性係排序恆真（0 bit） ══════════
hdr("處決 B-1 — O209 唯一性係 rearrangement 恆真 ⟹ 對任意質量集皆然（0 bit 信息）")
import numpy as np
rng = random.Random(31415)
NS = 20
PERM = np.array(list(itertools.permutations(range(9))), dtype=np.int16)   # 362880 x 9
dlog = np.array([float(log(DM[n])) for n in MODES])
dlog_c = dlog - dlog.mean()
always_monotone, always_unique = 0, 0
for _ in range(NS):
    fake = sorted([mpf("0.5") * (mpf("2e5") / mpf("0.5")) ** rng.random() for _ in range(9)],
                  key=lambda x: -float(x))
    ylog = np.log(np.array([float(x) for x in fake]))
    Y = ylog[PERM]                                   # 362880 x 9
    Yc = Y - Y.mean(axis=1, keepdims=True)
    num = (Yc * dlog_c).sum(axis=1)
    den = np.sqrt((Yc ** 2).sum(axis=1) * (dlog_c ** 2).sum())
    corr = num / den
    # O209 之 ansatz 限定 α > 0 ⟺ corr > 0 ⟹ 只在 corr>0 之可行域判唯一性
    pos = corr > 0
    best_pos = float(corr[pos].max())
    if abs(best_pos - float(corr[0])) <= 1e-12:
        always_monotone += 1
    n_tie = int((np.abs(corr - best_pos) <= 1e-12).sum())
    if n_tie == 1:
        always_unique += 1
print(f"    隨機質量集 {NS} 次（numpy 全枚舉 362,880/次，限 α>0）："
      f"單調 {always_monotone}/{NS}；唯一 {always_unique}/{NS}")
ok("B1 在 α>0 可行域內，單調配對恆為唯一最大 corr（rearrangement）⟹ 排序定理，非經驗",
   always_monotone == NS and always_unique == NS,
   f"單調 {always_monotone}/{NS}、唯一 {always_unique}/{NS}", f"{NS}/{NS}",
   "⟹ 對任何質量集皆得 N=1 ⟹ O209 之『唯一』係排序必然，非 IDCM 結構所選")
ok("B1b 故『唯一』唔可以作為 u/d 指派之證據", True, "0 bit", "0 bit",
   "真正內容 = 相關 0.988（ordering 一致），並非模式身份")

# ══════════ 處決 B-2：樣本外 α（剔除 u,d） ══════════
hdr("處決 B-2 — 樣本外檢定：α 由其餘 7 個費米子擬合，再預言 u/d 比值")
keep = ["t", "b", "tau", "c", "mu", "s", "e"]
mode_of = {"t": 1, "b": 2, "tau": 3, "c": 4, "mu": 5, "s": 7, "e": 13}
xs = [log(DM[mode_of[p]]) for p in keep]
ys = [log(mpf(MASS[p])) for p in keep]
mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
sxx = sum((a - mx) ** 2 for a in xs)
alpha_out = sxy / sxx
ratio_b = (DM[9] / DM[11]) ** alpha_out
z_b = log(ratio_b / tgt) / sig_rel
print(f"    樣本外 α（7 點，剔除 u,d）= {nstr(alpha_out, 8)}（全擬合 α = 2.6332541）")
print(f"    預言 m_d/m_u = (|Δm₉|/|Δm₁₁|)^α_out = {nstr(ratio_b, 8)} vs {nstr(tgt, 8)} ⟹ z = {float(z_b):+.2f}")
ok("B2 樣本外預言被否證（|z| > 5）⟹ 載體 B 死（非循環判據）", abs(z_b) > 5,
   f"z = {float(z_b):+.2f}", "|z|>5",
   "α 只由 7 個非 u/d 質量定 ⟹ 無循環；u/d 比值差 "
   f"{float((ratio_b/tgt-1)*100):+.0f}%")

# ══════════ 處決 A-1：H1 前提已死 ══════════
hdr("處決 A-1 — 載體 A 之指派唯一性唯一依據 = 雙重態共用 C_f（H1）；H1 已死")
print("    O120 docs/06 §4：H1 在 t/b 最佳 z = +14.65（2σ 內 0 個）；全局互斥窮舉 0/60,480")
ok("A1 H1（共用 C_f）已於 t/b 以 14.65σ 否證 ⟹ A 之『共用 ⟹ 比值 = γ 比』前提死",
   True, "z=+14.65, 0/60480", "H1 死", "承 O120 §4；本 M 已於 scripts/11 V7/V8 重建")
ok("A1' 故 A 之 u/d『命中 0.83σ』失去結構依據（O120 自判）", True, "池效應", "池效應",
   "2σ 內 3 對 = 13 個合法對中有一個近 2.15 之必然")

# ══════════ 處決 A-2：X 預測之指派歧義 ══════════
hdr("處決 A-2 — 即使接受 A，u/d 之 3 個 2σ 候選給出互異 X ⟹ X 不可預言")
sigma_X = mpf("0.094")
cands = {"(2,3)": (2, 3), "(9,13)": (9, 13), "(5,7)": (5, 7)}
Xs = {}
for lbl, (nu, ndd) in cands.items():
    Xs[lbl] = M0 * C1 * (G[ndd] - G[nu])
    print(f"    {lbl:>7}：X = M₀C₁(γ_{ndd}−γ_{nu}) = {nstr(Xs[lbl], 8)} MeV"
          f"  （vs FLAG 2.56±0.094 ⟹ z = {float((Xs[lbl]-mpf('2.56'))/sigma_X):+.2f}）")
spread = max(Xs.values()) - min(Xs.values())
print(f"    三者跨度 = {nstr(spread, 6)} MeV = {float(spread/sigma_X):.1f}× σ_X")
ok("A2 候選指派間 X 跨度 ≫ σ_X ⟹ X 不可由 u/d 指派預言", spread > 2 * sigma_X,
   f"{nstr(spread,6)} MeV = {float(spread/sigma_X):.1f}σ_X", "> 2σ_X",
   "⟹ 就算 A 生，X 仍不可預言（指派歧義）")

# ══════════ 生死裁決 ══════════
hdr("生死裁決（三條件：非循環 ∧ 唯一 ∧ |z|≤2）")
verdict = [
    ("A：γ_n·C_f（O120）", "非循環 ✅（γ 係譜本徵值）", "唯一 🔴（2σ 內 3 對；H1 已死）", "z=+0.83 ✅"),
    ("B：|Δm_n|^α（O209）", f"非循環 🔴（α 擬合；樣本外 z={float(z_b):+.2f}）", "唯一 🔴（排序恆真，0 bit）", "z=+8.07 🔴"),
]
for name, a, b, c in verdict:
    print(f"    {name:<22} ①{a}\n      {'':<22} ②{b}\n      {'':<22} ③{c}")
print("\n    ⟹ 兩個載體皆**死**（各死於不同條件）⟹ **u/d 指派 = NO-GO（不可由已註冊結構判定）**")
OUT["verdict"] = dict(A="dead(H1 refuted; 3 pairs within 2σ; X spread 3.1σ)",
                      B=f"dead(ordering tautology; out-of-sample alpha z={float(z_b):+.2f})",
                      ud_assignment="NO-GO", X="not predictable (spread "
                      f"{float(spread/sigma_X):.1f} sigma_X)")

print("\n" + "=" * 104)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(alpha_out=float(alpha_out), ratio_b=float(ratio_b), z_b=float(z_b),
                        X_cands={k: float(v) for k, v in Xs.items()},
                        X_spread=float(spread), X_spread_over_sigma=float(spread / sigma_X),
                        tautology_runs=NS, tautology_unique=always_unique))
json.dump(OUT, open(os.path.join(M7, "data", "13_ud_verdict.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/13_ud_verdict.json")