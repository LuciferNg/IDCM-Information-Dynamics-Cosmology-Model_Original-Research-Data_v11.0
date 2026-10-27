"""M7 scripts/12 — O209「唯一分配定理」之獨立重建 ＋ 與 O120 載體衝突 ＋ 對 M7 之影響
================================================================================================================
背景：O209 register `conclusion`（2026-10-24）稱「分配由冪律 ansatz 唯一確定 ⟹ 全局方案數 0→1」，
      即 u/d 指派=u@11, d@9（與 O120 之 γ 升序 u@2, d@3 **極性相反**）。
依據 skill `cross-record-result-import-verification`：§1 由文檔重建、§2 約定審計、§3 判詞三級。

輸入宣告：
  · O136 `data/63_odd_v_factor_theorem.json` 之 `dm2_M9`（|Δm^(2)| 全精度）——**取為輸入**，
    本檔不自建 88×88 鏈算子（重建成本高，且該值已屬 O136 已註冊產物；標 ⚠️ 未獨立重算）
  · O209 `scripts/03` 之 MASS 表（9 個帶電費米子，MeV）
  · O120/O234 凍結層（m_d/m_u = 2.150538 ± 0.110995）
輸出：data/12_o209_assignment.json
"""
import itertools
import json
import os
from mpmath import mp, mpf, sin, pi, log, sqrt, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
REPO = os.path.abspath(os.path.join(M7, "..", "..", ".."))
OUT = {"module": "M7", "script": "12_o209_assignment.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<52} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 106 + f"\n{s}\n" + "=" * 106)


# ── 輸入 ──────────────────────────────────────────────────────────────
DM = {int(k): abs(mpf(str(v))) for k, v in json.load(
    open(os.path.join(REPO, "IDCM-v10.0--It-From-Beat/open_problems_fresh/"
                            "o136_quark_mass_spectrum/data/63_odd_v_factor_theorem.json"))
)["dm2_M9"].items()}
MASS = {"e": 0.51099895, "mu": 105.6583755, "tau": 1776.86,
        "u": 2.16, "d": 4.67, "s": 93.4, "c": 1270.0, "b": 4180.0, "t": 172690.0}
MODES = sorted(DM)                                   # [1,2,3,4,5,7,9,11,13]
tgt_ud, sig_ud_rel = mpf("2.150538"), mpf("0.110995") / mpf("2.150538")
M0 = mpf("27.900915")
ND = 88
G = {n: 4 * sin(n * pi / ND) ** 2 for n in range(1, 45)}

print("=" * 106)
print("M7 scripts/12 — O209 唯一分配定理：重建 ＋ 載體衝突 ＋ 對 M7 之影響")
print("=" * 106)

# ══════════ Part A：重建 O209 §C1／C2 ══════════
hdr("Part A §1 — 模式集唯一性（C(9,9)=1）")
G2 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
import math
C44_9 = math.comb(44, 9)
C39_9 = math.comb(39, 9)
print(f"    C(44,9) = {C44_9}；C(39,9) = {C39_9}；|G2| = {len(G2)} ⟹ C(|G2|,9) = {math.comb(len(G2), 9)}")
ok("A1 C(44,9) = 708930508（對 O209 data）", C44_9 == 708930508, C44_9, 708930508)
ok("A1 C(39,9) = 211915132", C39_9 == 211915132, C39_9, 211915132)
ok("A1 G2 通過集 = 9 個 ⟹ 9-子集數 = 1（組合論證，無需暴力）",
   len(G2) == 9 and math.comb(len(G2), 9) == 1, f"C({len(G2)},9) = 1", 1)

hdr("Part A §2 — 分配唯一性（窮舉 9! = 362,880）")
parts = sorted(MASS, key=lambda p: -MASS[p])
logs_m = [log(mpf(MASS[p])) for p in parts]
logs_d = [log(DM[n]) for n in MODES]
nm, nd_ = len(parts), len(MODES)
best, min_rms, n_tie = None, None, 0
rows = []
for perm in itertools.permutations(range(nm)):        # perm[i] = 派給 MODES[i] 之粒子索引
    xs = [logs_d[i] for i in range(nm)]
    ys = [logs_m[perm[i]] for i in range(nm)]
    mx, my = sum(xs) / nm, sum(ys) / nm
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    corr = sxy / sqrt(sxx * syy)
    rms = sqrt(syy / nm) * sqrt(1 - corr ** 2)        # 最佳擬合下之殘餘 RMS（ln 空間）
    if min_rms is None or rms < min_rms - mpf("1e-9"):
        min_rms, best, n_tie = rms, perm, 1
    elif abs(rms - min_rms) <= mpf("1e-9"):
        n_tie += 1

pair = {MODES[i]: parts[best[i]] for i in range(nm)}
print(f"    最小 RMS = {nstr(min_rms, 10)}；達最小（±1e−9）之配對數 = {n_tie}")
print("    唯一配對：" + "、".join(f"{k}→{pair[k]}" for k in MODES))
ok("A2 唯一最小配對數 = 1（窮舉 362,880）", n_tie == 1, n_tie, 1)
ok("A2 最小 RMS = 0.587228（對 O209 0.5872281）", abs(min_rms - mpf("0.587228")) < mpf("1e-6"),
   nstr(min_rms, 9), "0.587228")
ok("A2 配對 = 單調（|Δm| 遞減 ↔ 質量遞增）",
   [pair[n] for n in MODES] == ["t", "b", "tau", "c", "mu", "s", "d", "u", "e"],
   [pair[n] for n in MODES], "['t','b','tau','c','mu','s','d','u','e']")
ok("A2 ★ O209 之 u/d 指派 = u@11、d@9（與 O120 之 2/3 極性相反）",
   pair[11] == "u" and pair[9] == "d", f"u@{[k for k, v in pair.items() if v == 'u'][0]}, "
   f"d@{[k for k, v in pair.items() if v == 'd'][0]}", "u@11, d@9")

hdr("Part A §3 — 冪律指數 α 與其預言")
mx, my = sum(logs_d) / nm, sum(logs_m) / nm
sxy = sum((a - mx) * (b - my) for a, b in zip(logs_d, logs_m))
sxx = sum((a - mx) ** 2 for a in logs_d)
alpha = sxy / sxx
syy = sum((b - my) ** 2 for b in logs_m)
corr = sxy / sqrt(sxx * syy)
print(f"    α（ln m 對 ln|Δm| 斜率）= {nstr(alpha, 8)}；corr = {nstr(corr, 8)}")
ok("A3 corr = 0.988353（對 O209 記載）", abs(corr - mpf("0.988353")) < mpf("1e-5"),
   nstr(corr, 9), "0.988353")
ratio_pl = (DM[9] / DM[11]) ** alpha
z_pl = log(ratio_pl / tgt_ud) / sig_ud_rel
print(f"    冪律載體下 u/d：m_d/m_u = (|Δm_9|/|Δm_11|)^α = {nstr(ratio_pl, 8)}"
      f"  vs {nstr(tgt_ud, 8)} ⟹ z = {float(z_pl):+.2f}")
ok("A3 ★判詞：O209 唯一指派在自身載體下 u/d 預言被否證（|z|>5）", abs(z_pl) > 5,
   f"z = {float(z_pl):+.2f}", "|z|>5",
   "記負：唯一性（排序）穩健，但定量比值 3.26 vs 2.15 ⟹ 指派內容不可用；"
   "注意唯一性係 α-無關（rearrangement inequality），故此否證與 α 取值無關")

# ══════════ Part B：載體衝突 ══════════
hdr("Part B — 載體衝突：γ_n（O120）vs |Δm|^α（O209）")
r_gamma_o209 = G[9] / G[11]
z_go = log(r_gamma_o209 / tgt_ud) / sig_ud_rel
r_gamma_o120 = G[3] / G[2]
z_g1 = log(r_gamma_o120 / tgt_ud) / sig_ud_rel
print(f"    γ 載體 + O209 指派（d@9, u@11）：γ₉/γ₁₁ = {nstr(r_gamma_o209, 8)}"
      f"  ⟹ z = {float(z_go):+.2f}（**m_u > m_d，符號反**）")
print(f"    γ 載體 + O120 指派（u@2, d@3）：γ₃/γ₂ = {nstr(r_gamma_o120, 8)}"
      f"  ⟹ z = {float(z_g1):+.2f}")
r_dl_o120 = (DM[2] / DM[3]) ** alpha
z_dl = log(r_dl_o120 / tgt_ud) / sig_ud_rel
print(f"    |Δm|^α 載體 + O120 指派（u@2, d@3）：(|Δm_2|/|Δm_3|)^α = {nstr(r_dl_o120, 8)}"
      f"  ⟹ z = {float(z_dl):+.2f}")
ok("B1 γ 載體下 O209 指派被否證（|z| > 5）", abs(z_go) > 5, f"z = {float(z_go):+.2f}", "|z|>5",
   "⟹ O209 之 u/d 指派**唔可以**搬去 γ 載體")
ok("B2 |Δm|^α 載體下 O120 指派被否證（|z| > 5）", abs(z_dl) > 5, f"z = {float(z_dl):+.2f}", "|z|>5",
   "⟹ 兩完成品互斥：唯一性各自只在**自己載體**內成立")

# ══════════ Part C：對 M7 之影響 ══════════
hdr("Part C — 對 M7 之影響")
AN = json.load(open(os.path.join(M7, "data", "01_reproduction_gate.json")))["anchors"]
AK = {"N": "N", "Σ": "Sigma", "Ξ": "Xi"}
HADR = [("N", {"u": 2, "d": 1}, 1), ("Σ", {"u": 2, "s": 1}, 2), ("Ξ", {"u": 1, "s": 2}, 1)]


def chi2_assign(amap):
    tot = mpf(0)
    for nm, nq, coeff in HADR:
        w = {}
        for f, c in nq.items():
            w[amap[f]] = w.get(amap[f], 0) + c
        val = mpf(w.get(amap["u"], 0) - w.get(amap["d"], 0)) / coeff
        tot += ((val - mpf(AN[AK[nm]][0])) / mpf(AN[AK[nm]][1])) ** 2
    return tot


a_o120 = {"u": 2, "d": 3, "s": 4}
a_o209 = {"u": 11, "d": 9, "s": 7}
c120, c209 = chi2_assign(a_o120), chi2_assign(a_o209)
print(f"    ⟨O⟩≡1 之 χ²：O120 指派 {nstr(c120, 8)}；O209 指派 {nstr(c209, 8)}（應為同一值）")
ok("C1 兩個『已完成』指派皆為單射 ⟹ M7 之 ⟨O⟩≡1 與 χ² 不變",
   abs(float(c120) - 18.220726674) < 1e-6 and abs(float(c209) - float(c120)) < 1e-9,
   f"{nstr(c120,8)} / {nstr(c209,8)}", "18.220727（同一）",
   "⟹ M7 主判對『邊個完成品』零依賴")

X120 = M0 * mpf(18) / 5 * (G[3] - G[2])
X209 = M0 * mpf(18) / 5 * (G[9] - G[11])
print(f"    X = m_d − m_u（γ 載體、C₁=18/5）：O120 指派 {nstr(X120,8)} MeV"
      f"；O209 指派 {nstr(X209,8)} MeV（**負號**）")
print(f"    M7／FLAG 錨：X = 2.56 ± 0.094 MeV")
ok(f"C2 X 分支依賴：O120 +2.546（z=−0.15）、O209 {nstr(X209,6)}（符號反）",
   X120 > 0 and X209 < 0, f"{nstr(X120,6)} / {nstr(X209,6)}", "正 / 負",
   "⟹ A3 之 X 預測**不可升格**；須先裁『哪個載體／指派為正典』")

print("\n" + "=" * 106)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(min_rms=float(min_rms), n_tie=n_tie, corr=float(corr), alpha=float(alpha),
                        pairing={str(k): v for k, v in pair.items()},
                        ratio_pl_ud=float(ratio_pl), z_pl_ud=float(z_pl),
                        z_gamma_o209=float(z_go), z_gamma_o120=float(z_g1), z_dl_o120=float(z_dl),
                        chi2_o120=float(c120), chi2_o209=float(c209),
                        X_o120=float(X120), X_o209=float(X209)))
json.dump(OUT, open(os.path.join(M7, "data", "12_o209_assignment.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/12_o209_assignment.json")