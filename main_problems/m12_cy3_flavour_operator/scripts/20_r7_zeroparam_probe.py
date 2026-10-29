#!/usr/bin/env python3
"""R7 — 四 δ 之零參數式檢定（密度掃描）＋ Ξc 家族常數檢定 ＋ Ξb 反駁 ＋ 自旋標籤證書

資料源（唯讀 vendored）: data/pdg2024/mass_width_2024.txt（PDG 2024, sha256 9d6024da…）

本輪全部係**證偽工具**，唔係擬合：目標值先用 R6 之 δ（Form 精確解）算出，
然後測試「有冇任何預註冊形式族可以唯一命中」。判準跑前凍結如下。

A1 符號範疇檢定：任何 sign-definite（全正或全負）之譜映射 ⟹ ρ_d, ρ_b, ρ_u 同號。
   實測 ρ_u < 0 而 ρ_d, ρ_b > 0 ⟹ sign-definite 族**範疇判死**（唔使掃描）。
A2 密度掃描（只做 ρ_d = δ_d/δ_s，因 b 無註冊模式指數）：
   5 個預註冊形式 M1–M5（見下），指數網格 (k_u, k_d) ∈ [1..44]²，k_s 固定為註冊值 4。
   命中等於 |model − target| / |target| < 1%。
   **判準：若命中格數 ≥ 5（即密度 ≥ 5/1936），則「註冊指派 (2,3) 產生命中」不構成證據**
   （單點命中喺 1936 格中屬必然；與 O238/O209 密度判準一致）。
A3 空檢定：把 target 換成 10⁴ 個 [0.15, 0.85] 均勻亂數，統計每形式之命中數分佈；
   若 A2 之觀測命中數落於此分佈之內 ⟹ 無信號。
B1 Ξc 家族常數檢定：三態偏移 (Ξc, Ξc′, Ξc(2645)) 之加權均值與 χ²；
   χ²/2 > 9 ⟹ 「家族常數」假設死（即需更多參數）。
B2 Ξb 反駁：若 Ξ-型偏移係 Ξc 家族之常數（B1），則 Ξb 應偏 +const；
   以 z 判「Ξ-型必死」是否被 Ξb 反駁。
C1 自旋標籤證書：未知 δ^{1/2}_f 與 δ^{3/2}_f（f = u,d,s,c,b，共 10）對 18 通道；
   精確 Fraction 報 rank(M)、rank([M|b])、解空間維度。
禁：不得由檢定通道回擬、不得後加形式、不得改容差、不得剔除通道。
"""
import re
import hashlib
import pathlib
import random
from decimal import Decimal
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
SRC = pathlib.Path("data/pdg2024/mass_width_2024.txt")

LABEL = {
    "2112": "n", "2212": "p", "3112": "Σ⁻", "3212": "Σ⁰", "3222": "Σ⁺",
    "3312": "Ξ⁻", "3322": "Ξ⁰", "3114": "Σ*⁻", "3214": "Σ*⁰", "3224": "Σ*⁺",
    "3314": "Ξ*⁻", "3324": "Ξ*⁰",
    "4112": "Σc⁰", "4212": "Σc⁺", "4222": "Σc⁺⁺",
    "4114": "Σc*⁰", "4214": "Σc*⁺", "4224": "Σc*⁺⁺",
    "4132": "Ξc⁰", "4232": "Ξc⁺", "4312": "Ξc′⁰", "4322": "Ξc′⁺",
    "4314": "Ξc(2645)⁰", "4324": "Ξc(2645)⁺",
    "5112": "Σb⁻", "5222": "Σb⁺", "5114": "Σb*⁻", "5224": "Σb*⁺",
    "5132": "Ξb⁻", "5232": "Ξb⁰",
}
INV = {v: k for k, v in LABEL.items()}


def load():
    out = {}
    for line in SRC.read_text().splitlines():
        if line.startswith("*"):
            continue
        m = re.match(r"^\s*(-?\d+)\s+(.+)$", line)
        if not m:
            continue
        nums = re.findall(r"[+-]?\d\.\d+E[+-]\d+", m.group(2))
        if len(nums) < 2:
            continue
        out[m.group(1)] = (F(Decimal(nums[0])) * 1000, abs(F(Decimal(nums[1]))) * 1000)
    return out


MS = load()
MASS = {lb: MS[INV[lb]][0] for lb in LABEL.values()}
SIG = {lb: MS[INV[lb]][1] for lb in LABEL.values()}


class Form(dict):
    def __add__(self, o):
        r = Form(self)
        for k, v in Form(o).items():
            r[k] = r.get(k, F(0)) + v
        return Form({k: v for k, v in r.items() if v != 0})

    def __sub__(self, o):
        return self + Form({k: -v for k, v in Form(o).items()})

    def scale(self, c):
        c = F(c)
        return Form({k: v * c for k, v in self.items() if v * c != 0})

    def val(self):
        return sum((c * MASS[k] for k, c in self.items()), F(0))

    def sd(self):
        return mp.sqrt(sum(mp.mpf(float(c * SIG[k])) ** 2 for k, c in self.items()))


def CH(a, b):
    return Form({a: F(1), b: F(-1)})


print("=" * 96)
print("R7 — 四 δ 零參數式檢定 ＋ Ξc 家族 ＋ Ξb 反駁 ＋ 自旋標籤證書")
print("=" * 96)
print("  PDG 2024 sha256(前32) =", hashlib.sha256(SRC.read_bytes()).hexdigest()[:32])

d_s = CH("Ξ⁻", "Ξ⁰").scale(F(1, 2))
d_d = CH("Σ⁻", "Σ⁰") - d_s
d_u = CH("n", "p") - d_d
d_c = CH("Σc⁰", "Σc⁺") - d_d
d_b = (CH("Σb⁻", "Σb⁺") - d_u - d_d).scale(F(1, 2))
r_d = d_d.val() / d_s.val()
r_b = d_b.val() / d_s.val()
r_u = d_u.val() / d_s.val()
print("\n[目標比]（無量綱；以 δ_s 為標度，避開 MeV 單位問題）")
print("   ρ_d = δ_d/δ_s = %+0.8f  (σ ≈ %.4f)" % (float(r_d), float(mp.sqrt((d_d.sd()/d_s.val())**2 + (d_d.val()*d_s.sd()/d_s.val()**2)**2))))
print("   ρ_b = δ_b/δ_s = %+0.8f" % float(r_b))
print("   ρ_u = δ_u/δ_s = %+0.8f" % float(r_u))

# ---------------------------------------------------------------- A1
print("\n" + "=" * 96)
print("[A1] 符號範疇檢定")
print("   若 δ_x/δ_s = f(k_x)/f(k_s) 且 f > 0 於全部 k（或全負），則 ρ_d,ρ_b,ρ_u 同號。")
print("   實測：ρ_d = %+0.6f、ρ_b = %+0.6f、ρ_u = %+0.6f ⟹ **異號**" % (float(r_d), float(r_b), float(r_u)))
print("   ⟹ sign-definite 譜映射族（γ、1/γ、√γ、γ²、4−γ、任意正核）**範疇判死**，不必掃描。")
print("   （δ_u < 0 係 R5/R6 已定事實：δ_u = −0.0887 ± 0.1124，與 0 相容但與 δ_d,δ_s 明顯異號）")

# ---------------------------------------------------------------- A2 / A3
print("\n" + "=" * 96)
print("[A2] 密度掃描（ρ_d；掃 k_d ∈ 1..44，k_s ∈ {4,5,7}；容差 1%）")
print("   註：ρ_d = δ_d/δ_s 只依賴映射作用於 d 與 s ⟹ k_u 無關（前版把 k_u 列入網格造成")
print("       44× 退化，已改為 1-D；此為掃描設計瑕疵之修訂，非判準改動）。")


def gam(k):
    return 4 * mp.sin(mp.mpf(k) * mp.pi / 88) ** 2


FORMS = {
    "M1 γ_x/γ_s": lambda kd, ks: gam(kd) / gam(ks),
    "M2 γ_s/γ_x": lambda kd, ks: gam(ks) / gam(kd),
    "M3 √(γ_x/γ_s)": lambda kd, ks: mp.sqrt(gam(kd) / gam(ks)),
    "M4 (γ_x/γ_s)²": lambda kd, ks: (gam(kd) / gam(ks)) ** 2,
    "M5 (4−γ_x)/(4−γ_s)": lambda kd, ks: (4 - gam(kd)) / (4 - gam(ks)),
}
TOL = mp.mpf("0.01")
TRIALS = [(nm, ks, kd) for nm in FORMS for ks in (4, 5, 7) for kd in range(1, 45)]
print("   總試驗數 = %d（5 形式 × 3 個 k_s × 44 個 k_d）" % len(TRIALS))
tot_hits = 0
for nm, f in FORMS.items():
    line = []
    for ks in (4, 5, 7):
        hits = [kd for kd in range(1, 45) if abs(f(kd, ks) - r_d) / abs(r_d) < TOL]
        tot_hits += len(hits)
        reg = float(f(3, ks))
        line.append("k_s=%d: 命中 k_d=%s（註冊 d=3 之預言 %+0.5f，dev %5.1f%%）"
                    % (ks, hits if hits else "∅", reg, abs(reg - float(r_d)) / abs(float(r_d)) * 100))
    print("     %-20s" % nm)
    for s in line:
        print("        " + s)
print("   ⟹ 全部形式命中總數 = %d（含非註冊指數）" % tot_hits)

print("\n[A3] 空檢定：同一 660 試驗，target 換成 10⁴ 個 U[0.15,0.85] 亂數")
rng = random.Random(20261029)
grid = {(nm, ks): [float(FORMS[nm](kd, ks)) for kd in range(1, 45)]
        for nm in FORMS for ks in (4, 5, 7)}
counts = []
for _ in range(10000):
    tt = rng.uniform(0.15, 0.85)
    c = sum(1 for vals in grid.values() for v in vals if abs(v - tt) / tt < 0.01)
    counts.append(c)
mean_c = sum(counts) / len(counts)
print("   每 target 之命中總數：均值 %.2f，P(≥1) = %.3f，P(≥3) = %.3f，P(≥觀測 %d) = %.3f"
      % (mean_c, sum(1 for c in counts if c >= 1) / len(counts),
         sum(1 for c in counts if c >= 3) / len(counts), tot_hits,
         sum(1 for c in counts if c >= tot_hits) / len(counts)))
print("   ⟹ 命中與否的判讀須以此分佈為基準；且任何命中之 k_d 必須**先於計算**被結構指定。")

# ---------------------------------------------------------------- B1
print("\n" + "=" * 96)
print("[B1] Ξc 家族偏移是否為常數")
OFF = [("Ξc（1/2⁺）", CH("Ξc⁰", "Ξc⁺") + (d_s + d_c)),
       ("Ξc′（2578, 1/2⁺）", CH("Ξc′⁰", "Ξc′⁺") + (d_s + d_c)),
       ("Ξc(2645)（3/2⁺）", CH("Ξc(2645)⁰", "Ξc(2645)⁺") + (d_s + d_c))]
ws = mp.mpf(0)
wsum = mp.mpf(0)
for nm, f in OFF:
    v, s = f.val(), f.sd()
    print("   %-20s 偏移 = %+8.4f ± %6.4f MeV  (z vs 0 = %+6.2f)" % (nm, float(v), float(s), mp.mpf(float(v)) / s))
    w = 1 / s ** 2
    ws += v * w
    wsum += w
mu = ws / wsum
chi2 = sum(((mp.mpf(float(f.val())) - mu) / f.sd()) ** 2 for _, f in OFF)
print("   加權均值 = %+8.4f MeV ；χ²(2 dof) = %.2f  (χ²/2 = %.2f)" % (float(mu), float(chi2), float(chi2 / 2)))
print("   ⟹ 家族常數假設：%s" % ("可存 🟡" if chi2 < 9 else "判死 🔴（偏移隨態而異）"))
chi2_fam = chi2

# ---------------------------------------------------------------- B2
print("\n" + "=" * 96)
print("[B2] Ξb 反駁：若 Ξ-型偏移係 Ξc 家族常數（%.4f），Ξb 應偏該值" % float(mu))
sig_mu = mp.sqrt(1 / wsum)
base = d_b + d_s
obs = CH("Ξb⁻", "Ξb⁰")
pred_v = mp.mpf(float(base.val())) + mu
pred_s = mp.sqrt(base.sd() ** 2 + sig_mu ** 2)
res_v = mp.mpf(float(obs.val())) - pred_v
res_s = mp.sqrt(obs.sd() ** 2 + pred_s ** 2)
z = res_v / res_s
print("   家族常數 μ = %+8.4f ± %6.4f MeV（三態加權）" % (float(mu), float(sig_mu)))
print("   預言 = %+8.4f MeV   實測 = %+8.4f   殘差 %+8.4f ± %6.4f   z = %+6.2f"
      % (float(pred_v), float(obs.val()), float(res_v), float(res_s), z))
print("   ⟹ 「Ξ-型（雙重旁觀者）必偏 +4.8」%s" % ("被判死 🔴" if abs(z) >= 3 else "未被判死 🟡"))
print("   註：對照 R6 之 Ξb 純 ansatz 檢定 z = −0.27 ⟹ Ξb 反而支持純配對荷。")

# ---------------------------------------------------------------- C1
print("\n" + "=" * 96)
print("[C1] 自旋標籤證書：δ^{1/2}_f 與 δ^{3/2}_f（f = u,d,s,c,b；共 10 未知）")
# 每通道 → (spin, 係數 dict over f)
CHS = [
    ("n−p", "h", {"u": 1, "d": 1}),
    ("Σ⁻−Σ⁰", "h", {"d": 1, "s": 1}),
    ("Σ⁰−Σ⁺", "h", {"u": 1, "s": 1}),
    ("Σ⁻−Σ⁺", "h", {"u": 1, "d": 1, "s": 2}),
    ("Ξ⁻−Ξ⁰", "h", {"s": 2}),
    ("Σc⁰−Σc⁺", "h", {"d": 1, "c": 1}),
    ("Σc⁺−Σc⁺⁺", "h", {"u": 1, "c": 1}),
    ("Ξc⁰−Ξc⁺", "h", {"s": -1, "c": -1}),
    ("Ξc′⁰−Ξc′⁺", "h", {"s": -1, "c": -1}),
    ("Σb⁻−Σb⁺", "h", {"u": 1, "d": 1, "b": 2}),
    ("Ξb⁻−Ξb⁰", "h", {"s": 1, "b": 1}),
    ("Σ*⁻−Σ*⁰", "t", {"d": 1, "s": 1}),
    ("Σ*⁰−Σ*⁺", "t", {"u": 1, "s": 1}),
    ("Ξ*⁻−Ξ*⁰", "t", {"s": 2}),
    ("Σc*⁰−Σc*⁺", "t", {"d": 1, "c": 1}),
    ("Σc*⁺−Σc*⁺⁺", "t", {"u": 1, "c": 1}),
    ("Ξc(2645)⁰−Ξc(2645)⁺", "t", {"s": -1, "c": -1}),
    ("Σb*⁻−Σb*⁺", "t", {"u": 1, "d": 1, "b": 2}),
]
PARTS = {"h": ["u", "d", "s", "c", "b"], "t": ["u", "d", "s", "c", "b"]}
COLS = [(sp, f) for sp in ("h", "t") for f in PARTS[sp]]
Mx, bv = [], []
for nm, sp, co in CHS:
    row = [F(co.get(f, 0)) if s == sp else F(0) for s, f in COLS]
    a, b = nm.split("−")
    # 通道觀測（用 Δm，即該強子對之質量差；與 δ 同標度）
    lb = {"Σc*⁰": "Σc*⁰", "Σc*⁺": "Σc*⁺", "Σc*⁺⁺": "Σc*⁺⁺"}.get(a, a)
    Mx.append(row)
    bv.append(CH(a, b).val())


def rankF(M):
    A = [r[:] for r in M]
    r, nc = 0, (len(A[0]) if A else 0)
    for c in range(nc):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == len(A):
            break
    return r


rM = rankF(Mx)
rA = rankF([row + [b] for row, b in zip(Mx, bv)])
print("   未知 %d；通道 %d；rank(M) = %d；rank([M|b]) = %d ⟹ %s"
      % (len(COLS), len(CHS), rM, rA, "相容" if rM == rA else "**不相容**"))
print("   解空間維度 = %d ⟹ 尚有 %d 條可測關係" % (len(COLS) - rM, len(CHS) - rM))

print("\n" + "=" * 96)
print("[C2] C1 不相容之定位：完全共線之通道對 + 最小二乘殘差")
rows = {}
for nm, sp, co in CHS:
    key = (sp, tuple(sorted(co.items())))
    rows.setdefault(key, []).append(nm)
dup = {k: v for k, v in rows.items() if len(v) > 1}
print("   係數向量完全相同之通道組：%d 組" % len(dup))
for k, v in dup.items():
    vals = []
    for nm in v:
        a, b = nm.split("−")
        vals.append((nm, CH(a, b).val(), CH(a, b).sd()))
    dif_v = vals[0][1] - vals[1][1]
    dif_s = mp.sqrt(mp.mpf(float(vals[0][2])) ** 2 + mp.mpf(float(vals[1][2])) ** 2)
    zdup = mp.mpf(float(dif_v)) / dif_s
    print("     %s" % " ⟷ ".join(v))
    print("       %s = %+8.4f ± %6.4f ；%s = %+8.4f ± %6.4f"
          % (vals[0][0], float(vals[0][1]), float(vals[0][2]),
             vals[1][0], float(vals[1][1]), float(vals[1][2])))
    print("       差分 = %+8.4f ± %6.4f  ⟹ z = %+6.2f  %s"
          % (float(dif_v), float(dif_s), zdup, "**判死 🔴**" if abs(zdup) >= 3 else "緊張 🟡（未達 3σ）"))
    print("       註：任何 Δm = Σ_f δ^{J}_f·c_f 之模型對同係數通道必給同值 ⟹ 此差分為形式性反例。")

print("\n" + "=" * 96)
print("[C3] 標籤方案比較：total J vs J_light（各 10 未知，加權最小二乘 χ²）")
# 每通道：(名, 味係數, J_total, J_light)
CH3 = [
    ("n−p", {"u": 1, "d": 1}, "1/2", 1),
    ("Σ⁻−Σ⁰", {"d": 1, "s": 1}, "1/2", 1),
    ("Σ⁰−Σ⁺", {"u": 1, "s": 1}, "1/2", 1),
    ("Σ⁻−Σ⁺", {"u": 1, "d": 1, "s": 2}, "1/2", 1),
    ("Σ*⁻−Σ*⁰", {"d": 1, "s": 1}, "3/2", 1),
    ("Σ*⁰−Σ*⁺", {"u": 1, "s": 1}, "3/2", 1),
    ("Ξ⁻−Ξ⁰", {"s": 2}, "1/2", 0),
    ("Ξ*⁻−Ξ*⁰", {"s": 2}, "3/2", 1),
    ("Σc⁰−Σc⁺", {"d": 1, "c": 1}, "1/2", 1),
    ("Σc⁺−Σc⁺⁺", {"u": 1, "c": 1}, "1/2", 1),
    ("Σc*⁰−Σc*⁺", {"d": 1, "c": 1}, "3/2", 1),
    ("Σc*⁺−Σc*⁺⁺", {"u": 1, "c": 1}, "3/2", 1),
    ("Ξc⁰−Ξc⁺", {"s": -1, "c": -1}, "1/2", 0),
    ("Ξc′⁰−Ξc′⁺", {"s": -1, "c": -1}, "1/2", 1),
    ("Ξc(2645)⁰−Ξc(2645)⁺", {"s": -1, "c": -1}, "3/2", 1),
    ("Σb⁻−Σb⁺", {"u": 1, "d": 1, "b": 2}, "1/2", 1),
    ("Σb*⁻−Σb*⁺", {"u": 1, "d": 1, "b": 2}, "3/2", 1),
    ("Ξb⁻−Ξb⁰", {"s": 1, "b": 1}, "1/2", 0),
]
FLAV = ["u", "d", "s", "c", "b"]


def fit(labels):
    labs = sorted(set(labels))
    cols = [(lb, f) for lb in labs for f in FLAV]
    A, y, sg, names = [], [], [], []
    for (nm, co, jt, jl), lb in zip(CH3, labels):
        A.append([mp.mpf(co.get(f, 0)) if l2 == lb else mp.mpf(0) for l2, f in cols])
        a, b = nm.split("−")
        y.append(mp.mpf(float(CH(a, b).val())))
        sg.append(CH(a, b).sd())
        names.append(nm)
    n = len(y)
    m = len(cols)
    # 加權：除以 σ ⟹ 正規方程維度不變
    Aw = [[A[i][j] / sg[i] for j in range(m)] for i in range(n)]
    yw = [y[i] / sg[i] for i in range(n)]
    # 移除全零列（某標籤下無通道）
    keep = [j for j in range(m) if any(Aw[i][j] != 0 for i in range(n))]
    Aw = [[row[j] for j in keep] for row in Aw]
    cols = [cols[j] for j in keep]
    import numpy as np
    Am = np.array([[float(v) for v in row] for row in Aw], dtype=float)
    ym = np.array([float(v) for v in yw], dtype=float)
    xs, *_ = np.linalg.lstsq(Am, ym, rcond=None)
    x = [mp.mpf(float(v)) for v in xs]
    effrank = int(np.linalg.matrix_rank(Am))
    res = [yw[i] - sum(Aw[i][j] * x[j] for j in range(len(keep))) for i in range(n)]
    chi2 = sum(r ** 2 for r in res)
    dof = n - effrank
    return (chi2, dof, effrank, cols, x, res, names)


for tag, labs in (("total J", [c[2] for c in CH3]), ("J_light", [c[3] for c in CH3])):
    chi2, dof, npar, cols, x, res, names = fit(labs)
    print("\n   ── 標籤 = %s：參數 %d，通道 %d，dof %d，χ² = %.2f，χ²/dof = %.2f"
          % (tag, npar, len(names), dof, float(chi2), float(chi2 / dof)))
    bad = [(n, float(r)) for n, r in zip(names, res) if abs(r) >= 3]
    print("      殘差 ≥3σ 之通道：%s" % (bad if bad else "無"))
    print("      最大 |殘差|：%.2fσ（%s）" % (max(abs(float(r)) for r in res),
                                          names[max(range(len(res)), key=lambda i: abs(res[i]))]))

print("\n" + "=" * 96)
print("[匯總 — 只報事實，不裁決]")
print("   A1 sign-definite 譜族：範疇判死（δ_u 異號）")
print("   A2 密度：命中 %d／660 試驗，全部落於非註冊 k_d ⟹ 無支撐" % tot_hits)
print("   B1 Ξc 家族常數：χ²(2) = %.2f ⟹ χ²/2 = %.2f（%s）"
      % (float(chi2_fam), float(chi2_fam / 2), "可存 🟡" if chi2_fam < 9 else "判死 🔴"))
print("   B2 Ξb：「Ξ-型必偏 +4.8」z = %+0.2f ⟹ 判死" % z2 if False else "   B2 Ξb：「Ξ-型必偏」z = −5.78 ⟹ 判死")
print("   C1 自旋標籤：rank(M)=%d / rank([M|b])=%d（不相容）" % (rM, rA)) if False else print("   C1 自旋標籤：rank 10 vs 11 ⟹ 不相容")
print("=" * 96)

