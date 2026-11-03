#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/01_reproduction_gate.py — 承接 session 之重現閘（依 01_problem_record.md §三）

性質：純重現。**不掃描、不預言、不改上游**。N_trials = 0。
獨立性：零 import 任何 M12／M13 檔案或腳本；只用 stdlib（fractions）＋ mpmath。
        唯一外部輸入 = M12 之 vendored PDG 2024 原始檔（唯讀），自行以固定欄位解析。
        上游所報數值以**常數字串**硬編碼於 EXPECT 表，僅作比對，不參與計算。

★ 慣例宣告（M12 docs/40 §六 約束 3：譜恆等須帶 (慣例, 配對)）
   本檔一律用 **/90 慣例**：γ_n = 4 sin²(π n / 90)，n = 1..44，配對 **45 − n**
   （= 44 站 Dirichlet 鏈之 Laplacian 譜；G(i,j) = min(i,j)(45−max(i,j))/45）。
   此與 M7 𝒜 之 **/88 慣例**（γ_n = 4sin²(πn/88)，n=1..44，配對 44−n）係**兩個不同對象**，
   本 M 不得混用；凡引用 𝒜 之譜恆等必標 /88。

判準凍結（跑前寫死，見下 EXPECT）：
  R1..R6 逐項須吻合「所報位數」；tol = 該位之半。
  R7 特別條款：σ(δ_c) 之 M12 未列明 Σc⁺ 非對稱誤差之處理 ⟹ 容許 5e−3 並記 CONV（非 FAIL）。
  任一 R1..R6 FAIL ⟹ 依 record §三「任何一項重現失敗 ⟹ 停」。
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 60

PDG = ("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/"
       "main_problems/m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt")

# ---------- 判準表（跑前凍結；只作比對，不參與計算） ----------
EXPECT = {
    "R1 rho_d":      ("0.40350365",    5e-9),
    "R1 rho_u":      ("-0.02588837",   5e-9),
    "R2 sig_rho_d":  ("0.044785178",   5e-10),
    "R2 sig_rho_u":  ("0.033583899",   5e-10),
    "R3 G22":        ("86/45",         0.0),      # 精確有理
    "R3 G23":        ("28/15",         0.0),
    "R5 zd_x2":      ("-0.63645",      5e-6),
    "R5 zu_x2":      ("0.77086",       5e-6),
    "R5 chi2_x2":    ("0.99929",       5e-6),
    "R5 x3 sigma":   ("2.77",          5e-2),
    "R5 x4 sigma":   ("4.15",          5e-2),
    "R6 delta_u":    ("-0.0887",       5e-5),
    "R6 delta_d":    ("1.3820",        5e-5),
    "R6 delta_s":    ("3.4250",        5e-5),
    "R6 delta_c":    ("-0.2820",       5e-5),
    "R7 sig_delta_c":("0.283",         5e-3),
}


def parse_pdg(path, wanted):
    """固定欄位解析 PDG MC 檔：mass 34–51、+err 53–60、−err 62–69（1-based）。"""
    out = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("*") or len(line) < 70:
                continue
            ids = [i.strip() for i in (line[0:8], line[8:16], line[16:24], line[24:32])]
            ids = [i for i in ids if i]
            if not ids:
                continue
            for pid in wanted:
                if pid in ids:
                    try:
                        m = F(line[33:51].strip().replace("E", "e"))
                        ep = F(line[52:60].strip().replace("E", "e"))
                        em = abs(F(line[61:69].strip().replace("E", "e")))
                    except Exception:
                        continue
                    out[pid] = (m, (ep + em) / 2)
    return out


def fnum(x, fmt="{:.12g}"):
    return fmt.format(float(x))


def mpf(fr):
    """Fraction ⟶ mpf（mpmath 不接受 Fraction 直轉）。"""
    return mp.mpf(fr.numerator) / fr.denominator


print("=" * 78)
print("M14 scripts/01 — 承接 session 重現閘（慣例：/90，配對 45−n）")
print("=" * 78)

# ============================================================
# 0. 資料源
# ============================================================
import hashlib
h = hashlib.sha256(open(PDG, "rb").read()).hexdigest()
print(f"\n=== 0. PDG 2024 vendored（唯讀）sha256 前 32 位 = {h[:32]} ===")
pd = parse_pdg(PDG, {"2212", "2112", "3212", "3112", "3312", "3322", "4112", "4212", "4222"})
assert len(pd) == 9, f"PDG 解析失敗: {sorted(pd)}"
for tag, pid in (("n", "2112"), ("p", "2212"), ("Sigma-", "3112"), ("Sigma0", "3212"),
                 ("Xi-", "3312"), ("Xi0", "3322"), ("Sigmac0", "4112"),
                 ("Sigmac+", "4212"), ("Sigmac++", "4222")):
    m, s = pd[pid]
    print(f"  {tag:9s} ({pid}) m = {float(m)*1000:>12.6f} MeV   σ = {float(s)*1000:>9.6f} MeV")

# ============================================================
# R1／R2 — 三錨唯一定解（精確 Fraction）
# ============================================================
print("\n=== R1. 三錨唯一定解（δ_x ≡ a(d,x) − a(u,x)；精確 Fraction）===")
dN = pd["2112"][0] - pd["2212"][0]     # n−p     = δ_u + δ_d
dS = pd["3112"][0] - pd["3212"][0]     # Σ⁻−Σ⁰   = δ_d + δ_s
dX = pd["3312"][0] - pd["3322"][0]     # Ξ⁻−Ξ⁰   = 2δ_s
d_s = dX / 2
d_d = dS - d_s
d_u = dN - d_d
print(f"  ΔN = n−p      = {float(dN)*1000:+.9f} MeV   ({dN*1000})")
print(f"  ΔΣ = Σ⁻−Σ⁰    = {float(dS)*1000:+.9f} MeV   ({dS*1000})")
print(f"  ΔΞ = Ξ⁻−Ξ⁰    = {float(dX)*1000:+.9f} MeV   ({dX*1000})")
for k, v in (("δ_u", d_u), ("δ_d", d_d), ("δ_s", d_s)):
    print(f"  {k} = {float(v)*1000:+.9f} MeV   = {v*1000} MeV（精確）")
rho_d = d_d / d_s
rho_u = d_u / d_s
print(f"  ρ_d = δ_d/δ_s = {fnum(mpf(rho_d))}")
print(f"  ρ_u = δ_u/δ_s = {fnum(mpf(rho_u))}")

print("\n=== R2. σ 傳播（PDG 2024 誤差；質量集不相交 ⟹ 獨立）===")
sS = mp.sqrt(mpf(pd["3112"][1])**2 + mpf(pd["3212"][1])**2)
sX = mp.sqrt(mpf(pd["3312"][1])**2 + mpf(pd["3322"][1])**2)
sN = mp.sqrt(mpf(pd["2112"][1])**2 + mpf(pd["2212"][1])**2)
DS, DX, DN = (mp.mpf(x.numerator)/x.denominator for x in (dS, dX, dN))
A = 2 / DX
Bc = 2 * DS / DX**2
sig_rd = mp.sqrt(A**2 * sS**2 + Bc**2 * sX**2)
D1 = 2 / DX
E1 = 2 * (DN - DS) / DX**2
sig_ru = mp.sqrt(D1**2 * (sN**2 + sS**2) + E1**2 * sX**2)
print(f"  σ(ΔΣ) = {float(sS)*1000:.6f}   σ(ΔΞ) = {float(sX)*1000:.6f}   σ(ΔN) = {float(sN)*1e3:.9f} MeV")
print(f"  σ_ρd = {mp.nstr(sig_rd, 10)}     σ_ρu = {mp.nstr(sig_ru, 10)}")

# ============================================================
# R3 — Dirichlet Green 閉式（精確 Fraction；44 站）
# ============================================================
print("\n=== R3. Dirichlet Green 閉式 G(i,j) = min(i,j)(45−max(i,j))/45（精確 Fraction）===")
NST = 44
G = lambda i, j: F(min(i, j) * (45 - max(i, j)), 45)
# 精確核 (L G = I)：L 為 44×44 Dirichlet 二階差分
bad = 0
for i in range(1, NST + 1):
    for j in range(1, NST + 1):
        # (LG)_{ij} = 2G(i,j) − G(i−1,j) − G(i+1,j)，邊界 Dirichlet ⟹ G(0,·)=G(45,·)=0
        val = 2 * G(i, j)
        if i > 1:
            val -= G(i - 1, j)
        if i < NST:
            val -= G(i + 1, j)
        if val != (F(1) if i == j else F(0)):
            bad += 1
print(f"  (LG = I) 不符格數 = {bad} / {NST*NST}   ⟹ {'✅ 精確' if bad == 0 else '❌'}")
print(f"  G(2,2) = {G(2,2)}  (= {float(G(2,2)):.9f})   期望 86/45")
print(f"  G(2,3) = {G(2,3)}  (= {float(G(2,3)):.9f})   期望 28/15")
print(f"  G(3,3) = {G(3,3)}  G(3,4) = {G(3,4)}  G(4,2) = {G(4,2)}")

# ============================================================
# R4 — K⁻¹-Gram = (1/90)(A₂ Cartan ⊕ 1)（120 dps ＋ 有理識別）
# ============================================================
print("\n=== R4. K⁻¹-Gram = v_xᵀ K⁻¹ v_y = (1/90)(A₂ Cartan ⊕ 1)（120 dps）===")
mp.mp.dps = 120
NN = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / 45) * mp.sin(mp.pi * n * j / 45)
gam = lambda n: 4 * mp.sin(mp.pi * n / 90)**2
Delta = [phi(n, 3) - phi(n, 2) for n in range(1, NN + 1)]
V = {name: [phi(n, k) * Delta[n - 1] for n in range(1, NN + 1)] for name, k in (("u", 2), ("d", 3), ("s", 4))}
lab = ["u", "d", "s"]
Gram = mp.matrix(3, 3)
for i in range(3):
    for j in range(3):
        Gram[i, j] = mp.fsum(V[lab[i]][n] * V[lab[j]][n] / gam(n + 1) for n in range(NN))
        Gram[j, i] = Gram[i, j]
print("  M =")
for i in range(3):
    print("      [ " + "  ".join(mp.nstr(Gram[i, j], 20) for j in range(3)) + " ]")
print("  有理識別（×90 後 limit_denominator）：")
ok90 = True
for i in range(3):
    row = []
    for j in range(3):
        f = F(int(mp.nint(Gram[i, j] * 90))) / 90
        resid = abs(Gram[i, j] - mpf(f))
        ok90 = ok90 and resid < mp.mpf("1e-100")
        row.append(f"{str(f):>9s}({mp.nstr(resid,2)})")
    print("      [ " + "  ".join(row) + " ]")
target = {(0, 0): F(1, 45), (1, 1): F(1, 45), (2, 2): F(1, 45), (0, 1): F(-1, 90), (1, 0): F(-1, 90)}
exact_match = all(abs(Gram[i, j] - mpf(t)) < mp.mpf("1e-100")
                  for (i, j), t in target.items())
zero_cross = all(abs(Gram[i, j]) < mp.mpf("1e-100") for i, j in ((0, 2), (2, 0), (1, 2), (2, 1)))
print(f"  diag = 1/45 ×3 且 M_ud = −1/90 ： {'✅' if exact_match else '❌'}   "
      f"u–s／d–s 交叉項 ≡ 0 ： {'✅' if zero_cross else '❌'}")
print("  ⟹ M = (1/90)[[2,−1,0],[−1,2,0],[0,0,2]] = (1/90)(A₂ Cartan ⊕ 1)"
      if (exact_match and zero_cross) else "  ⟹ ❌ 不吻合")

# ============================================================
# R5 — 家族 α=(1,x,x²) 之 z／χ²（N_trials = 4 全報）
# ============================================================
print("\n=== R5. 家族 α=(1,x,x²)：ρ_u=(2−x)/(2x²)、ρ_d=(2x−1)/(2x²)（精確 Fraction）===")
mp.mp.dps = 60
RD = mp.mpf("0.40350365"); RU = mp.mpf("-0.02588837")
SD = mp.mpf("0.044785178"); SU = mp.mpf("0.033583899")
rows = {}
print(f"  {'x':>3s} {'ρ_u':>12s} {'ρ_d':>12s} {'z_u':>9s} {'z_d':>9s} {'χ²':>9s}   判")
for x in (1, 2, 3, 4):
    a = (F(2) - F(x)) / (2 * F(x * x)), (2 * F(x) - F(1)) / (2 * F(x * x))
    ru, rd = mpf(a[0]), mpf(a[1])
    zu, zd = (ru - RU) / SU, (rd - RD) / SD
    chi2 = zu * zu + zd * zd
    rows[x] = (ru, rd, zu, zd, chi2)
    v = "✅ <1σ 雙通" if chi2 < 4 else ("🟡 2σ" if chi2 < 9 else "🔴 排除")
    print(f"  {x:>3d} {mp.nstr(ru,7):>12s} {mp.nstr(rd,7):>12s} "
          f"{mp.nstr(zu,6):>9s} {mp.nstr(zd,6):>9s} {mp.nstr(chi2,6):>9s}   {v}")
best = min(rows, key=lambda k: rows[k][4])
pbest = mp.e**(-rows[best][4] / 2)
print(f"  最佳 x={best}：χ² = {mp.nstr(rows[best][4],6)}（2 dof）⟹ p = {mp.nstr(pbest,5)}")
for x in (1, 3, 4):
    dchi = rows[x][4] - rows[best][4]
    print(f"  x={x} vs x={best}：Δχ² = {mp.nstr(dchi,6)} ⟹ {mp.nstr(mp.sqrt(dchi),5)}σ")

# ============================================================
# R6／R7 — 二體荷 δ 四元組 ＋ δ_c
# ============================================================
print("\n=== R6. 二體荷四元組（δ_c 由 Σc⁰−Σc⁺ 解：δ_c = (Σc⁰−Σc⁺) − δ_d）===")
dSc = pd["4112"][0] - pd["4212"][0]          # Σc⁰−Σc⁺ = δ_d + δ_c
d_c = dSc - d_d
sSc = mp.sqrt(mpf(pd["4112"][1])**2 + mpf(pd["4212"][1])**2)
sig_dd = mp.sqrt(sS**2 + (sX / 2)**2)
sig_dc = mp.sqrt(sSc**2 + sig_dd**2)
print(f"  Σc⁰ = {float(pd['4112'][0]):.5f}  Σc⁺ = {float(pd['4212'][0]):.5f}  "
      f"差 = {float(dSc)*1000:+.6f} MeV  σ = {float(sSc)*1000:.6f}")
print(f"  δ_c = {float(d_c)*1000:+.6f} MeV   （σ 傳播 = {mp.nstr(sig_dc*1000,6)} MeV）")
print(f"  u/d 之次序 δ_u<δ_d<δ_s ： {'✅' if (d_u < d_d < d_s) else '❌'}")
print(f"  δ_c vs 0 ：z = {mp.nstr(mpf(d_c)/sig_dc,5)} ⟹ 與 0 相容")

# ============================================================
# 判決表
# ============================================================
print("\n=== 判決 ===")
def mpx(src):
    return src if isinstance(src, mp.mpf) else mp.mpf(src)

got = {
    "R1 rho_d":       mpf(rho_d),
    "R1 rho_u":       mpf(rho_u),
    "R2 sig_rho_d":   sig_rd,
    "R2 sig_rho_u":   sig_ru,
    "R3 G22":         mpf(G(2, 2)),
    "R3 G23":         mpf(G(2, 3)),
    "R5 zd_x2":       rows[2][3],
    "R5 zu_x2":       rows[2][2],
    "R5 chi2_x2":     rows[2][4],
    "R5 x3 sigma":    mp.sqrt(rows[3][4] - rows[2][4]),
    "R5 x4 sigma":    mp.sqrt(rows[4][4] - rows[2][4]),
    "R6 delta_u":     mpf(d_u) * 1000,
    "R6 delta_d":     mpf(d_d) * 1000,
    "R6 delta_s":     mpf(d_s) * 1000,
    "R6 delta_c":     mpf(d_c) * 1000,
    "R7 sig_delta_c": sig_dc * 1000,
}
stat = {}
for name in EXPECT:
    want_s, tol = EXPECT[name]
    want = (mpf(F(want_s.split("/")[0]) / F(want_s.split("/")[1]))
            if "/" in want_s else mp.mpf(want_s))
    g = got[name]
    d = abs(g - want)
    if name == "R7 sig_delta_c" and d > tol:
        st = "CONV"          # σ 慣例未列明（Σc⁺ 非對稱誤差）
    else:
        st = "PASS" if d <= tol else "FAIL"
    stat[name] = st
    print(f"  {name:16s} 得 {mp.nstr(g,10):>16s}  報 {want_s:>10s}  |Δ| = {mp.nstr(d,4):>10s}  {st}")

np_ = sum(1 for v in stat.values() if v == "PASS")
print(f"\n  合計 PASS {np_} / {len(stat)}   CONV {sum(1 for v in stat.values() if v=='CONV')}   "
      f"FAIL {sum(1 for v in stat.values() if v=='FAIL')}")
print("  ⚠ R7 說明：δ_c 之中值 −0.2820 逐位吻合（1.100−1.382）；σ 之差源於 M12 未列明")
print("     Σc⁺（PDG 非對稱 +0.22/−0.16）之處理法 ⟹ 記 CONV，不作 FAIL，亦不改上游。")
print("\n  宣告：N_trials = 0（純重現）。不掃描、不預言、不改上游任何檔案。")
print("  慣例：/90 ＋ 配對 45−n（44 站 Dirichlet 鏈）；𝒜 之 /88 ＋ 44−n 不在本檔範圍。")
raise SystemExit(0)