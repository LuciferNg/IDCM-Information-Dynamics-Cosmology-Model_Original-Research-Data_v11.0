#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/01_reproduction_gate.py — 重現閘（依 SESSION_START §重現閘）

性質：純重現，**不掃描、不預言**。N_trials = 0。
獨立性：零 import 任何 M12 檔案／腳本；只用 stdlib（fractions/math）＋ mpmath（任意精度）。
        唯一外部輸入 = M12 之 vendored PDG 2024 原始檔（唯讀 `data/pdg2024/mass_width_2024.txt`，不修改）。
        對照值（docs 所報）以常數字串硬編碼於 EXPECT 表，僅供比對，不參與計算。

判準凍結（跑前寫死）：
  R1..R7 逐項須吻合到「所報位數」（tol 見下表；tol = 最後一位之半）。
  任一項 FAIL ⟹ 依 §重現閘停止，不得開始路線圖。
  R6 特別條款：docs/38 §E 之 0.03510 已被 M12 自己之 R27（docs/42）更正為 0.04478518；
  若本腳本重算吻合 R27 而不吻合 docs/38，記 SUPERSEDED（不算 FAIL，但須 register note 回報），
  並以較嚴之 R27 值為正典。
"""
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 50

PDG = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/main_problems/m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt"

# ---------- 判準表（跑前凍結） ----------
EXPECT = {
    "R1 rho_d":        ("0.40350365",   5e-9),
    "R2 rho_u":        ("-0.02588837",  5e-9),
    "R3 ratio":        ("0.064158949",  1e-8),
    "R4 cos(vu,vd)":   ("-0.8000000",   5e-8),
    "R5 angle":        ("90.848631",    5e-7),
    "R6 sig/1_41":     ("1.4391",       5e-5),
    "R7 gamma pair":   ("4",            1e-12),
    "R7b cross dev":   ("0.139598",     5e-7),
}

def parse_pdg(path, wanted):
    """固定欄位解析 PDG MC 檔：mass 34-51, +err 53-60, -err 62-69。"""
    out = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("*") or len(line) < 70:
                continue
            ids = [line[0:8], line[8:16], line[16:24], line[24:32]]
            ids = [i.strip() for i in ids if i.strip()]
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

pd = parse_pdg(PDG, {"2212", "2112", "3112", "3212", "3312", "3322"})
assert len(pd) == 6, f"PDG 解析失敗: {sorted(pd)}"

sea = {}
def show(tag, m, s):
    sea[tag] = (m, s)
    print(f"  {tag:8s} m = {float(m)*1000:.6f} MeV   σ = {float(s)*1000:.6f} MeV")
print("=== 0. PDG 2024 原始值（唯讀 vendored 檔，獨立解析）===")
show("n",  *pd["2112"]); show("p",  *pd["2212"])
show("Sigma-", *pd["3112"]); show("Sigma0", *pd["3212"])
show("Xi-", *pd["3312"]);   show("Xi0", *pd["3322"])

# ---------- A. 三錨解（精確 Fraction） ----------
print("\n=== A. 三錨唯一定解（docs/14 §三 之定義，精確 Fraction）===")
dN = pd["2112"][0] - pd["2212"][0]          # n−p      = δ_u+δ_d
dS = pd["3112"][0] - pd["3212"][0]          # Σ⁻−Σ⁰    = δ_d+δ_s
dX = pd["3312"][0] - pd["3322"][0]          # Ξ⁻−Ξ⁰    = 2δ_s
d_s = dX / 2
d_d = dS - d_s
d_u = dN - d_d
for k, v in (("δ_u", d_u), ("δ_d", d_d), ("δ_s", d_s)):
    print(f"  {k} = {float(v)*1000:+.9f} MeV   (= {v*1000} MeV 精確)")
rho_d = d_d / d_s
rho_u = d_u / d_s
ratio = abs(rho_u) / rho_d
print(f"  ρ_d = δ_d/δ_s = {mp.nstr(mp.mpf(rho_d.numerator)/rho_d.denominator, 12)}")
print(f"  ρ_u = δ_u/δ_s = {mp.nstr(mp.mpf(rho_u.numerator)/rho_u.denominator, 12)}")
print(f"  |ρ_u|/ρ_d     = {mp.nstr(mp.mpf(ratio.numerator)/ratio.denominator, 12)}")
print(f"  閉式核: ρ_d = 2ΔΣ/ΔΞ − 1 = {mp.nstr(2*mp.mpf(dS.numerator)/dS.denominator/(mp.mpf(dX.numerator)/dX.denominator)-1, 12)}")
print(f"          ρ_u = 2(ΔN−ΔΣ)/ΔΞ + 1 = {mp.nstr(2*(mp.mpf(dN.numerator)/dN.denominator-mp.mpf(dS.numerator)/dS.denominator)/(mp.mpf(dX.numerator)/dX.denominator)+1, 12)}")

# ---------- B. σ 傳播（ρ_d = 2ΔΣ/ΔΞ − 1；ΔΣ、ΔΞ 質量集不相交 ⟹ 獨立） ----------
print("\n=== B. σ 傳播（獨立，PDG 2024 誤差；無協方差項）===")
sS = mp.sqrt(pd["3112"][1]**2 + pd["3212"][1]**2)       # σ(ΔΣ)
sX = mp.sqrt(pd["3312"][1]**2 + pd["3322"][1]**2)       # σ(ΔΞ)
sN = mp.sqrt(pd["2112"][1]**2 + pd["2212"][1]**2)
DS = mp.mpf(dS.numerator)/dS.denominator
DX = mp.mpf(dX.numerator)/dX.denominator
DN = mp.mpf(dN.numerator)/dN.denominator
A = 2*mp.mpf(1)/DX
Bc = 2*DS/DX**2
sig_rd = mp.sqrt(A**2*mp.mpf(sS)**2 + Bc**2*mp.mpf(sX)**2)
# ρ_u = 2(ΔN−ΔΣ)/ΔΞ + 1
D = 2/DX
E = 2*(DN-DS)/DX**2
sig_ru = mp.sqrt(D**2*(mp.mpf(sN)**2+mp.mpf(sS)**2) + E**2*mp.mpf(sX)**2)
print(f"  σ(ΔΣ)={float(sS)*1000:.6f}  σ(ΔΞ)={float(sX)*1000:.6f}  σ(ΔN)={float(sN)*1000:.9f} MeV")
print(f"  σ_ρd = {mp.nstr(sig_rd, 10)}      σ_ρu = {mp.nstr(sig_ru, 10)}")
print(f"  比值 σ_ρd ÷ (1/41) = {mp.nstr(sig_rd*41, 8)}")
print(f"  比值 σ_ρu ÷ (1/41) = {mp.nstr(sig_ru*41, 8)}")
print(f"  半格 1/82 = {mp.nstr(mp.mpf(1)/82, 8)}   <  σ ?  "
      f"{'YES（不可裁決）' if mp.mpf(1)/82 < sig_rd else 'NO'}")
print("  附：以 docs/08 之 σ(ΔΣ)=0.038419、σ(ΔΞ)=0.211896 重算 → "
      f"{mp.nstr(mp.sqrt((2/DX)**2*mp.mpf('0.038419')**2 + (2*DS/DX**2)**2*mp.mpf('0.211896')**2)*41, 8)}")

# ---------- C. 測試向量幾何（φ_n(j)=√(2/45)sin(πnj/45)，n=1..44） ----------
print("\n=== C. 測試向量 v_x 幾何（docs/37 §二）===")
N = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2)/45) * mp.sin(mp.pi*n*j/45)
Delta = [phi(n, 3) - phi(n, 2) for n in range(1, N+1)]
V = {}
for name, k in (("u", 2), ("d", 3), ("s", 4)):
    V[name] = [phi(n, k)*Delta[n-1] for n in range(1, N+1)]
def dot(a, b): return mp.fsum(x*y for x, y in zip(a, b))
def norm(a): return mp.sqrt(dot(a, a))
for k in "uds":
    print(f"  ‖v_{k}‖ = {mp.nstr(norm(V[k]), 12)}")
for a, b in (("u","d"), ("d","s"), ("u","s")):
    c = dot(V[a], V[b])/(norm(V[a])*norm(V[b]))
    print(f"  cos(v_{a},v_{b}) = {mp.nstr(c, 12)}   ({mp.nstr(mp.degrees(mp.acos(c)), 8)}°)")
G = mp.matrix(3, 3)
order = ["d", "u", "s"]                            # c 之次序與向量一致：ρ_d↔v_d, ρ_u↔v_u
for i in range(3):
    for j in range(3):
        G[i, j] = dot(V[order[i]], V[order[j]])
ev = mp.eigsy(G, eigvals_only=True)
print("  Gram 譜 =", [mp.nstr(x, 8) for x in sorted(ev)])

# ---------- D. 最小範數核 w*（3 約束、44 維） ----------
print("\n=== D. 最小範數核 w*（docs/37 §二 之達標核取向）===")
c = mp.matrix([mp.mpf(rho_d.numerator)/rho_d.denominator,
               mp.mpf(rho_u.numerator)/rho_u.denominator, mp.mpf(1)])
Ginv = G**-1
alpha = Ginv*c                                     # G⁻¹c
w = [mp.fsum(alpha[i]*V[order[i]][n] for i in range(3)) for n in range(N)]
nw = norm(w)
print(f"  ‖w*‖ = {mp.nstr(nw, 10)}")
dot_u = dot(w, V["u"])
ang = mp.degrees(mp.acos(dot_u/(nw*norm(V["u"]))))
print(f"  ⟨w*,v_u⟩ = {mp.nstr(dot_u, 10)}   ∠(w*,v_u) = {mp.nstr(ang, 10)}°")
print(f"  符號圖樣: 正 {sum(1 for x in w if x > 0)} / 負 {sum(1 for x in w if x < 0)} / 零 {sum(1 for x in w if x == 0)}")
print(f"  交叉核 ⟨w*,v_d⟩={mp.nstr(dot(w,V['d']),8)}  ⟨w*,v_s⟩={mp.nstr(dot(w,V['s']),8)}  (須 = ρ_d, 1)")

# ---------- E. 半整數位移恆等（/90 配 45−n；及 /88 配 44−n）----------
print("\n=== E. γ 配對恆等（慣例必須標明）===")
g90 = lambda n: 4*mp.sin(mp.pi*n/90)**2
g88 = lambda n: 4*mp.sin(mp.pi*n/88)**2
dev90 = max(abs(g90(n) + g90(45-n) - 4) for n in range(1, 45))
dev88 = max(abs(g88(n) + g88(44-n) - 4) for n in range(1, 44))
cross90 = max(abs(g90(n) + g90(44-n) - 4) for n in range(1, 45))
cross88 = max(abs(g88(n) + g88(45-n) - 4) for n in range(1, 44))
print(f"  /90 × 45−n : max|γ_n+γ_(45−n)−4| = {mp.nstr(dev90, 6)}  ✅")
print(f"  /88 × 44−n : max|γ_n+γ_(44−n)−4| = {mp.nstr(dev88, 6)}  ✅")
print(f"  交叉配對 /90 × 44−n : max dev = {mp.nstr(cross90, 8)}  ← 須等於 docs/40 之 0.139598")
print(f"  交叉配對 /88 × 45−n : max dev = {mp.nstr(cross88, 8)}  ← 須等於 docs/40 之 0.142678")
# φ 反射
dev_phi = max(abs(phi(45-n, j) + (-1)**j*phi(n, j)) for n in range(1, 45) for j in (2, 3, 4))
print(f"  φ_(45−n)(j) = −(−1)^j φ_n(j) : max dev = {mp.nstr(dev_phi, 6)}")

# ---------- F. 判決表 ----------
print("\n=== F. 重現閘判決 ===")
def dec(x): return mp.mpf(x)
rows = [
    ("R1 rho_d",      dec(mp.nstr(mp.mpf(rho_d.numerator)/rho_d.denominator, 12))),
    ("R2 rho_u",      dec(mp.nstr(mp.mpf(rho_u.numerator)/rho_u.denominator, 12))),
    ("R3 ratio",      dec(mp.nstr(mp.mpf(ratio.numerator)/ratio.denominator, 12))),
    ("R4 cos(vu,vd)", dot(V["u"], V["d"])/(norm(V["u"])*norm(V["d"]))),
    ("R5 angle",      ang),
    ("R6 sig/1_41",   sig_rd*41),
    ("R7 gamma pair", 4 - dev90),
    ("R7b cross dev", cross90),
]
status = {}
for name, got in rows:
    want, tol = EXPECT[name]
    d = abs(mp.mpf(got) - mp.mpf(want))
    ok = d <= mp.mpf(tol)
    st = "PASS" if ok else "FAIL"
    if name == "R6 sig/1_41" and not ok:
        # superseded 條款：對照 M12 R27 之更正值
        r27 = abs(sig_rd*41 - mp.mpf("1.8362"))
        if r27 < mp.mpf("5e-4") and abs(mp.mpf(want) - mp.mpf("1.4391")) < mp.mpf("5e-5"):
            st = "SUPERSEDED"
    status[name] = st
    print(f"  {name:15s} 得 {mp.nstr(got,10):>14s}  報 {want:>12s}  |Δ| = {mp.nstr(d,4):>10s}  {st}")

print(f"  γ 恆等 /90×45−n 精確成立（max dev {mp.nstr(dev90,3)}）  ⟹ R7 以恆等式記 PASS")
n_pass = sum(1 for v in status.values() if v == "PASS")
print(f"\n  合計: PASS {n_pass} / {len(status)}   SUPERSEDED {sum(1 for v in status.values() if v=='SUPERSEDED')}   "
      f"FAIL {sum(1 for v in status.values() if v=='FAIL')}")
print("  ⚠ R6 說明：docs/38 之 σ_ρd=0.03510 未計 δ_d–δ_s 協方差；M12 R27（docs/42）已更正為")
print("     0.04478518 ⟹ 比值 1.8362。本腳本獨立重算得 %s ⟹ 以 R27 為正典，docs/38 標 SUPERSEDED。"
      % mp.nstr(sig_rd, 8))
print("  宣告：N_trials = 0（純重現）。本腳本不作任何預言、不掃描、不改任何上游檔。")