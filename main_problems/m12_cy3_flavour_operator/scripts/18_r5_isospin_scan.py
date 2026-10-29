#!/usr/bin/env python3
"""R5 — 味對荷之全通道掃描（加大範圍；本輪不裁決）

資料源（vendored, 唯讀）: data/pdg2024/mass_width_2024.txt
  PDG 2024 RPP masses/widths/MC-IDs, sha256 9d6024da6705fcab556f1136053c802f…
  下載日 2026-10-29；單位 GeV，本腳本換算為 MeV。
  夸克內容以 PDG MC-ID 標準對應硬編碼（見 CONTENT 表）——內容非擬合參數。

預註冊判準（2026-10-29 跑前凍結，不得事後改）
--------------------------------------------------
S1 一體荷判別（每夸克荷）：Δm = Σ_i a(f_i) ⟹ 光夸克四通道全等
       n−p = Σ⁻−Σ⁰ = Σ⁰−Σ⁺ = Ξ⁻−Ξ⁰
   逐通道 |z| ≥ 3 判死。

S2 二體荷（每價對荷）：Δm = Σ_(價對) a(f_i,f_j)。定義
       δ_x ≡ a(d,x) − a(u,x)   (x = u,d,s,c,b)
   則光夸克八重態
       n−p   = δ_u + δ_d
       Σ⁻−Σ⁰ = δ_d + δ_s
       Ξ⁻−Ξ⁰ = 2 δ_s
   ⟹ 3 未知 3 方程，**δ 唯一解**（無自由度）。其餘全部為樣本外預言。

S3 樣本外預言清單（跑前寫死，共 8 條，禁增刪）
   [A1] Σ⁰−Σ⁺          = δ_u + δ_s
   [A2] Σ*⁻−Σ*⁰        = δ_d + δ_s        （自旋盲：decuplet 同價對）
   [A3] Σ*⁰−Σ*⁺        = δ_u + δ_s        （自旋盲）
   [A4] Ξ*⁻−Ξ*⁰        = 2 δ_s            （自旋盲）
   [B1] (Σc⁰−Σc⁺)−(Σc⁺−Σc⁺⁺) = δ_d − δ_u
   [B2] (Σc⁰−Σc⁺)+(Ξc⁰−Ξc⁺)  = δ_d − δ_s
   [C1] (Σb⁻−Σb⁺) − 2(Ξb⁻−Ξb⁰) = (n−p) − (Ξ⁻−Ξ⁰)
   [C2] δ_c := (Σc⁰−Σc⁺) − δ_d ；預言 Ξc⁰−Ξc⁺ = −(δ_s + δ_c)

S4 判死門檻 |z| ≥ 3；本輪**只報通過/總數，不作裁決**（依指示：未出結論前不裁定）。
S5 紀律：禁由任一被測通道回擬 δ；禁新增參數；禁剔除通道；禁改質量表。
S6 介子扇區（K, K*, D, D*, B）：每型只有一個價對 ⟹ 對配對 ansatz 零約束，
   只報 δ^q q̄ 值與其與 δ^qq 之比（正反不對稱），不作判據。
"""
import re
import pathlib
import mpmath as mp

mp.mp.dps = 40

SRC = pathlib.Path("data/pdg2024/mass_width_2024.txt")

# PDG MC-ID → 夸克內容（標準 MC 編號約定；非擬合參數）
CONTENT = {
    "2112": "udd",   # n
    "2212": "uud",   # p
    "3112": "dds",   # Sigma-
    "3212": "uds",   # Sigma0
    "3222": "uus",   # Sigma+
    "3312": "dss",   # Xi-
    "3322": "uss",   # Xi0
    "3114": "dds",   # Sigma(1385)-
    "3214": "uds",   # Sigma(1385)0
    "3224": "uus",   # Sigma(1385)+
    "3314": "dss",   # Xi(1530)-
    "3324": "uss",   # Xi(1530)0
    "3334": "sss",   # Omega-
    "4112": "ddc",   # Sigma_c(2455)0
    "4212": "udc",   # Sigma_c(2455)+
    "4222": "uuc",   # Sigma_c(2455)++
    "4132": "usc",   # Xi_c0
    "4232": "dsc",   # Xi_c+
    "5112": "bdd",   # Sigma_b-
    "5222": "buu",   # Sigma_b+
    "5132": "bds",   # Xi_b-
    "5232": "bus",   # Xi_b0
    "321":  "us~",   # K+
    "311":  "ds~",   # K0
    "323":  "us~",   # K*(892)+
    "313":  "ds~",   # K*(892)0
    "411":  "cd~",   # D+
    "421":  "cu~",   # D0
    "413":  "cd~",   # D*(2010)+
    "423":  "cu~",   # D*(2007)0
    "511":  "db~",   # B0
    "521":  "ub~",   # B+
}

LABEL = {
    "2112": "n", "2212": "p", "3112": "Σ⁻", "3212": "Σ⁰", "3222": "Σ⁺",
    "3312": "Ξ⁻", "3322": "Ξ⁰", "3114": "Σ*⁻", "3214": "Σ*⁰", "3224": "Σ*⁺",
    "3314": "Ξ*⁻", "3324": "Ξ*⁰", "3334": "Ω⁻",
    "4112": "Σc⁰", "4212": "Σc⁺", "4222": "Σc⁺⁺", "4132": "Ξc⁰", "4232": "Ξc⁺",
    "5112": "Σb⁻", "5222": "Σb⁺", "5132": "Ξb⁻", "5232": "Ξb⁰",
    "321": "K⁺", "311": "K⁰", "323": "K*⁺", "313": "K*⁰",
    "411": "D⁺", "421": "D⁰", "413": "D*⁺", "423": "D*⁰",
    "511": "B⁰", "521": "B⁺",
}


def load():
    """解析 PDG 檔 → {id: (m_MeV, sigma_MeV)}"""
    out = {}
    for line in SRC.read_text().splitlines():
        if line.startswith("*"):
            continue
        m = re.match(r"^\s*(-?\d+)\s+(.+)$", line)
        if not m:
            continue
        pid, rest = m.group(1), m.group(2)
        nums = re.findall(r"[+-]?\d\.\d+E[+-]\d+", rest)
        if len(nums) < 2:
            continue
        out[pid] = (mp.mpf(nums[0]) * 1000, abs(mp.mpf(nums[1])) * 1000)
    return out


MS = load()
miss = [p for p in CONTENT if p not in MS]
assert not miss, "缺 PDG ID: %s" % miss


def m(pid):
    return MS[pid][0]


def s(pid):
    return MS[pid][1]


def dmm(a, b):
    """Δm(a−b) 與 σ（獨立誤差）"""
    return m(a) - m(b), mp.sqrt(s(a) ** 2 + s(b) ** 2)


def fmt(x, w=10, d=4):
    return ("%" + str(w) + "." + str(d) + "f") % x


print("=" * 92)
print("R5 — 味對荷全通道掃描（PDG 2024 vendored；本輪只報不裁決）")
print("=" * 92)
print("  資料源 sha256(前32) =", __import__("hashlib").sha256(SRC.read_bytes()).hexdigest()[:32])
print("  通道數 = %d（夸克內容硬編碼，非擬合參數）" % len(CONTENT))

# ------------------------------------------------ §1 原始分割表
print("\n[1] 同味多重態 isospin 分割（Δm = m_高 − m_低 按表定方向）")
CH = [("n", "p"), ("Σ⁻", "Σ⁰"), ("Σ⁰", "Σ⁺"), ("Ξ⁻", "Ξ⁰"),
      ("Σ*⁻", "Σ*⁰"), ("Σ*⁰", "Σ*⁺"), ("Ξ*⁻", "Ξ*⁰"),
      ("Σc⁰", "Σc⁺"), ("Σc⁺", "Σc⁺⁺"), ("Ξc⁰", "Ξc⁺"),
      ("Σb⁻", "Σb⁺"), ("Ξb⁻", "Ξb⁰"),
      ("K⁺", "K⁰"), ("K*⁺", "K*⁰"), ("D⁺", "D⁰"), ("D*⁺", "D*⁰"), ("B⁰", "B⁺")]
INV = {}
for pid, lb in LABEL.items():
    INV[lb] = pid
RAW = {}
for a, b in CH:
    d, sg = dmm(INV[a], INV[b])
    RAW[(a, b)] = (d, sg)
    print("   %-10s − %-6s = %+10.4f ± %7.4f MeV" % (a, b, d, sg))

# ------------------------------------------------ §2 一體 vs 二體
print("\n" + "=" * 92)
print("[S1] 一體荷：光夸克四通道全等 ⟹ 逐通道判死")
anchors = [("n", "p"), ("Σ⁻", "Σ⁰"), ("Σ⁰", "Σ⁺"), ("Ξ⁻", "Ξ⁰")]
ref, ref_s = RAW[("n", "p")]
s1_kill = 0
for a, b in anchors:
    d, sg = RAW[(a, b)]
    z = (d - ref) / mp.sqrt(sg ** 2 + ref_s ** 2)
    bad = abs(z) >= 3
    s1_kill += bad
    print("   %-16s 預言 %+9.4f   實測 %+9.4f   z = %+8.2f   %s"
          % ("%s−%s" % (a, b), ref, d, z, "判死 🔴" if bad else "通 🟢"))
print("   ⟹ 一體荷：%d/%d 通道判死" % (s1_kill, len(anchors)))

print("\n[S2] 二體荷：由三錨唯一解 δ（無自由度）")
d_np, s_np = RAW[("n", "p")]
d_SS, s_SS = RAW[("Σ⁻", "Σ⁰")]
d_XX, s_XX = RAW[("Ξ⁻", "Ξ⁰")]
d_s = d_XX / 2
d_d = d_SS - d_s
d_u = d_np - d_d
print("   δ_s = (Ξ⁻−Ξ⁰)/2      = %+10.4f MeV" % d_s)
print("   δ_d = (Σ⁻−Σ⁰) − δ_s  = %+10.4f MeV" % d_d)
print("   δ_u = (n−p) − δ_d    = %+10.4f MeV   ⟹ 輕夸克 δ 並非單調、且 δ_u < 0" % d_u)
print("   註：δ_x ≡ a(d,x) − a(u,x)；1 體荷之對應量係 δ_x 全等（即 S1 之預言）")

# ------------------------------------------------ §3 樣本外預言
print("\n" + "=" * 92)
print("[S3] 樣本外預言 8 條（判死門檻 |z| ≥ 3）")

TESTS = []
# A1
d, sg = RAW[("Σ⁰", "Σ⁺")]
TESTS.append(("[A1] Σ⁰−Σ⁺", d_u + d_s, d, sg))
# A2-A4 自旋盲
for nm, key, pred in [("[A2] Σ*⁻−Σ*⁰", ("Σ*⁻", "Σ*⁰"), d_d + d_s),
                      ("[A3] Σ*⁰−Σ*⁺", ("Σ*⁰", "Σ*⁺"), d_u + d_s),
                      ("[A4] Ξ*⁻−Ξ*⁰", ("Ξ*⁻", "Ξ*⁰"), 2 * d_s)]:
    d, sg = RAW[key]
    TESTS.append((nm, pred, d, sg))
# B1 / B2（粲：δ_c 消去）
dS0c, sS0c = RAW[("Σc⁰", "Σc⁺")]
dScp, sScp = RAW[("Σc⁺", "Σc⁺⁺")]
dXc, sXc = RAW[("Ξc⁰", "Ξc⁺")]
TESTS.append(("[B1] (Σc⁰−Σc⁺)−(Σc⁺−Σc⁺⁺)", d_d - d_u, dS0c - dScp,
              mp.sqrt(sS0c ** 2 + sScp ** 2)))
TESTS.append(("[B2] (Σc⁰−Σc⁺)+(Ξc⁰−Ξc⁺)", d_d - d_s, dS0c + dXc,
              mp.sqrt(sS0c ** 2 + sXc ** 2)))
# C1（底）
dSb, sSb = RAW[("Σb⁻", "Σb⁺")]
dXb, sXb = RAW[("Ξb⁻", "Ξb⁰")]
TESTS.append(("[C1] (Σb⁻−Σb⁺)−2(Ξb⁻−Ξb⁰)", d_np - d_XX, dSb - 2 * dXb,
              mp.sqrt(sSb ** 2 + (2 * sXb) ** 2)))
# C2（δ_c 由 Σc 解，再預言 Ξc）
d_c = dS0c - d_d
TESTS.append(("[C2] Ξc⁰−Ξc⁺（δ_c 自 Σc 解）", -(d_s + d_c), dXc, sXc))
print("\n   δ_c = (Σc⁰−Σc⁺) − δ_d = %+10.4f MeV" % d_c)

npass = nkill = 0
for nm, pred, obs, sg in TESTS:
    res = obs - pred
    z = res / sg
    bad = abs(z) >= 3
    nkill += bad
    npass += (not bad)
    print("   %-32s 預言 %+10.4f  實測 %+10.4f  殘差 %+9.4f ± %6.4f  z = %+7.2f  %s"
          % (nm, pred, obs, res, sg, z, "判死 🔴" if bad else "通 🟢"))

print("\n   掃描 tally：通 %d / 判死 %d / 總 %d（不裁決）" % (npass, nkill, len(TESTS)))
print("   ⚠ B2 與 C2 代數上係同一方程（δ_c 由 Σc⁰−Σc⁺ 定義 ⟹ 殘差恆同）")
print("     ⟹ 有效獨立檢定 = %d 條，判死 %d 條" % (len(TESTS) - 1, nkill - 1))

# ------------------------------------------------ §3b 通道相依 δ^eff
print("\n" + "=" * 92)
print("[S2b] 通道相依 δ^eff（把每個通道反解成本 ansatz 下之 δ，看是否為常數）")
print("   δ_s^eff(通道) ≡ 該通道實測值減去已解出之其餘 δ；若 ansatz 真，全部應等於上表 δ")
rows_eff = [
    ("Ξ⁻−Ξ⁰（錨, 1/2⁺, ss）", d_XX / 2, d_s),
    ("Σ⁰−Σ⁺（1/2⁺, s）", RAW[("Σ⁰", "Σ⁺")][0] - d_u, d_s),
    ("Σ*⁻−Σ*⁰（3/2⁺, s）", RAW[("Σ*⁻", "Σ*⁰")][0] - d_d, d_s),
    ("Σ*⁰−Σ*⁺（3/2⁺, s）", RAW[("Σ*⁰", "Σ*⁺")][0] - d_u, d_s),
    ("Ξ*⁻−Ξ*⁰（3/2⁺, ss）", RAW[("Ξ*⁻", "Ξ*⁰")][0] / 2, d_s),
]
for nm, eff, ref_s_val in rows_eff:
    print("   %-26s δ_s^eff = %+9.4f   （基準 δ_s = %+9.4f，差 %+9.4f MeV）"
          % (nm, eff, ref_s_val, eff - ref_s_val))
rows_eff_d = [
    ("Σ⁻−Σ⁰（錨, 1/2⁺, s）", d_SS - d_s),
    ("Σ*⁻−Σ*⁰（3/2⁺, s）", RAW[("Σ*⁻", "Σ*⁰")][0] - d_s),
]
print("   ── δ_d^eff ──")
for nm, eff in rows_eff_d:
    print("   %-26s δ_d^eff = %+9.4f   （基準 δ_d = %+9.4f，差 %+9.4f MeV）"
          % (nm, eff, d_d, eff - d_d))

# ------------------------------------------------ §4 介子扇區
print("\n" + "=" * 92)
print("[S6] 介子扇區：每型單一價對 ⟹ 對配對 ansatz 零約束；只報值")
mes = [("K⁺−K⁰", ("K⁺", "K⁰")), ("K*⁺−K*⁰", ("K*⁺", "K*⁰")),
       ("D⁺−D⁰", ("D⁺", "D⁰")), ("D*⁺−D*⁰", ("D*⁺", "D*⁰")),
       ("B⁰−B⁺", ("B⁰", "B⁺"))]
mesv = {}
for nm, key in mes:
    d, sg = RAW[key]
    mesv[nm] = d
    print("   %-12s = %+10.4f ± %7.4f MeV" % (nm, d, sg))
ratio = abs(mesv["K⁺−K⁰"]) / (abs(d_s))
print("   正反不對稱 ρ_s ≡ |δ^q q̄_s| / |δ^qq_s| = %.6f  （前輪 R4-P2 由 r 比得 1.148613）" % ratio)
print("   ρ_c ≡ |D⁺−D⁰| / |δ_c| = %.6f" % (abs(mesv["D⁺−D⁰"]) / abs(d_c) if d_c else mp.nan))

# ------------------------------------------------ §5 匯總（不裁決）
print("\n" + "=" * 92)
print("[匯總 — 只報事實，不作裁定]")
print("   光夸克 δ：δ_u = %+8.4f  δ_d = %+8.4f  δ_s = %+8.4f MeV" % (d_u, d_d, d_s))
print("   一體荷：%d/%d 通道 |z| ≥ 3" % (s1_kill, len(anchors)))
print("   二體荷：樣本外 8 條中 通 %d / 判死 %d" % (npass, nkill))
print("   介子扇區：5 值，零約束")
print("=" * 92)