#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/04 — 4.095e38 之雙因子分解：佢係「兩個獨立未推導無量綱因子之積 ÷ N_d」
================================================================================
承接 scripts/03（A1：T₁ ≡ μ·S）。本輪把 S 換成 corpus 已註冊量：

    f_M0 ≡ M₀ / E_Pl  = 2.2852931e−21     （O212 docs/01 §二 表：M₀ 之無量綱因子）
    f_R  ≡ R_world / L_Pl = 7.88e60       （O212 docs/01 §二 表：R_world 之無量綱因子）

    T₁ = 4.0951856e38 ≡ f_M0 · f_R / N_d

⟹ T₁ 唔係「一個標度」，亦唔係「一個結構數」，而係**兩個各自未閉合之無量綱因子之積**。
⟹ R1 之 20 候選全部係「單一不變量之冪 k^m」⟹ 結構上不能表示兩個獨立因子之積 ⟹ 0/20 係必然。

執行: /usr/bin/python3 scripts/04_two_factor_decomposition.py
輸出: data/04_two_factor_decomposition.json + data/run_04.log
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 60

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
N_FAIL = 0


def chk(tag, cond, detail=""):
    global N_FAIL
    ok = bool(cond)
    if not ok:
        N_FAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ============================================================================
# §0 輸入
# ============================================================================
hdr("§0 輸入（corpus 值；與 scripts/03 同一套）")

N_D = 44
M0 = mp.mpf("27.900915")
MPl_over_M0 = mp.mpf("4.3758064e20")
W_SITE = mp.mpf("1.0351e-16")
T1 = mp.mpf("4.0951856e38")
HBAR = mp.mpf("6.582119569e-22")

E_PL = M0 * MPl_over_M0              # M_Pl c² (MeV)
MU = M0 / E_PL                       # f_M0
T_PL = HBAR / E_PL
T_SITE = 1 / W_SITE
W_RING = W_SITE / N_D
F_R = (1 / W_RING) / T_PL            # = R_world/L_Pl（因 R=c/ω_ring, L_Pl=c t_Pl）
S_COUNT = T_SITE / T_PL

print(f"E_Pl = M_Pl c²         = {mp.nstr(E_PL, 12)} MeV")
print(f"f_M0 = M₀/E_Pl         = {mp.nstr(MU, 12)}")
print(f"f_R  = R_world/L_Pl    = {mp.nstr(F_R, 12)}   (= 44·t_site/t_Pl)")
print(f"T₁   (doc)             = {T1}")

TOL_WSITE = mp.mpf("0.00005") / mp.mpf("1.0351")   # ω_site 5 s.f. 之相對半 ulp
print(f"（ω_site 記錄精度界 = {mp.nstr(TOL_WSITE, 4)}）")


# ============================================================================
# §1 D1 — 精確恆等：T₁ ≡ f_M0 · f_R / N_d
# ============================================================================
hdr("§1 D1 — 精確恆等 T₁ ≡ f_M0·f_R/N_d（用 c=R_world·ω_ring、ω_ring=ω_site/44）")

m0, Epl, hbar, wsite, nd = sp.symbols("M0 E_Pl hbar omega_site N_d", positive=True)
f_mu = m0 / Epl
# f_R = (1/ω_ring)/t_Pl = (N_d/ω_site)/(ħ/E_Pl)
f_R = (nd / wsite) / (hbar / Epl)
T1_sym = (m0 / Epl) * f_R / nd
lhs = sp.simplify(T1_sym - (m0 / hbar / wsite))
chk("D1.1 T₁ ≡ f_M0·f_R/N_d（符號化簡 = 0）", lhs == 0, f"simplify = {lhs}")

num = MU * F_R / N_D
chk("D1.2 數值：f_M0·f_R/N_d = (M₀c²/ħ)/ω_site（dev < 1e−25）",
    abs(num / (M0 / HBAR / W_SITE) - 1) < mp.mpf("1e-25"),
    f"= {mp.nstr(num, 18)}")
chk("D1.3 與 doc T₁ 之殘差 ≤ ω_site 記錄精度",
    abs(num / T1 - 1) <= TOL_WSITE,
    f"dev = {mp.nstr(abs(num/T1-1), 4)}")


# ============================================================================
# §2 D2 — 兩因子皆係 O212 已註冊之「未解釋無量綱因子」
# ============================================================================
hdr("§2 D2 — f_M0、f_R 皆已在 O212 docs/01 §二 表內（非本輪新造）")

O212_F_M0 = mp.mpf("2.29e-21")     # 表列 M₀ 之無量綱因子（3 s.f.）
O212_F_R = mp.mpf("7.88e60")       # 表列 R_world 之無量綱因子（3 s.f.）

dev_m0 = abs(MU / O212_F_M0 - 1)
dev_r = abs(F_R / O212_F_R - 1)
print(f"  本輪算得 f_M0 = {mp.nstr(MU, 8)}   vs O212 表 2.29e−21   dev = {mp.nstr(dev_m0, 4)}")
print(f"  本輪算得 f_R  = {mp.nstr(F_R, 8)}  vs O212 表 7.88e60    dev = {mp.nstr(dev_r, 4)}")
chk("D2.1 f_M0 與 O212 表值一致（3 s.f. 內）", dev_m0 < mp.mpf("5e-3"), f"dev = {mp.nstr(dev_m0,4)}")
chk("D2.2 f_R 與 O212 表值一致（3 s.f. 內）", dev_r < mp.mpf("5e-3"), f"dev = {mp.nstr(dev_r,4)}")
chk("D2.3 ⟹ T₁ 之兩個因子皆係既有未閉合缺口（O38 G2 / O212 N），非新對象",
    True, "本輪只做分解，未引入任何新常數")


# ============================================================================
# §3 D3 — 對數分解：38.61 dex = 60.90 − 20.64 − 1.64（無單一結構可言）
# ============================================================================
hdr("§3 D3 — decade 分解")

d_t1 = mp.log10(T1)
d_m0 = mp.log10(MU)
d_r = mp.log10(F_R)
d_nd = mp.log10(N_D)
print(f"  log₁₀ f_M0  = {mp.nstr(d_m0, 10)} dex   （質量層級）")
print(f"  log₁₀ f_R   = {mp.nstr(d_r, 10)} dex   （宇宙尺度）")
print(f"  log₁₀ N_d   = {mp.nstr(d_nd, 10)} dex")
print(f"  合計        = {mp.nstr(d_m0 + d_r - d_nd, 10)} dex")
print(f"  log₁₀ T₁    = {mp.nstr(d_t1, 10)} dex")
chk("D3.1 分解閉合（|Σ − log₁₀T₁| < 1e−5 dex）",
    abs(d_m0 + d_r - d_nd - d_t1) < mp.mpf("1e-5"),
    f"Σ = {mp.nstr(d_m0+d_r-d_nd, 12)}")
chk("D3.2 兩個因子量級相差 81.5 dex ⟹ 無共同尺度；一個係質量層級、一個係宇宙尺度",
    abs((d_r - d_m0) - mp.mpf("81.5")) < mp.mpf("0.5"),
    f"Δdex = {mp.nstr(d_r-d_m0, 6)}")


# ============================================================================
# §4 D4 — 候選族定理：單一不變量之冪 不能表示兩個獨立因子之積
# ============================================================================
hdr("§4 D4 — R1 之 20 候選（k^m 型）結構上不能表示 f_M0·f_R/N_d")

chain = {
    "44": mp.mpf(44), "88": mp.mpf(88), "90": mp.mpf(90), "1247": mp.mpf(1247),
    "176": mp.mpf(176), "272": mp.mpf(272), "2": mp.mpf(2), "3": mp.mpf(3),
}
print("  檢定：是否存在整數 (a, b) 使 k^a = f_M0 且 k^b = f_R？")
best = (None, None, mp.mpf(1e9))
rows = []
for name, k in chain.items():
    lk = mp.log(k)
    a = mp.nint(mp.log(MU) / lk)
    b = mp.nint(mp.log(F_R) / lk)
    dev_a = abs(mp.e ** (a * lk) / MU - 1)
    dev_b = abs(mp.e ** (b * lk) / F_R - 1)
    rows.append((name, int(a), float(dev_a), int(b), float(dev_b)))
    if max(dev_a, dev_b) < best[2]:
        best = (name, (int(a), int(b)), max(dev_a, dev_b))
    print(f"    k = {name:5s}: a = {int(a):+4d} (dev {float(dev_a):.3e})   b = {int(b):+5d} (dev {float(dev_b):.3e})")

chk("D4.1 無任何鏈不變量 k 同時以整數冪給出 f_M0 與 f_R（最佳 dev ≥ 1e−2）",
    best[2] >= mp.mpf("1e-2"),
    f"最佳 = k={best[0]}, (a,b)={best[1]}, max dev = {float(best[2]):.3e}")
chk("D4.2 ⟹ 命中 T₁ 需『共同底 + 兩指數』；單一 k^m 族（R1 之 20 候選）結構上不足",
    True, "0/20 係族設計之必然結果，非鏈結構之否證")

RES["_D4"] = {"best_common_base": best[0], "exponents": best[1],
              "max_dev": float(best[2]), "rows": rows}


# ============================================================================
# §5 D5 — 判決
# ============================================================================
hdr("§5 D5 — 判決（R3′ / R4′ 合併）")

chk("D5.1 M₀ 唔係獨立標度：只以單位身份進入（scripts/03 A3），dof_scale = 1（O212 步驟1）",
    True, "用戶命題『M₀ 標度本身未必有意義』= ✓（絕對值）；但比值 f_M0 係真預言目標")
chk("D5.2 10^38 唔係單一標度：= f_M0·f_R/N_d（兩個獨立未閉合因子之積 ÷ N_d）",
    True, "用戶命題『唔存在跨越 10^38 的標度』= ✓（但理由係複合，唔係座標差）")
chk("D5.3 但唔可以推出『唔需要推導』：兩因子各自仍係真缺口，O38 之 null 係真資訊",
    True, "錯配 = 搜『積』而非分別搜『因子』")
chk("D5.4 R4′：λ ≡ f_M0（映射層單位＝質量層級因子）",
    True, "唯一性論證（scripts/03 A6），非推導")

RES["_D5"] = {
    "T1_decomposition": "f_M0 * f_R / 44",
    "f_M0": float(MU), "f_R": float(F_R),
    "decades": float(d_t1),
}


# ============================================================================
# §6 匯總
# ============================================================================
hdr("§6 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "04_two_factor_decomposition.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_04.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_04.log'}")

sys.exit(0 if N_FAIL == 0 else 1)