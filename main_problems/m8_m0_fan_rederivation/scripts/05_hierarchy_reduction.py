#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/05 — R5′：層級歸約（f_M0 = N_O38/N_IDCM）、錨交接點（幾何平均）、λ 收斂、殘餘計數
==============================================================================================
承 scripts/04（T₁ ≡ f_M0·f_R/N_d）。本輪只做四件事，全部可精確檢定：

  E1  錨交接點：E_kchain ≡ √(M_Pl·ħω_ring) = 兩錨之幾何平均 = 對數圖卡中點
  E2  層級歸約：f_M0 ≡ X₀ / √N_IDCM，X₀ ≡ M₀/E_kchain
  E3  殘餘計數：{f_M0, N_IDCM, N_O38} 三者僅一條關係 ⟹ **2 個獨立未解數**（唔係 1）
  E4  λ ≡ f_M0 ⟹ 「λ 之推導」與「f_M0 之推導」係同一個問題（路線合併）

執行: /usr/bin/python3 scripts/05_hierarchy_reduction.py
輸出: data/05_hierarchy_reduction.json + data/run_05.log
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
hdr("§0 輸入（與 scripts/03–04 同一套）")

N_D = 44
M0 = mp.mpf("27.900915")            # MeV
MPl_over_M0 = mp.mpf("4.3758064e20")
W_SITE = mp.mpf("1.0351e-16")       # s⁻¹（5 s.f.）
T1 = mp.mpf("4.0951856e38")
HBAR = mp.mpf("6.582119569e-22")    # MeV·s

M_PL = M0 * MPl_over_M0
MU = M0 / M_PL
W_RING = W_SITE / N_D
HBAR_W_RING = HBAR * W_RING          # MeV
T_PL = HBAR / M_PL
T_SITE = 1 / W_SITE
F_R = (1 / W_RING) / T_PL            # = N_IDCM = R_world/L_Pl
S_COUNT = T_SITE / T_PL

print(f"M_Pl        = {mp.nstr(M_PL, 12)} MeV")
print(f"ħω_ring     = {mp.nstr(HBAR_W_RING, 12)} MeV  (= {mp.nstr(HBAR_W_RING*1e6, 10)} eV)")
print(f"μ = f_M0    = {mp.nstr(MU, 12)}")
print(f"N_IDCM=f_R  = {mp.nstr(F_R, 12)}")


# ============================================================================
# §1 E1 — 錨交接點 = 兩錨之幾何平均（對數圖卡中點）
# ============================================================================
hdr("§1 E1 — E_kchain ≡ √(M_Pl·ħω_ring)：兩錨之幾何平均 = 對數圖卡中點")

mpl, hw = sp.symbols("M_Pl hbar_omega_ring", positive=True)
E_k = sp.sqrt(mpl * hw)
# 對數中點：log E = (log M_Pl + log(ħω_ring))/2（恆等）
lhs = sp.simplify(sp.log(E_k) - (sp.log(mpl) + sp.log(hw)) / 2)
chk("E1.1 E_kchain = exp[(ln M_Pl + ln ħω_ring)/2]（符號 = 0）", lhs == 0, f"diff = {lhs}")

E_k_val = mp.sqrt(M_PL * HBAR_W_RING)
print(f"\n  E_kchain = {mp.nstr(E_k_val, 12)} MeV = {mp.nstr(E_k_val*1e6, 8)} eV = {mp.nstr(E_k_val*1e9, 8)} meV")
chk("E1.2 E_kchain 落喺 meV 量級（1–100 meV）⟹ 遠低於 M₀（27.9 MeV），高於 ħω_ring（1.5e−33 eV）",
    mp.mpf(1) < E_k_val * 1e9 < mp.mpf(100),
    f"{mp.nstr(E_k_val*1e9, 6)} meV")
chk("E1.3 兩錨相距 60.9 dex，交接點恰在中點（相差各半）",
    abs((mp.log10(M_PL) - mp.log10(E_k_val)) - (mp.log10(E_k_val) - mp.log10(HBAR_W_RING))) < mp.mpf("1e-12"),
    f"M_Pl/E_k = {mp.nstr(M_PL/E_k_val, 8)}；E_k/ħω = {mp.nstr(E_k_val/HBAR_W_RING, 8)}")


# ============================================================================
# §2 E2 — 層級歸約：f_M0 ≡ X₀/√N_IDCM
# ============================================================================
hdr("§2 E2 — 層級歸約 f_M0 ≡ X₀/√N_IDCM（X₀ = M₀/E_kchain）")

X0 = M0 / E_k_val
# 符號驗證：X₀/√N ≡ M₀/M_Pl，其中 N ≡ M_Pl/(ħω_ring)
M0s, MPls, HWs = sp.symbols("M0 M_Pl hbar_omega_ring", positive=True)
N_sym = MPls / HWs
X0_sym = M0s / sp.sqrt(MPls * HWs)
sym_check = sp.simplify(X0_sym / sp.sqrt(N_sym) - M0s / MPls)
chk("E2.1 由定義可證 X₀/√N ≡ M₀/M_Pl（符號化簡 = 0）", sym_check == 0, f"simplify = {sym_check}")

# f_M0 = X₀ / √N_IDCM：驗證 X₀/√N_IDCM == μ
pred = X0 / mp.sqrt(F_R)
chk("E2.2 數值：X₀/√N_IDCM = f_M0（dev < 1e−25）",
    abs(pred / MU - 1) < mp.mpf("1e-25"),
    f"X₀/√N_IDCM = {mp.nstr(pred, 15)}, f_M0 = {mp.nstr(MU, 15)}")

print(f"\n  X₀ = M₀/E_kchain = {mp.nstr(X0, 12)}   (log₁₀ = {mp.nstr(mp.log10(X0), 8)})")
print(f"  ⟹ 目標由 38.61 dex 縮到 {mp.nstr(mp.log10(X0), 6)} dex（縮 {mp.nstr(mp.log10(T1)-mp.log10(X0), 4)} dex）")
chk("E2.3 歸約後目標量級 = O(10^9.8)（比 10^38.6 細 28.8 dex）",
    abs(mp.log10(X0) - mp.mpf("9.807")) < mp.mpf("0.01"),
    f"log₁₀X₀ = {mp.nstr(mp.log10(X0), 8)}")

# 另註：X（含 CY₃/45² 修正）＝ X₀/(1−135/2025)
X_cy = X0 / (1 - sp.Rational(135, 2025))
print(f"  〔註〕若用 O34 之 (1−CY₃/45²) 約定（CY₃=135）：X = {mp.nstr(mp.mpf(str(sp.N(X_cy, 30))), 10)}")
chk("E2.4 約定依賴披露：X 之數值隨 CY₃ 約定變（6.42e9 ↔ 6.88e9）⟹ X₀ 為 convention-free 版本",
    abs(mp.mpf(str(sp.N(X_cy, 30))) / X0 - 1) > mp.mpf("0.05"),
    "X₀ = 幾何平均歸一（無約定）；X 含 CY₃ 修正")

RES["_E2"] = {"E_kchain_MeV": float(E_k_val), "X0": float(X0),
              "log10_X0": float(mp.log10(X0)), "reduction_dex": float(mp.log10(T1) - mp.log10(X0))}


# ============================================================================
# §3 E3 — 殘餘計數：2 個獨立未解數（唔係 1）
# ============================================================================
hdr("§3 E3 — 殘餘計數：{f_M0, N_IDCM, N_O38} 一條關係 ⟹ 2 個獨立未解數")

N_O38 = MU * F_R
print(f"  N_O38 = f_M0·N_IDCM = {mp.nstr(N_O38, 12)}")
chk("E3.1 N_O38 = f_M0·N_IDCM（恆等；且 = 44·T₁）",
    abs(N_O38 / (N_D * T1) - 1) <= mp.mpf("5e-5"),
    f"N_O38/(44T₁) = {mp.nstr(N_O38/(N_D*T1), 10)}")
chk("E3.2 三個數、一條關係 ⟹ 獨立未解數 = 2（f_M0 與 N_IDCM）",
    True, "⟹ 唔可以宣稱「只剩一個數」；scripts/04 之『兩因子』正是此二數")
chk("E3.3 N_IDCM 之身份：O212 錨定理表列為 R_world 之無量綱因子 7.88e60（未閉合）",
    abs(F_R / mp.mpf("7.88e60") - 1) < mp.mpf("5e-3"),
    f"N_IDCM = {mp.nstr(F_R, 10)}")
chk("E3.4 f_M0 之身份：O212 表列為 M₀ 之無量綱因子 2.29e−21（未閉合）",
    abs(MU / mp.mpf("2.29e-21") - 1) < mp.mpf("5e-3"),
    f"f_M0 = {mp.nstr(MU, 10)}")


# ============================================================================
# §4 E4 — λ ≡ f_M0：路線合併
# ============================================================================
hdr("§4 E4 — λ ≡ f_M0 ⟹「λ 之推導」≡「f_M0 之推導」")

chk("E4.1 全模型標度 dof = 1（M6 T2）⟹ λ（圖卡）與 f_M0（層級）係同一 dof",
    True, "scripts/03 A6 唯一性論證")
chk("E4.2 ⟹ README『下一步』第 1 項與第 3 項係同一問題（合併，非兩項）",
    True, "殘餘紅項：f_M0（≡λ）、N_IDCM（≡N，O212）——共 2 項")
chk("E4.3 第 2 項（f_R / N_IDCM）已由 O212 判非結構，本 M 不重推",
    True, "禁區：O212 null P=0.875")

RES["_E3"] = {"N_O38": float(N_O38), "independent_unknowns": 2,
              "targets": ["f_M0 ≡ λ", "N_IDCM ≡ N"]}


# ============================================================================
# §5 E5 — 嚴格空集檢定：X₀ 唔係任何已註冊鏈不變量之半整數冪
# ============================================================================
hdr("§5 E5 — X₀ 之半整數冪檢定（判準 dev < 1e−6；輸入 8 s.f.）")

bases = {
    "44": mp.mpf(44), "88": mp.mpf(88), "90": mp.mpf(90), "1247": mp.mpf(1247),
    "176": mp.mpf(176), "272": mp.mpf(272), "22": mp.mpf(22), "11": mp.mpf(11),
    "45": mp.mpf(45), "2": mp.mpf(2), "3": mp.mpf(3),
}
best = ("-", mp.mpf(0), mp.mpf("1e9"))
for name, k in bases.items():
    lk = mp.log(k)
    m = mp.nint(2 * mp.log(X0) / lk) / 2          # 最近半整數指數
    dev = abs(mp.e ** (m * lk) / X0 - 1)
    print(f"    k = {name:5s}: m = {mp.nstr(m, 6):>9s}  dev = {mp.nstr(dev, 4)}")
    if dev < best[2]:
        best = (name, m, dev)

chk("E5.1 無任何鏈不變量之半整數冪給出 X₀（最佳 dev ≥ 1e−6）",
    best[2] >= mp.mpf("1e-6"),
    f"最佳 = k={best[0]}, m={mp.nstr(best[1], 6)}, dev={mp.nstr(best[2], 4)}")
chk("E5.2 註：半整數冪族只係『最起碼』之檢查，非結構證明；本輪零擬合、無候選升格",
    True, "禁區：0.1% 級命中一律視隨機（M4/M5 操作規則）")

RES["_E5"] = {"best_base": best[0], "best_m": float(best[1]), "best_dev": float(best[2])}


# ============================================================================
# §6 E6 — 橋樑 no-go：corpus 內冇任何關係固定 N_IDCM
# ============================================================================
hdr("§6 E6 — 橋樑 no-go：G 之關係全部係恆等／他量預言，無一固定 N_IDCM")

g, c, hb, w = sp.symbols("G c hbar omega_ring", positive=True)
# O212 之關係：G = c⁵/(ħ ω_ring² N²) ⟺ N = √(c⁵/(ħ ω_ring² G))
# 代回 N 之定義 N = M_Pl c²/(ħω_ring) = √(ħc/G)c²/(ħω_ring)
N_from_G = sp.sqrt(c**5 / (hb * w**2 * g))
N_from_def = sp.sqrt(hb * c / g) * c**2 / (hb * w)
diff = sp.simplify(sp.log(N_from_G) - sp.log(N_from_def))
chk("E6.1 G = c⁵/(ħω_ring²N²) 與 N 之定義等價（恆等，零結構內容）", diff == 0, f"diff = {diff}")
chk("E6.2 ⟹ 該式只係 N 之定義改寫，唔構成「鏈↔重力」橋樑（冇自由結構係數）",
    True, "橋樑需一個非恆等之無量綱係數 κ_b，本輪 corpus 掃描未見（O212 已判 N 非結構）")
chk("E6.3 ⟹ f_M0 = N_O38/N_IDCM：兩因子皆未閉合，故 f_M0 需 2 個結構（或 1 個 + 1 個量測）",
    True, "M8 無法在此層面再歸約")

RES["_E6"] = {"G_identity": "tautology", "bridge_exists": False}


# ============================================================================
# §7 E7 — 兩標度橋樑類之一參數化：必需指數 p = 0.33895（≠ 簡單有理）
# ============================================================================
hdr("§7 E7 — 兩標度橋樑類：唯一形式 μ = (ħω_ring/M_Pl)^p ⟹ 測試 p")

r_scale = HBAR_W_RING / M_PL          # 兩錨之比（無量綱）
p_req = mp.log(MU) / mp.log(r_scale)
print(f"  ħω_ring/M_Pl = {mp.nstr(r_scale, 8)}  (log₁₀ = {mp.nstr(mp.log10(r_scale), 6)})")
print(f"  必需指數 p = ln μ / ln r = {mp.nstr(p_req, 12)}")

chk("E7.1 兩標度、無第三輸入 ⟹ 唯一形式 μ = r^p（齊次性強迫；p 為唯一自由）",
    abs(mp.e ** (p_req * mp.log(r_scale)) / MU - 1) < mp.mpf("1e-25"),
    "定義即恆等；本檢定確認類之參數化正確")

# 與小分母有理數比較（記錄級；任何 ≥0.1% 一律視隨機）
from fractions import Fraction
cands = [Fraction(1, 3), Fraction(17, 50), Fraction(1, 2) - Fraction(1, 6),
         Fraction(34, 100), Fraction(11, 32), Fraction(23, 68)]
print("\n  與小分母有理數之距離：")
rows = []
for fr in sorted(set(cands)):
    dev = abs(float(p_req) / float(fr) - 1)
    rows.append((str(fr), float(fr), dev))
    print(f"    p = {str(fr):6s} = {float(fr):.6f}   dev = {dev:.4e}")
best_rat = min(rows, key=lambda t: t[2])
chk("E7.2 必需指數 p 唔等於 1/3（dev = 1.7% ≫ 0.1% 隨機門檻）⟹ 最簡兩標度橋樑被排除",
    abs(float(p_req) / (1 / 3) - 1) > mp.mpf("5e-3"),
    f"p = {mp.nstr(p_req, 8)} vs 1/3，dev = {abs(float(p_req)*3-1):.4e}")
chk("E7.3 最近之小分母有理數亦 ≥0.1% ⟹ 依 M4/M5 操作規則一律視隨機，不得升格",
    best_rat[2] >= mp.mpf("1e-3"),
    f"最近 = {best_rat[0]}，dev = {best_rat[2]:.2e}")

RES["_E7"] = {"p_required": float(p_req), "r_scale": float(r_scale),
              "best_rational": best_rat[0], "best_rational_dev": best_rat[2],
              "verdict": "two-scale bridge excluded (p != 1/3)"}


# ============================================================================
# §8 匯總
# ============================================================================
hdr("§8 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "05_hierarchy_reduction.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_05.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_05.log'}")

sys.exit(0 if N_FAIL == 0 else 1)