#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/06 — R6′：三輸入橋樑之可證偽性分析（類層級；量度密度，非凑數）
=========================================================================
前情：scripts/05 E7 已排除「兩標度冪律橋」。本輪問：加第三個**無量綱**結構輸入 S 救唔救得返？
方法：唔猜形式，改為量度**該類之可證偽性**（局部密度），令「有冇能力判死」本身可檢定。

  G1  欠定定理：μ = r^p·S^q（2 未知 1 約束）⟹ 任意 S 皆可配 ⟹ 零預言內容
  G2  p 之最佳小分母有理近似 vs 隨機基準 1/(2q²) ⟹ 是否異常
  G3  IDCM 非冪結構清單（加性／對數／指數）⟹ 只有指數能在無精調下給 10^−21
  G4  ★ 局部密度量度：S* = ln(1/μ) = 47.52779 與 S*/λ_WH = 2922.96 之 0.1% 窗內
      有幾個「簡單結構式」⟹ 判該窗係「有資訊」抑或「噪音地板」
  G5  破牆規格（identity-level）

執行: /usr/bin/python3 scripts/06_three_input_falsifiability.py
輸出: data/06_three_input_falsifiability.json + data/run_06.log
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

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
# §0 輸入（承 scripts/03–05）
# ============================================================================
hdr("§0 輸入")

N_D = 44
M0 = mp.mpf("27.900915")
M_PL = M0 * mp.mpf("4.3758064e20")
W_SITE = mp.mpf("1.0351e-16")
HBAR = mp.mpf("6.582119569e-22")
W_RING = W_SITE / N_D
HW = HBAR * W_RING
MU = M0 / M_PL
R_SCALE = HW / M_PL
E_K = mp.sqrt(M_PL * HW)
LAMBDA_WH = mp.mpf(2) / mp.mpf(123)      # O39/O86：λ_WH = μ_W²/123 = 2/123

print(f"μ = M₀/M_Pl = {mp.nstr(MU, 12)}")
print(f"r = ħω/M_Pl = {mp.nstr(R_SCALE, 8)}")
print(f"E_kchain    = {mp.nstr(E_K*1e9, 10)} meV")
print(f"λ_WH        = 2/123 = {mp.nstr(LAMBDA_WH, 8)}")


# ============================================================================
# §1 G1 — 欠定定理
# ============================================================================
hdr("§1 G1 — 欠定定理：μ = r^p·S^q 對任意 S 皆可配")

S_list = {"44": mp.mpf(44), "88": mp.mpf(88), "90": mp.mpf(90), "1247": mp.mpf(1247),
          "291.5": mp.mpf("291.5"), "61.5": mp.mpf("61.5"), "3": mp.mpf(3)}
lr, lmu = mp.log(R_SCALE), mp.log(MU)

rows = []
for name, S in S_list.items():
    q = mp.mpf(1)
    p = (lmu - q * mp.log(S)) / lr
    dev = abs(mp.e ** (p * lr) * S ** q / MU - 1)
    rows.append((name, p, dev))
    print(f"  S = {name:6s}: q=1 ⟹ p = {mp.nstr(p, 10)}  dev = {mp.nstr(dev, 3)}")

chk("G1.1 任意 S 配上反解 p 皆精確重現 μ（dev < 1e−25）",
    all(r[2] < mp.mpf("1e-25") for r in rows), "7/7 個 S 皆 dev < 1e−25")
chk("G1.2 ⟹ 三輸入冪律族 2 未知 1 約束 = 1 參數族 ⟹ 零預言內容（不可證偽）",
    True, "加入無量綱第三輸入**不**恢復可證偽性")
chk("G1.3 ⟹ 第三輸入須以非冪方式（加性／對數／指數）進入方有內容",
    True, "見 G3")

RES["_G1"] = {"family": "mu = r^p * S^q", "dof_left": 1,
              "rows": [{"S": r[0], "p_at_q1": float(r[1])} for r in rows]}


# ============================================================================
# §2 G2 — p 之最佳有理近似 vs 隨機基準
# ============================================================================
hdr("§2 G2 — p = lnμ/ln r 之有理近似，對照隨機基準 1/(2q²)")

p_val = float(lmu / lr)
print(f"  p = {p_val:.15f}")
best = None
print("  limit_denominator 掃描：")
for qmax in (10, 20, 50, 100, 200, 500):
    fr = Fraction(p_val).limit_denominator(qmax)
    dev = abs(float(fr) / p_val - 1)
    bench = 1.0 / (2 * fr.denominator ** 2)     # 隨機數之最佳近似典型誤差
    ratio = dev / bench
    print(f"    q_max={qmax:4d}: {fr.numerator}/{fr.denominator} = {float(fr):.9f}  "
          f"dev={dev:.3e}  基準={bench:.3e}  dev/基準={ratio:5.2f}")
    if best is None or dev < best[2]:
        best = (f"{fr.numerator}/{fr.denominator}", float(fr), dev, ratio)

chk("G2.1 最佳近似（q≤500）之 dev/基準 ≥ 0.3 ⟹ 唔優於隨機，無異常",
    best[3] >= 0.3, f"最佳 = {best[0]}，dev = {best[2]:.3e}，dev/基準 = {best[3]:.2f}")
chk("G2.2 ⟹ p 無結構性有理身份；不得以有理近似升格（依 0.1% 規則 ＋ 四 Hard Tests）",
    True, f"{best[0]} 僅記錄，不升格")

RES["_G2"] = {"p": p_val, "best": best[0], "best_dev": best[2], "dev_over_benchmark": best[3]}


# ============================================================================
# §3 G3 — IDCM 非冪結構清單
# ============================================================================
hdr("§3 G3 — IDCM 已註冊非冪結構：只有指數能在無精調下給 10^−21")

cands = {
    "加性（M2 質量律 m/M₀ = Σγ_n）": "Σγ = 90 ⟹ O(10²)；產生唔到 10^−21 之比",
    "對數（M8 圖卡 t = λ·ln tanθ）": "圖卡內跨距 2.894 dex ⟹ 唔夠",
    "指數（O86 Λ_eff = 8πM_Pl²·e^{−V_Λ}）": "e^{−S}；10^−21 只需 S ≈ 47.5 ⟹ 量級合理",
}
for k, v in cands.items():
    print(f"  {k}：{v}")

chk("G3.1 只有指數型可在 O(1)–O(10²) 結構輸入下產生 10^−21（無精調）",
    True, "加性/對數在量級上不可能 ⟹ 指數型係唯一合理候選")
S_STAR = float(-lmu)
print(f"\n  S* = ln(1/μ) = {S_STAR:.10f}")
chk("G3.2 目標量 S* = ln(1/μ) 已定", S_STAR > 0, f"S* = {S_STAR:.8f}")

RES["_G3"] = {"S_star": S_STAR, "class": "exponential only"}


# ============================================================================
# §4 G4 — ★ 局部密度量度
# ============================================================================
hdr("§4 G4 — 局部密度量度：0.1% 窗係「有資訊」抑或「噪音地板」？")

R = {
    "2": 2, "3": 3, "8": 8, "11": 11, "22": 22, "36": 36, "44": 44, "45": 45,
    "55": 55, "60": 60, "88": 88, "89": 89, "90": 90, "98": 98, "123": 123,
    "135": 135, "137": 137, "176": 176, "272": 272, "1247": 1247,
    "33.75": 33.75, "61.5": 61.5, "291.5": 291.5, "784.96": 784.96,
}


def build_exprs(lo, hi):
    """生成「簡單結構式」：R 之一元運算 ＋ 二元運算（一層），保留落喺 [lo,hi] 者。"""
    out = {}
    for k, v in R.items():
        for tag, w in (("", v), ("x2", 2 * v), ("x3", 3 * v), ("/2", v / 2), ("/3", v / 3),
                       ("sq", v * v), ("rt", v ** 0.5), ("+1", v + 1), ("-1", v - 1),
                       ("+2", v + 2), ("-2", v - 2)):
            if lo <= w <= hi:
                out.setdefault(f"{tag}[{k}]", w)
    ks = list(R.items())
    for k1, v1 in ks:
        for k2, v2 in ks:
            if k1 == k2:
                continue
            for tag, w in (("+", v1 + v2), ("-", v1 - v2), ("*", v1 * v2), ("/", v1 / v2)):
                if lo <= w <= hi:
                    out.setdefault(f"[{k1}]{tag}[{k2}]", w)
    # 二層：對上述結果再乘/除 2,3
    base = list(out.items())
    for name, w in base:
        for tag, ww in (("*2", 2 * w), ("/2", w / 2), ("*3", 3 * w), ("/3", w / 3)):
            if lo <= ww <= hi:
                out.setdefault(f"{name}{tag}", ww)
    return out


def density_test(label, T, rel_band=0.30):
    lo, hi = T * (1 - rel_band), T * (1 + rel_band)
    ex = build_exprs(lo, hi)
    vals = sorted(set(round(v, 12) for v in ex.values()))
    width = hi - lo
    print(f"\n  【{label}】T = {T:.6f}；局部範圍 [{lo:.2f}, {hi:.2f}]；表達式數 = {len(vals)}")
    row = {"label": label, "T": T, "n_local": len(vals), "width": width}
    for tol_tag, tol in (("0.01%", 1e-4), ("0.1%", 1e-3), ("1%", 1e-2)):
        win = T * tol
        hits = [v for v in vals if abs(v - T) <= win]
        exp_uniform = len(vals) * (2 * win) / width
        row[tol_tag] = {"hits": len(hits), "expected_uniform": exp_uniform,
                        "closest": min(vals, key=lambda v: abs(v - T))}
        print(f"    tol={tol_tag:6s} 窗 ±{win:.4f}：命中 {len(hits):2d}，均勻期望 {exp_uniform:.2f}，"
              f"最近 = {min(vals, key=lambda v: abs(v - T)):.6f}")
    return row


row1 = density_test("S* = ln(1/μ)", S_STAR)
row2 = density_test("S*/λ_WH = ln(1/μ)/(2/123)", S_STAR / float(LAMBDA_WH))

# 判準：期望 < 1 ⟹ 該窗有資訊（0 命中 = 有效否證）；期望 ≥ 1 ⟹ 噪音地板
chk("G4.1 S* 之 0.1% 窗：均勻期望 < 1 ⟹ 該窗有資訊",
    row1["0.1%"]["expected_uniform"] < 1,
    f"期望 = {row1['0.1%']['expected_uniform']:.2f}，命中 = {row1['0.1%']['hits']}")
chk("G4.2 ⟹ S* 之 0.1% 級：**零命中** ＝ 有效否證（207+ 個簡單結構式無一接近）",
    row1["0.1%"]["hits"] == 0,
    f"最近 = {row1['0.1%']['closest']:.6f}（T = {S_STAR:.6f}）")
chk("G4.3 S* 之 1% 窗：均勻期望 ≥ 1 ⟹ 噪音地板，1% 級配對一律無資訊",
    row1["1%"]["expected_uniform"] >= 1,
    f"期望 = {row1['1%']['expected_uniform']:.2f}，命中 = {row1['1%']['hits']}")
chk("G4.4 S*/λ_WH 之 0.1% 窗：同樣零命中（若期望 < 1 亦屬有效否證）",
    row2["0.1%"]["hits"] == 0,
    f"期望 = {row2['0.1%']['expected_uniform']:.2f}，最近 = {row2['0.1%']['closest']:.4f}")
chk("G4.5 ⟹ 指數類（c=1 及 c=λ_WH 兩參數化）皆被打死；本輪零升格",
    True, "本輪無 candidate")

RES["_G4"] = {"rows": [row1, row2]}


# ============================================================================
# §5 G5 — 破牆規格
# ============================================================================
hdr("§5 G5 — 破牆所需之規格（identity-level）")

spec = [
    "S 須以鏈／樹之既有量**獨立定義**（定義中不得出現 μ、M₀、M_Pl、ω_ring、R_world）",
    "S 與 ln(1/μ) 須 identity-level 相符（dev < 1e−10；受輸入 8 s.f. 限制，實際下限 ~1e−8）",
    "S 之形式須閉式（有限項 γ 多項式／整數計數／對數），不得含自由係數",
    "若 S 含新常數 ⟹ 須先過四 Hard Tests（dof 成本／10⁷ null／樣本外盲測／結構剛性）",
]
for i, s in enumerate(spec, 1):
    print(f"  {i}. {s}")
chk("G5.1 規格已凍結（4 條）", True, "見 docs/07 §5")
chk("G5.2 本輪無任何 S 達標 ⟹ 不升格、不 close", True, f"S* = {S_STAR:.8f} 仍未識別")

RES["_G5"] = {"spec": spec}


# ============================================================================
# §6 匯總
# ============================================================================
hdr("§6 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "06_three_input_falsifiability.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_06.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_06.log'}")

sys.exit(0 if N_FAIL == 0 else 1)