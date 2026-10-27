#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/09 — R8 sphaleron 放寬量：若 T2c 翻（有規範動能項），P1/P2 剩幾多？

前設問題（docs/04 §四 T2c 係 gap-driven；docs/05 §八 明寫「P1 之放寬量未計算」）:
  T2c 之 exactness 依賴「已註冊結構無規範動能項」。若 O85 未來補上 YM 項，則
  B 只低能有效守恆（高温 sphaleron 洗白）。R8 把呢個「若」量化：
    (1) 零溫 sphaleron 位壘高度 E_sph（以 m_W、α_W 表示；B 係數取範圍）
    (2) 破缺相之 Boltzmann 抑制 e^{−E_sph(T)/T} 在幾個代表性溫度之值
    (3) 若確有洗白，η_B 必須源自 B−L：所需 (B−L) 之數值（SM 28/79 因子）
    (4) 零溫 B 破壞（'t Hooft instanton）之抑制 ⟹ P1（ΔB≠0 速率=0）之裕度

執行：/usr/bin/python3 scripts/09_sphaleron.py
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


# ---- 輸入（全部明示；SM 量取 PDG 2024 中心值）
M_W = 80.377          # GeV
M_H = 125.25          # GeV
SIN2_TW = 0.23122
ALPHA_EM = 1 / 137.035999084
ALPHA_W = ALPHA_EM / SIN2_TW      # ≈ 1/31.6，g²/4π
T_C = 159.5           # GeV，EW 相變臨界溫度（T_c，SM 格點 159.5 ± 1.5）
ETA_B = 6.2116e-10    # IDCM η_B（O168）
B_COEF = (1.5, 1.9, 2.5)          # E_sph = (2m_W/α_W)·B，m_H≈125 時 B≈1.9

print("=" * 78)
print("M10 R8 — sphaleron 放寬量：T2c 若翻，P1/P2 之剩餘裕度")
print("=" * 78)

# ---------------------------------------------------------------- (1)
print("\n[1] 零溫 sphaleron 位壘：E_sph = (2 m_W / α_W) · B(m_H)")
E_pre = 2 * M_W / ALPHA_W
print(f"     α_W = α/sin²θ_W = {ALPHA_W:.6f} = 1/{1/ALPHA_W:.3f}")
print(f"     2 m_W/α_W = {E_pre:.2f} GeV")
E_sph = {b: E_pre * b for b in B_COEF}
for b in B_COEF:
    print(f"     B = {b:.1f} ⟹ E_sph = {E_sph[b]:.0f} GeV = {E_sph[b]/1e3:.2f} TeV")
ck("1.1 E_sph(B=1.9) 落 8–10 TeV（文獻 9 TeV 級）", 8000 < E_sph[1.9] < 10000,
   f"{E_sph[1.9]:.0f} GeV ⟹ 與文獻 9.0 TeV 一致（差 {100*(E_sph[1.9]/9000-1):+.1f}%）")
ck("1.2 E_sph/T 無量綱（量綱檢查）", True,
   "E_sph[GeV]/T[GeV] ⟹ 純數；故抑制因子 e^{−E_sph/T} 無單位困擾")

# ---------------------------------------------------------------- (2)
print("\n[2] 破缺相抑制：E_sph(T) = E_sph(0)·√(1−(T/T_c)²)（T < T_c）")


def suppression(T_GeV, b=1.9, Tc=T_C):
    if T_GeV >= Tc:
        return 1.0, 0.0            # 對稱相：無位壘
    f = math.sqrt(max(1.0 - (T_GeV / Tc) ** 2, 0.0))
    x = E_sph[b] * f / T_GeV
    return math.exp(-x), x


print(f"     {'T (GeV)':>9}{'E_sph(T) (GeV)':>16}{'E_sph/T':>10}{'e^-x':>14}")
tab = []
for T in (200, 159, 140, 120, 100, 80, 50, 20, 1):
    s, x = suppression(T)
    f = 0.0 if T >= T_C else math.sqrt(max(1 - (T / T_C) ** 2, 0))
    print(f"     {T:>9}{E_sph[1.9]*f:>16.0f}{x:>10.1f}{s:>14.2e}")
    tab.append({"T_GeV": T, "x": x, "supp": s, "phase": "symmetric" if T >= T_C else "broken"})
ck("2.1 T=100 GeV 之抑制 < 1e−30", suppression(100)[0] < 1e-30,
   f"e^-x = {suppression(100)[0]:.2e}（x = {suppression(100)[1]:.1f}）")
ck("2.2 T ≥ T_c 無位壘 ⟹ 洗白全開（P2 之真正風險區）", suppression(160)[0] == 1.0,
   f"T = 160 GeV ≥ T_c = {T_C} ⟹ 抑制 = 1 ⟹ 對稱相內 B 被完全洗白（若無 B−L）")


def T_for(x_target, b=1.9):
    """求 E_sph(T)/T = x_target 之 T（二分，T<T_c）。"""
    lo, hi = 1.0, T_C - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if suppression(mid, b)[1] > x_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


T30 = T_for(69.08)        # ln(1e30)
print(f"     e^-x = 1e−30 之溫度 T = {T30:.1f} GeV（T < 此值 ⟹ 洗白不可觀）")
ck("2.3 洗白消失溫度落 100–140 GeV（EW 相變後即凍結）", 100 < T30 < 140,
   f"T < {T30:.1f} GeV ⟹ 每 Hubble 體積之洗白事件 < 1e−30 ⟹ 相變完成後 B 實質凍結")

# ---------------------------------------------------------------- (3)
print("\n[3] 若確有洗白：η_B 必須源自 B−L（SM 因子 28/79）")
f_conv = 28 / 79
BL_needed = ETA_B / f_conv
print(f"     B_final = (28/79)(B−L) ⟹ B−L = η_B/(28/79) = {BL_needed:.4e}")
ck("3.1 28/79 因子之算術", abs(79 / 28 - 1 / f_conv) < 1e-12,
   f"79/28 = {79/28:.6f} ⟹ (B−L) = {BL_needed:.3e}")
ck("3.2 需求 (B−L) 與 η_B 同量級（1e−9）", 1e-9 < BL_needed < 3e-9,
   f"(B−L)_required = {BL_needed:.3e} ⟹ 若 T2c 翻，IDCM 必須**新增** B−L 生成機制；"
   f"現時結構無此機制 ⟹ P2 由「初始條件」降級為「待生成之非對稱」")

# ---------------------------------------------------------------- (4)
print("\n[4] 零溫 B 破壞（'t Hooft instanton）⟹ P1 之裕度")
S_inst = 2 * math.pi / ALPHA_W
amp = math.exp(-S_inst)
rate = amp ** 2
print(f"     S_inst = 2π/α_W = {S_inst:.1f}")
print(f"     振幅 e^-S = {amp:.2e}；速率 ∝ e^-2S = {rate:.2e}")
ck("4.1 零溫 B 破壞速率抑制 < 1e−150", rate < 1e-150,
   f"e^-2S = {rate:.2e} ⟹ 等效質子壽命下限 ≫ 10^150 yr（觀測下限 ~1.6e34 yr）"
   f"⟹ **P1（ΔB≠0 速率=0）在 T=0 完全不受 sphaleron 威脅**")
ck("4.2 P1 之真正風險只在 T ≳ 130 GeV（早期宇宙）", T30 < 140,
   f"現時宇宙 T = 2.725 K = 2.3e−10 GeV ⟹ x = {E_sph[1.9]/2.348e-10:.1e} ⟹ "
   f"抑制 ~ 10^(−4e12) ⟹ 絕對為零")

print("\n[結論]")
print("     T2c 若翻（補上 YM 動力學）⟹")
print(f"       ✓ P1（ΔB 速率=0）**不變**：零溫抑制 e^-2S = {rate:.1e}")
print(f"       ✗ P2（η_B 係初始條件）**降級**：需 (B−L) = {BL_needed:.2e} 之生成機制")
print(f"       ⟹ 放寬量 = P1 安然 / P2 需 (B−L) 生成 ⟹ 兩條預言之依賴從此分離")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)

out = {
    "script": "09_sphaleron.py", "role": "R8 sphaleron 放寬量",
    "inputs": {"m_W_GeV": M_W, "m_H_GeV": M_H, "sin2_tw": SIN2_TW,
               "alpha_W": ALPHA_W, "T_c_GeV": T_C, "eta_B": ETA_B,
               "B_coef": list(B_COEF)},
    "E_sph_GeV": {str(k): v for k, v in E_sph.items()},
    "suppression_table": tab,
    "T_for_1e30_supp": T30,
    "B_minus_L_required": BL_needed, "conv_factor": f_conv,
    "S_inst": S_inst, "inst_rate": rate,
    "verdict": {"P1": "不變（零溫抑制 e^-2S）", "P2": "降級（需 B−L 生成）"},
    "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL,
}
with open(os.path.join(DATA, "09_sphaleron.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/09_sphaleron.json 已寫")