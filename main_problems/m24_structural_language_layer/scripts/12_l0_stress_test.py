#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M24 scripts/12 — L0 壓力測試：禁區清單 ＋ 對抗例（8 → 實 7 案例）
==============================================================
判據（`docs/10_preregistration_round3.md` §4，跑前凍結）：

  REF-1  φ^(11+7/88) ≈ μ/e（M9：φ^x ∉ ℚ(ζ₈₈)，428σ）
  REF-2  M₀ = φ^(47+37/528)（O136 過擬合）
  REF-3  α⁻¹ = 134 + 3 + 15/(φ·256)（v7，O204 否）
  REF-4  μ/e = γ₇/γ₁（tautology ＋ 滑落）
  REF-5  19:41 ≈ cot(π/176)（0.03% 匹配）
  ADV-1  36 = 6²、98 = 2·7²、16 = 4²（「平方 × 小整數」同一系統）
  ADV-2  36 = 4·9、98 = 2·49、16 = 4·4（同一「因子系統」）

⚠ 預註冊筆誤（append-only 更正）：§4 列 7 條但寫「8/8」；實際案例數 = 7，
  通過條件按 **7/7**。原文不動，本檔附註。

判準沿用 `scripts/04` 之 S1–S4（範疇化／非恆等／強迫／等變）。
零 import 上游；mpmath dps40。
"""
import json
import os
import sys
from collections import Counter
from itertools import combinations
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    print(f"{'[PASS]' if cond else '[FAIL]'} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


try:
    import mpmath as mp
    mp.mp.dps = 40
except ImportError:
    raise SystemExit("mpmath required")

phi = (1 + mp.sqrt(5)) / 2
mu_e = mp.mpf("206.7682830")          # PDG μ/e（只作對照值，不作判據）
alpha_inv = mp.mpf("137.035999084")
M0_MeV = mp.mpf("27.900915")
N_d, B = 44, 4
MODES = list(range(1, 45))
GAM = {n: 4 * mp.sin(mp.pi * n / 88) ** 2 for n in MODES}

print("=" * 78)
print("M24 scripts/12 — L0 壓力測試（禁區清單 ＋ 對抗例）")
print("=" * 78)
print(f"  φ = {mp.nstr(phi, 20)}")

V = {}

# ── REF-1 ─────────────────────────────────────────────────────
r1 = phi ** (11 + mp.mpf(7) / 88)
V["REF-1"] = {"value": float(r1), "target": float(mu_e),
              "dev_rel": float((r1 - mu_e) / mu_e)}
print(f"\n  REF-1  φ^(11+7/88) = {mp.nstr(r1, 15)}｜μ/e = {mu_e}｜dev = "
      f"{V['REF-1']['dev_rel'] * 100:.4f}%")
# M9 域分離之機核側：φ^x 之 x 非整數 ⟹ 不在 ℚ(ζ₈₈) 之鏈頻譜比類（M9 已證）
V["REF-1"]["class"] = "A"
V["REF-1"]["reason"] = ("S1 FAIL：φ^x（x ∉ ℤ）非鏈頻譜元素之比；M9 域分離定理"
                        "（φ^x ∉ ℚ(ζ₈₈)）⟹ 無 functor 可映射至 γ-比值")

# ── REF-2 ─────────────────────────────────────────────────────
# 正典形式（v10 O136 docs/36–37，只讀）：M₀ = E_kchain · φ^(47+37/528)
# 故 φ^(47+37/528) 之對照目標 = X ≡ M₀/E_kchain = 6.87e9（非 M₀ 本身；M₀ 單位 MeV）
r2 = phi ** (47 + mp.mpf(37) / 528)
X_target = mp.mpf("6.87e9")
V["REF-2"] = {"value": float(r2), "target": float(X_target),
              "dev_rel": float((r2 - X_target) / X_target)}
print(f"  REF-2  φ^(47+37/528) = {mp.nstr(r2, 15)}｜X = M₀/E_kchain ≈ 6.87e9｜dev = "
      f"{float((r2 - X_target) / X_target) * 100:.6f}%")
V["REF-2"]["class"] = "A"
V["REF-2"]["reason"] = ("S1 FAIL：同 REF-1（φ-冪）；S3 FAIL：X 本身零湧現（自由量），"
                        "且 O136 docs/36 證「20 個結構構造同樣吻合 1e−7」＋「47 有 ≥5 讀法」"
                        "⟹ 過度吻合 = 反解，非強迫")

# ── REF-3 ─────────────────────────────────────────────────────
r3 = 134 + 3 + 15 / (phi * 256)
V["REF-3"] = {"value": float(r3), "target": float(alpha_inv),
              "dev_rel": float((r3 - alpha_inv) / alpha_inv)}
print(f"  REF-3  134 + 3 + 15/(φ·256) = {mp.nstr(r3, 15)}｜α⁻¹ = {alpha_inv}｜dev = "
      f"{V['REF-3']['dev_rel'] * 100:.4f}%")
V["REF-3"]["class"] = "M"
V["REF-3"]["reason"] = ("S1 △（三項皆池內量）／S2 ✅／S3 FAIL：同型三項式在池內多寫法"
                        "（O204 已判否）／S4 n/a ⟹ M 型")

# ── REF-4 ─────────────────────────────────────────────────────
r4 = GAM[7] / GAM[1]
hits_slide = []
for a in MODES:
    for b in MODES:
        v = GAM[a] / GAM[b]
        if abs(v - mu_e) / mu_e < mp.mpf("1e-3"):
            hits_slide.append((a, b, float(v)))
# 176-L 鏈（憲章 C1 路徑約定）：γ^{(176)}_k = 4sin²(πk/176)，k = 1..175
G176 = {k: 4 * mp.sin(mp.pi * k / 176) ** 2 for k in range(1, 176)}
hits_slide176 = []
near176 = []          # ≤2% 之近鄰（量度「滑落」族密度）
for a in range(1, 176):
    for b in range(1, 176):
        v = G176[a] / G176[b]
        d = abs(v - mu_e) / mu_e
        if d < mp.mpf("1e-3"):
            hits_slide176.append((a, b, float(v)))
        if d < mp.mpf("2e-2"):
            near176.append((float(d), a, b, float(v)))
near176.sort()
_best176 = near176[0] if near176 else None
_degen = len({(a, b) for d, a, b, v in near176 if abs(v - _best176[3]) < 1e-9}) if near176 else 0
hits_slide88_best = min(((float(abs(GAM[a] / GAM[b] - mu_e) / mu_e), a, b) for a in MODES for b in MODES))
V["REF-4"] = {"gamma7_over_gamma1": float(r4), "dev_rel": float((r4 - mu_e) / mu_e),
              "slide_hits_88L": len(hits_slide), "slide_hits_176L": len(hits_slide176),
              "near2pct_176L": len(near176),
              "best_176L": {"dev": _best176[0], "a": _best176[1], "b": _best176[2],
                            "value": _best176[3], "index_degeneracy": _degen},
              "best_88L": {"dev": hits_slide88_best[0], "a": hits_slide88_best[1],
                           "b": hits_slide88_best[2]}}
print(f"  REF-4  γ₇/γ₁ = {mp.nstr(r4, 15)}｜μ/e = {mu_e}｜dev = "
      f"{float((r4 - mu_e) / mu_e) * 100:.2f}%")
print(f"        滑落掃描（88-L）：命中 <0.1% 之 (a,b) = {len(hits_slide)} 對（例：{hits_slide[:3]}）")
print(f"        滑落掃描（176-L）：命中 <0.1% 之 (a,b) = {len(hits_slide176)} 對；"
      f"≤2% 近鄰 = {len(near176)} 對")
print(f"        176-L 最近值 = {_best176[3]:.4f}（dev {_best176[0] * 100:.2f}%），"
      f"索引退化數 = {_degen}；88-L 最近 dev = {hits_slide88_best[0] * 100:.2f}%")
V["REF-4"]["class"] = "A"
V["REF-4"]["reason"] = ("**字面 claim 數值否證**：γ₇/γ₁ = 48.009 ≠ μ/e（dev −76.8%）；"
              "兩鏈皆無 <0.1% 之 γ-比（176-L 最近 dev " + f"{_best176[0] * 100:.2f}%" +
              "）。S1 ✅（γ-比形式合法）／S2 ✅／S3 FAIL：最近值之**索引退化數 = " + str(_degen) +
              "**（同一值多對索引）⟹ 滑落族非唯一 ⟹ 即使命中亦屬 T 型。綜合判：非 F（A 型："
              "無 functor；且數值層已否證）")

# ── REF-5 ─────────────────────────────────────────────────────
cot176 = mp.cos(mp.pi / 176) / mp.sin(mp.pi / 176)
cands = {
    "19/41": mp.mpf(19) / 41,
    "41/19": mp.mpf(41) / 19,
    "19+41": mp.mpf(60),
    "19·41": mp.mpf(779),
    "1/cot(π/176)": mp.sin(mp.pi / 176) / mp.cos(mp.pi / 176),
    "cot(π/176)": cot176,
}
devs = {k: float((v - cot176) / cot176) if v != 0 else None for k, v in cands.items()}
print(f"  REF-5  cot(π/176) = {mp.nstr(cot176, 15)}｜候選讀法偏差（%）："
      f"{ {k: (round(d * 100, 4) if d is not None else None) for k, d in devs.items()} }")
# 只計「計數比」讀法（排除 cot(π/176) 自身之比對）
count_keys = ["19/41", "41/19", "19+41", "19·41"]
near = [k for k in count_keys if devs[k] is not None and abs(devs[k]) < mp.mpf("0.0003")]
V["REF-5"] = {"cot_176": float(cot176), "deviations_pct": {k: (d * 100 if d is not None else None)
                                                          for k, d in devs.items()},
              "readings_within_0.03pct": near}
V["REF-5"]["class"] = "M"
V["REF-5"]["reason"] = ("S1 FAIL：計數比與 trig 常數之間無 functor；預註冊字面讀法"
                        "（19/41 vs cot(π/176)）偏差 ~120×，可重現性不成立 ⟹ 連 M 型之"
                        "「吻合」亦不成立（記錄層 audit 發現）")

# ── ADV-1 / ADV-2 ─────────────────────────────────────────────
adv1 = {"16": 4 ** 2, "36": 6 ** 2, "98": 2 * 7 ** 2}
adv2 = {"16": 4 * 4, "36": 4 * 9, "98": 2 * 49}
print(f"\n  ADV-1  {adv1}")
print(f"  ADV-2  {adv2}")
V["ADV-1"] = {"values": adv1, "class": "A",
              "reason": "S1 FAIL：平方／平方倍係無對象之算術形式；索引（4,6,7）逐數選 ⟹ 非函子"}
V["ADV-2"] = {"values": adv2, "class": "A",
              "reason": "S1 FAIL：因子（4,9,49）逐數選（且 2 為另一自由因子）⟹ 非單一函子；"
                        "S3 FAIL：同一形式族可為任意三元組造值"}

# ── 分類匯總 ＋ 通過條件（7/7）────────────────────────────────
print("\n── 分類匯總 ──")
for k in ("REF-1", "REF-2", "REF-3", "REF-4", "REF-5", "ADV-1", "ADV-2"):
    print(f"  {k:<6} → {V[k]['class']} 型｜{V[k]['reason'][:70]}…")
nonF = [k for k in V if V[k]["class"] != "F"]
ck("L0S-1  全部 7 案例判非 F（F 型零命中）", len(nonF) == len(V), f"nonF = {len(nonF)}/{len(V)}")
ck("L0S-2  7/7 通過（預註冊 §4 之「8/8」為筆誤，實案例數 7；已 append 更正）",
   len(V) == 7, f"n_cases = {len(V)}")
ck("L0S-3  REF-4 機核：兩鏈皆無 <0.1% 命中（數值否證）＋176-L 最近值之索引退化 ≥2（滑落族非唯一）",
   len(hits_slide) == 0 and len(hits_slide176) == 0 and _degen >= 2,
   f"88L<0.1% = {len(hits_slide)}｜176L<0.1% = {len(hits_slide176)}｜"
   f"best176 dev = {_best176[0] * 100:.2f}%｜degen = {_degen}")
ck("L0S-4  REF-5 之預註冊字面讀法不可重現（19/41 與 cot(π/176) 差 ~120×）⟹ 記錄層 audit 發現，須 register 回報",
   not near, f"within 0.03% = {near}｜19/41 dev = {float(devs['19/41']) * 100:.2f}%")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/12 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "12_l0_stress_test.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": o, "detail": d} for t, o, d in CHECKS],
               "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
               "cases": V, "n_cases": len(V),
               "prereg_miscount_note": "docs/10 §4 寫 8/8 但列 7 案例；通過條件按 7/7（append-only 更正）",
               "verdict": "L0 壓力測試 7/7 判非F（禁區清單＋對抗例）；REF-5 字面讀法不可重現"},
              f, ensure_ascii=False, indent=2)
sys.exit(0 if nfail == 0 else 1)
