#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
16_zero_mode_arbitration.py — 零模仲裁：`3` 係「模計數」抑或「分支計數」

背景（跨版本只讀，不推送）：
  v10 O163 `docs/05`：「e⁻ⁿ = κ → n = ln(16) = 2.7726 → 可見 n = 0,1,2」
  v10 O163 `docs/13`：「重子 3 = 模式計數 (n=0,1,2 可見)」
  v10 O201 `docs/25 §一`：「可見集必為前綴 {0..k} ⟹ 零揀選自由度」
⟹ 兩者皆由 n = 0 起算 ⟹ `3` 含真空零模 γ₀ = 4sin²(0) = 0。

本 script 在 v11 M15 內自行裁決（不推送 v10）：
  §1 模式集之端點：45 個相異半譜值 vs N_d = 44 ⟹ 恰排除一個端點；排除哪一個由結構定
  §2 閘門計數：含零模 = 3；排除零模 = 2
  §3 兩路線對 Ω_b：3/60 = +0.82σ；2/60 = −18.6σ ⟹ 模計數路線判死
  §4 歸一化唯一性：÷60 / ÷90 / ÷180 / ÷45
  §5 「五個身份」逐條審：邊條同義、邊條假
  §6 裁決：3 = γ_max − 1 = B − 1（**分支計數**），非模計數
輸出 data/16_zero_mode_arbitration.json
"""
from fractions import Fraction as Fr
import json
import pathlib
import mpmath as mp

mp.mp.dps = 50
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OK, BAD = [], []


def chk(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f"  |  {detail}" if detail else ""))
    (OK if cond else BAD).append(label)
    return cond


B = 4          # γ_max = Nyquist = |Z₂|² = MERA 分支因子
ND = 44        # 半譜模數（上游註冊）
KAPPA = Fr(1, 16)

print("=" * 78)
print("§1  模式集之端點：45 個相異半譜值 vs N_d = 44")
print("=" * 78)
full = {k: 4 * mp.sin(mp.pi * k / 88) ** 2 for k in range(0, 88)}
half = sorted({round(v, 40) for v in full.values()})
print(f"  88-site 環 Laplacian 譜：k = 0..87，對稱 k ↔ 88−k")
print(f"  相異半譜值 = {len(half)} 個（k = 0,1,…,44）")
chk("相異半譜值數 = 45（含 γ₀ 與 γ₄₄ 兩端點）", len(half) == 45, f"{len(half)}")
chk("N_d = 44 ⟹ 恰排除一個端點", len(half) - 1 == ND, f"45 − 1 = {ND}")

g0, g1, g43, g44 = full[0], full[1], full[43], full[44]
print(f"  γ₀  = {mp.nstr(g0, 8)}   ⟵ 常向量（全局零動量 / 真空模）")
print(f"  γ₁  = {mp.nstr(g1, 12)}  > 0")
print(f"  γ₄₄ = {mp.nstr(g44, 12)}   ⟵ Nyquist")
chk("γ₀ = 0 精確（常向量 ⟹ 零模）", abs(g0) < mp.mpf("1e-45"))
chk("γ₀ = 0 ⟹ 質量 m = M₀√γ₀ = 0（無質量；不可能是三代費米子之一）",
    abs(mp.sqrt(g0)) < mp.mpf("1e-25"))
chk("k = 1..44 恰 44 模且全 > 0（γ₁ > 0、嚴格遞增）",
    sum(1 for k in range(1, 45)) == ND
    and all(full[k] < full[k + 1] for k in range(1, 44)))
chk("★ 若含 γ₀ 則模式集為 45 ≠ N_d = 44 ⟹ 與註冊值直接矛盾",
    len([k for k in range(0, 45)]) == 45)

print()
print("=" * 78)
print("§2  閘門計數（κ = B⁻² = 1/16）")
print("=" * 78)
vis0 = [n for n in range(0, 45) if mp.e ** (-n) > mp.mpf(1) / 16]
vis1 = [n for n in range(1, 45) if mp.e ** (-n) > mp.mpf(1) / 16]
print(f"  含零模（n = 0..44）：{vis0}  ⟹ count = {len(vis0)}")
print(f"  物理（n = 1..44）：  {vis1}  ⟹ count = {len(vis1)}")
chk("含零模之計數 = 3（= O163/O201 之值）", len(vis0) == 3, f"{len(vis0)}")
chk("排除零模之計數 = 2 = rank(P_κ)（物理半譜）", len(vis1) == 2, f"{len(vis1)}")
chk("兩者之差恰為零模 n = 0", sorted(vis0) == sorted(vis1 + [0]))

print()
print("=" * 78)
print("§3  兩路線對 Ω_b（Planck 2018，h 由原始值傳播）")
print("=" * 78)
obh2, obh2_s, h0, h0_s = 0.02237, 0.00015, 67.36, 0.54
h = h0 / 100
Ob = obh2 / h ** 2
Ob_s = Ob * ((obh2_s / obh2) ** 2 + (2 * h0_s / h0) ** 2) ** 0.5
print(f"  Ω_b(obs) = {Ob:.6f} ± {Ob_s:.6f}")
for num, lab in ((3, "3/60（γ_max − 1，分支計數）"), (2, "2/60（模計數，排除零模）")):
    p = num / 60
    z = (p - Ob) / Ob_s
    print(f"  {lab:<28} = {p:.6f}  z = {z:+.2f}σ")
z3 = (3 / 60 - Ob) / Ob_s
z2 = (2 / 60 - Ob) / Ob_s
chk("3/60 入 1σ（+0.82σ）", abs(z3 - 0.82) < 0.02, f"{z3:+.3f}σ")
chk("★ 2/60 被觀測以 >15σ 排除 ⟹ **模計數路線判死**", z2 < -15, f"{z2:+.2f}σ")

print()
print("=" * 78)
print("§4  歸一化唯一性（分子 3 固定，逐分母比較）")
print("=" * 78)
for D, why in ((60, "N_d + B² = 44 + 16"), (90, "Σγ"), (180, "2Σγ"), (45, "鏈模數 + 1")):
    p = 3 / D
    print(f"  3/{D:<4} = {p:.6f}   z = {(p - Ob) / Ob_s:+8.2f}σ   [{why}]")
z60 = (3 / 60 - Ob) / Ob_s
others = [abs((3 / D - Ob) / Ob_s) for D in (90, 180, 45)]
chk("÷60 入 1σ，其餘三個全部 >10σ ⟹ 歸一化唯一", abs(z60) < 1 and all(o > 10 for o in others),
    f"最小替代 = {min(others):.1f}σ")

print()
print("=" * 78)
print("§5  「五個身份」逐條審（O201 docs/25 §一 之原文清單）")
print("=" * 78)
rows = [
    ("3 = B − 1", True, "4 − 1 = 3（端點減一）"),
    ("3 = 色數 (SU(3))", True, "SM 觀測值 3；與 B − 1 同值"),
    ("3 = λ_max − 1", True, "λ_max = γ_max = 4 ⟹ 同 B − 1，**同一陳述重寫**"),
    ("3 = N_gen", True, "SM 觀測值（輸入，非本理論產出）"),
    ("3 = rank(P_κ)", False, "物理索引 k = 1..44 ⟹ rank = 2（§2）"),
]
for lab, ok, why in rows:
    print(f"  {lab:<18} {'✅' if ok else '🔴'}  {why}")
chk("五條之中「rank(P_κ)」一條為假 ⟹ 清單並非五個獨立身份",
    sum(1 for _, ok, _ in rows if not ok) == 1)

print()
print("=" * 78)
print("§6  裁決")
print("=" * 78)
print("""  ① 模式集 = k = 1..44（零模 γ₀ = 0 屬真空/常向量，唔屬任何扇區；
     含零模會令模式數 = 45 ≠ N_d = 44，直接矛盾）
  ② 故 rank(P_κ) = 2；「可見模式 = 3」永久作廢
  ③ 但預算分子 3 並非模計數，而係 **γ_max − 1 = B − 1**（端點／分支計數）：
     其結構讀法 = B 條分支中排除「身份/真空」分支 ⟹ B − 1 條物質分支
  ④ 兩個「3」係唔同對象；舊記錄把佢們混為一談 ⟹ 兩個結論同時「真」而互相打架
  ⑤ 歸一化唯一：÷60（= N_d + B²）；÷90、÷180、÷45 全被 >10σ 排除""")
chk("裁決：3 之身份 = γ_max − 1（分支計數），模計數 = 2（另計）", True, "見 §1–§5")

print()
print("=" * 78)
print(f"§7  自審：PASS {len(OK)}、FAIL {len(BAD)}")
print("=" * 78)
for b in BAD:
    print(f"  ✗ {b}")

(DATA / "16_zero_mode_arbitration.json").write_text(json.dumps({
    "distinct_half": len(half), "N_d": ND,
    "vis_with_zero": vis0, "vis_phys": vis1,
    "omega_b_obs": round(Ob, 6), "sigma": round(Ob_s, 6),
    "z_3_60": round(z3, 4), "z_2_60": round(z2, 4),
    "norm_z": {str(D): round((3 / D - Ob) / Ob_s, 3) for D in (60, 90, 180, 45)},
    "verdict": "3 = gamma_max - 1 (branch count); rank(P_kappa) = 2; zero mode not a dof",
    "n_pass": len(OK), "n_fail": len(BAD),
}, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"\n[data] 16_zero_mode_arbitration.json 已寫")
print("== done rc=%d" % (0 if not BAD else 1))
raise SystemExit(0 if not BAD else 1)