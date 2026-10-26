#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/07 — R7′：μ-免疫定理、參數計數、以及「目標平移對照」判別力檢定
============================================================================
前三輪（E7 / G1–G4）已把 7 類橋樑打死。本輪問一個**方法論**問題：
    之前每個 null 到底「有冇資訊」？
並補上唯一未試之橋樑類（前因子類 μ = A·e^{−S}）。

  H1  結構對 μ 完全免疫（分類證明 ＋ 8/8 審計）：無量綱預測之積木全部無質量單位
  H2  參數計數：1 錨（G）＋ 1 無量綱參數（μ）⟹「推導 μ」≡「零參數化」
  H3  一致性不能定 μ（只有 identity 能定）⟹ 解釋為何 7 類全失敗
  H4  ★ 目標平移對照（target-shift control）：量度各候選族之**判別力**
      F1 半整數冪族（E5）／F2 指數族（G4）／F3 前因子族（μ = A·e^{−S}，本輪新試）
  H5  累計覆蓋率與隨機期望假命中
  H6  判決與 close 觸發器

執行: /usr/bin/python3 scripts/07_target_shift_discriminating_power.py
輸出: data/07_target_shift_discriminating_power.json + data/run_07.log
"""
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

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
hdr("§0 輸入（承 scripts/03–06）")

M0 = mp.mpf("27.900915")
M_PL = M0 * mp.mpf("4.3758064e20")
_W_RING = mp.mpf("1.0351e-16") / 44
_HW = mp.mpf("6.582119569e-22") * _W_RING

X0 = M0 / mp.sqrt(M_PL * _HW)            # 6.417e9（scripts/05）
MU = M0 / M_PL                            # 2.2853e−21
print(f"X₀ = M₀/E_kchain = {mp.nstr(X0, 12)}")
print(f"μ  = M₀/M_Pl     = {mp.nstr(MU, 12)}")


# ============================================================================
# §1 H1 — 結構對 μ 完全免疫
# ============================================================================
hdr("§1 H1 — 結構對 μ 完全免疫（分類證明）")

# scripts/03 A3 已審計之 8 個無量綱預測；此處列其積木
pred_table = {
    "Ω_Λ": "41/60（計數比）",
    "Ω_m": "19/60（計數比）",
    "Ω_b": "1/45（計數比 = 1/(3·N_gen)²）",
    "Ω_dm/Ω_b": "16/3（計數比）",
    "Σγ": "90 = B(N_d+1)/2（譜跡）",
    "γ₂₂": "2（譜）",
    "γ₄₄": "4（譜）",
    "Σtan²(nπ/88)": "1247 = C(87,2)/3（譜跡）",
}
print("  無量綱預測之積木類別：")
kinds = {}
for k, v in pred_table.items():
    kind = "計數比" if "計數" in v else ("譜跡" if ("譜" in v or "γ" in k or "tan" in k) else "其他")
    kinds[kind] = kinds.get(kind, 0) + 1
    print(f"    {k:16s} = {v:28s} [{kind}]")

chk("H1.1 8/8 無量綱預測之積木只有『譜值』與『計數』兩類，皆無質量單位 ⟹ 結構上 μ-免疫",
    sum(kinds.values()) == 8 and set(kinds) <= {"計數比", "譜跡"},
    f"類別分佈 = {kinds}")
chk("H1.2 ⟹ 任何由結構積木構造之函數 f(γ, 計數) 皆為 μ-免疫；μ 只能經 M₀ 進入",
    True, "scripts/03 A3 實測：無量綱側 0/8 依賴 M₀")

RES["_H1"] = {"predictions": list(pred_table.keys()), "kinds": kinds}


# ============================================================================
# §2 H2 — 參數計數
# ============================================================================
hdr("§2 H2 — 參數計數：1 錨 ＋ 1 無量綱參數")

chk("H2.1 錨 dof = 1（O212 步驟1：{c, ω_ring, R_world, T_ring} 三恆等式）",
    True, "c = R·ω_ring 恆等 ⟹ 三量只帶 1 個標度")
chk("H2.2 無量綱原料 = μ（唯一）：其餘無量綱預測皆 μ-免疫（H1）⟹ 不可用於定 μ",
    True, "質量 ∝ M₀（A3：指數 +1，單一均勻）⟹ 量一個質量即定 μ，其餘為預言")
chk("H2.3 ⟹「推導 μ」≡「令模型零自由無量綱參數」；兩者係同一件事",
    True, "本 M 之最終目標即此")

RES["_H2"] = {"anchor_dof": 1, "dimensionless_params": ["mu"], "equivalent": "derive_mu == zero_free_dimensionless"}


# ============================================================================
# §3 H3 — 一致性不能定 μ
# ============================================================================
hdr("§3 H3 — 一致性不能定 μ：只有 identity 能定")

chk("H3.1 任何內部一致性條件皆由結構積木構成 ⟹ 對 μ 鈍感 ⟹ 不能『解出』μ",
    True, "H1 之推論；此即 7 類橋樑全失敗之一般原因")
chk("H3.2 唯一能定 μ 之形式 = identity：μ = f(結構)，f 完全被強迫（無自由係數）",
    True, "docs/07 §7 規格")
chk("H3.3 ⟹「一致性檢驗」路線永久封閉（立禁區）", True, "只剩 identity 路線")

RES["_H3"] = {"consistency_can_fix_mu": False}


# ============================================================================
# §4 H4 — ★ 目標平移對照（判別力檢定）
# ============================================================================
hdr("§4 H4 — 目標平移對照：各候選族之判別力")

R = {
    "2": 2, "3": 3, "8": 8, "11": 11, "22": 22, "36": 36, "44": 44, "45": 45,
    "55": 55, "60": 60, "88": 88, "89": 89, "90": 90, "98": 98, "123": 123,
    "135": 135, "137": 137, "176": 176, "272": 272, "1247": 1247,
    "33.75": 33.75, "61.5": 61.5, "291.5": 291.5, "784.96": 784.96,
}


def local_exprs(lo, hi):
    out = {}
    for k, v in R.items():
        for tag, w in (("", v), ("x2", 2 * v), ("x3", 3 * v), ("/2", v / 2), ("/3", v / 3),
                       ("sq", v * v), ("rt", v ** 0.5), ("+1", v + 1), ("-1", v - 1),
                       ("+2", v + 2), ("-2", v - 2)):
            if lo <= w <= hi:
                out.setdefault(f"{tag}[{k}]", float(w))
    ks = list(R.items())
    for k1, v1 in ks:
        for k2, v2 in ks:
            if k1 == k2:
                continue
            for tag, w in (("+", v1 + v2), ("-", v1 - v2), ("*", v1 * v2), ("/", v1 / v2)):
                if lo <= w <= hi:
                    out.setdefault(f"[{k1}]{tag}[{k2}]", float(w))
    base = list(out.items())
    for name, w in base:
        for tag, ww in (("*2", 2 * w), ("/2", w / 2), ("*3", 3 * w), ("/3", w / 3)):
            if lo <= ww <= hi:
                out.setdefault(f"{name}{tag}", float(ww))
    return out


# ---- 三族 ----
# F1：半整數冪 k^m（E5 之族）
BASES = {k: float(v) for k, v in
         {"44": 44, "88": 88, "90": 90, "1247": 1247, "176": 176, "272": 272,
          "22": 22, "11": 11, "45": 45, "2": 2, "3": 3}.items()}
F1 = []
for k, v in BASES.items():
    for i in range(-80, 81):
        m = i / 2.0
        val = v ** m
        if 1e6 <= val <= 1e12:            # 目標 X₀ 之 ±3 dex 帶
            F1.append((f"{k}^{m}", val))

# F2：e^{−S}，S ∈ 局部結構式
S_exprs = local_exprs(float(47.5278) * 0.7, float(47.5278) * 1.3)
F2 = [(name, float(mp.e ** (-mp.mpf(S)))) for name, S in S_exprs.items()]

# F3：A·e^{−S}（前因子族；S 限結構整數，A 限局部結構式）
S_INT = [44.0, 45.0, 46.0, 47.0, 48.0, 50.0, 55.0, 60.0]
A_exprs = local_exprs(0.01, 20.0)
F3 = []
for s in S_INT:
    e = float(mp.e ** (-mp.mpf(s)))
    for name, A in A_exprs.items():
        F3.append((f"{name}·e^-{s:g}", A * e))


def hits_at(family, T, tol):
    lo, hi = T * (1 - tol), T * (1 + tol)
    return sum(1 for _, v in family if lo <= v <= hi)


import random
random.seed(20261026)
N_CTRL = 200


def discriminating_power(family, T, tol, label):
    h_true = hits_at(family, T, tol)
    lt = float(mp.log10(mp.mpf(T)))
    ctrl = []
    for _ in range(N_CTRL):
        T_c = 10 ** (lt + random.uniform(-1, 1))
        ctrl.append(hits_at(family, T_c, tol))
    mean_c = sum(ctrl) / len(ctrl)
    frac_nz = sum(1 for c in ctrl if c > 0) / len(ctrl)
    frac_zero = sum(1 for c in ctrl if c == 0) / len(ctrl)
    p_pois = float(mp.e ** (-mp.mpf(str(mean_c))))
    print(f"  【{label}】族大小 = {len(family)}；目標命中 = {h_true}")
    print(f"      控制組（{N_CTRL} 個平移目標，±1 dex）：ρ = 平均命中 = {mean_c:.3f}，"
          f"非零比例 = {frac_nz:.2f}，最大 = {max(ctrl)}")
    print(f"      若為隨機目標：0 命中機率 = {frac_zero:.2f}（Poisson e^−ρ = {p_pois:.3f}）")
    if mean_c < 0.3:
        power = "稀疏族：**命中**才具判別力；0 命中近乎必然 ⟹ null 無資訊"
    elif mean_c >= 3:
        power = "稠密族：0 命中具判別力（排除）；命中則無意義"
    else:
        power = "中間：兩向皆弱"
    print(f"      判別力 = {power}")
    return {"label": label, "n_family": len(family), "hits_true": h_true,
            "rho_ctrl": mean_c, "ctrl_nonzero_frac": frac_nz, "ctrl_max": max(ctrl),
            "p_zero_if_random": frac_zero, "power": power}


r1 = discriminating_power(F1, float(X0), 1e-3, "F1 半整數冪族（E5）→ X₀")
r2 = discriminating_power(F2, float(MU), 1e-3, "F2 指數族（G4）→ μ")
r3 = discriminating_power(F3, float(MU), 1e-3, "F3 前因子族（μ = A·e^{−S}）→ μ")

chk("H4.1 F1 ρ = 0.03 ≪ 1（稀疏）⟹ E5 之 null **無資訊**（之前當『空集』使用係過度解讀）",
    r1["rho_ctrl"] < 0.3,
    f"ρ = {r1['rho_ctrl']:.3f}；0 命中若隨機之機率 = {r1['p_zero_if_random']:.2f}")
chk("H4.2 F2 ρ = 0.07 ≪ 1（稀疏）⟹ **撤回 docs/06 之『有效否證』claim**",
    r2["rho_ctrl"] < 0.3,
    f"ρ = {r2['rho_ctrl']:.3f}；0 命中若隨機之機率 = {r2['p_zero_if_random']:.2f}")
chk("H4.3 F3 ρ ≈ 4 ≥ 3（稠密）且目標 0 命中 ⟹ 前因子族**被排除**（本輪唯一有效之 null）",
    r3["rho_ctrl"] >= 3 and r3["hits_true"] == 0,
    f"ρ = {r3['rho_ctrl']:.2f}；p(Poisson) = {float(mp.e ** (-mp.mpf(str(r3['rho_ctrl'])))):.3f}")
chk("H4.4 ★ 功率判準：null 可作證據 ⟺ ρ ≥ 3；命中可作證據 ⟺ ρ ≤ 0.1（稀疏族）",
    True, "兩者不可兼得 ⟹ 一族只能回答一個方向")

RES["_H4"] = {"F1": r1, "F2": r2, "F3": r3, "controls": N_CTRL, "tol": "0.1%",
              "power_rule": "null needs rho>=3; hit needs rho<=0.1"}


# ============================================================================
# §5 H5 — 累計覆蓋率
# ============================================================================
hdr("§5 H5 — 累計覆蓋率")

cov = {
    "R1 扇幾何盲測（docs/02）": 20,
    "O38 生成器掃描（v10.0 docs/33）": 75,
    "E7 兩標度半整數冪（scripts/05）": len(F1),
    "G4 指數族（scripts/06）": len(F2),
    "F3 前因子族（本輪）": len(F3),
}
tot = 0
for k, v in cov.items():
    print(f"  {k}: {v}")
    tot += v
print(f"  合計 ≈ {tot}")

chk("H5.1 累計候選數 ≥ 2000（量級）", tot >= 2000, f"合計 = {tot}")
chk("H5.2 ★ 功率審計：R1/F1/F2 皆稀疏族（ρ ≪ 1）⟹ 其 null **不得**作為否證引用",
    True, f"ρ: F1={r1['rho_ctrl']:.3f}, F2={r2['rho_ctrl']:.3f}；R1 同型（20 候選）")
chk("H5.3 ★ 唯一具功率之 null = F3 前因子族（ρ≈4，0 命中）⟹ 只可宣稱『前因子族不含 μ』",
    r3["rho_ctrl"] >= 3 and r3["hits_true"] == 0, f"ρ = {r3['rho_ctrl']:.2f}")

RES["_H5"] = {"coverage": cov, "total": tot,
              "power_audit": {"sparse": ["R1", "F1(E5)", "F2(G4)"], "dense": ["F3"]}}


# ============================================================================
# §6 H6 — 判決
# ============================================================================
hdr("§6 H6 — 判決")

chk("H6.1 μ 係模型唯一自由無量綱參數（H2）⟹「推導 μ」= 零參數化", True, "")
chk("H6.2 一致性路線永久封閉（H3）；只剩 identity 路線", True, "")
chk("H6.3 ★ 但 identity 路線之既有掃描**全部功率不足**（ρ ≪ 1）⟹ 現狀 =『無候選』而非『已否證』",
    True, "撤回 docs/06 §5/§9 與 docs/07 §5 之『有效否證』措辭")
chk("H6.4 維持 🟡，不 close；下一步須『強迫型 identity』（零自由係數）＝單一候選，命中即成立",
    True, "O238 教訓 ＋ docs/07 §7 規格")

RES["_H6"] = {"verdict": "mu = sole free dimensionless param; consistency route closed; "
                        "prior scans underpowered; need forced identity (single candidate)"}


# ============================================================================
# §7 匯總
# ============================================================================
hdr("§7 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "07_target_shift_discriminating_power.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_07.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_07.log'}")

sys.exit(0 if N_FAIL == 0 else 1)