"""M5 docs/31 配套 — `a₀`、`t` 結構化之可行性（可識別性 no-go ＋ 預註冊掃描）
================================================================================
`docs/30` 已封死 contact 幾何；剩餘唯一全開之門 = 兩體尺度 `(a₀, t)` 有無結構式（M6 交接項）。
`docs/28` 校準：`ħc/a₀ = 0.513626232 GeV`、`t = 1.7711 MeV`；以 `M₀ = 0.027900914 GeV`（M6）除之得兩個無量綱數：

  u ≡ ħc/(a₀ M₀) = 18.408939…,    v ≡ t/M₀ = 0.063478206…

本輪唔係「搵條式出嚟湊」，而係先問：**喺可達精度下，封閉式能否唯一識別 (u, v)？**
做法 = 預註冊文法（只准 IDCM 已認證原子 ＋ 有限運算）＋ 全枚舉 ＋ **隨機目標零假設**（Monte Carlo 對照）。
若命中數 ≈ 隨機期望，則任何「命中」皆不可作為結構證據 ⟹ 封閉式路線 no-go，必須走算子推導。

判準（事前寫死）：
  W1 原子集只准已認證結構量（N_d、層計數、φ、π、√2/√3/√5、γ_n 譜、d_H、Σγ、F₁₀、L₇ 等），不得臨時加入。
  W2 表達式複雜度上限：≤ 3 個原子（含 1 個二元運算 ＋ 可選 √ 或整數次方 2,3,−1）。
  W3 命中定義：|expr/u − 1| ≤ 1e-3（0.1%）；另報 1e-4 嚴格命中。
  W4 顯著性：以 20,000 個隨機對數均勻目標（範圍 u ∈ [1,100]、v ∈ [1e-3,1e-1]）求每次期望命中數 E；
     須 obs/E ≥ 10 才算「非偶然」，否則記 no-go。
  W5 不得以命中為由改 atom 集重跑（p-hacking）；本腳本 atom 集寫死於碼內、單次枚舉。
輸出：`data/33_scale_identifiability.json`
"""
import itertools
import json

import numpy as np

HBARC = 0.1973269804
A0, T0 = 0.384184e-15, 1.7711e-3          # fm / GeV（v2 校準）
M0 = 0.027900914                            # GeV（M6）
U = HBARC / (A0 * 1e15) / M0                # ħc/(a₀ M₀)，fm→GeV
V = T0 / M0
ND, PHI = 88, (1 + 5 ** 0.5) / 2

print("=" * 104)
print("M5 docs/31 — (a₀, t) 封閉式可識別性（預註冊文法 ＋ 隨機目標零假設）")
print("=" * 104)
print(f"\n目標：u = ħc/(a₀M₀) = {U:.9f}   v = t/M₀ = {V:.9f}")
print(f"      （需 0.1% 精度 = ±{U*1e-3:.6f} / ±{V*1e-3:.9f}）")

# ---------------- W1 原子集（預註冊、寫死） ----------------
gam = [4 * np.sin(np.pi * n / ND) ** 2 for n in range(0, 45)]      # 鏈譜 γ_n = 4 sin²(πn/88)
atoms = {
    "N_d": ND, "N_d/2": ND / 2, "L7": 29, "F10": 55, "N_gen": 3, "B": 4, "b-1": 3,
    "Sig_g": 90, "d_H": np.log(3) / np.log(4), "pi": np.pi, "phi": PHI,
    "sqrt2": 2 ** 0.5, "sqrt3": 3 ** 0.5, "sqrt5": 5 ** 0.5, "e": np.e,
    "g44": gam[44], "g22": gam[22], "g11": gam[11], "g7": gam[7], "g1": gam[1],
    "mu2": 8.0, "sec88": 1 / np.cos(np.pi / ND), "F7": 13, "L6": 18, "L8": 47,
}
print(f"\n[W1] 原子集（{len(atoms)} 個，全部已認證）：{', '.join(atoms)}")

# ---------------- W2 表達式枚舉 ----------------
exprs = {}
for k, v in atoms.items():
    exprs[k] = v
keys = list(atoms)
for x, y in itertools.permutations(keys, 2):
    a, b = atoms[x], atoms[y]
    if b != 0:
        exprs[f"{x}/{y}"] = a / b
    exprs[f"{x}*{y}"] = a * b
    exprs[f"{x}+{y}"] = a + b
    exprs[f"{x}-{y}"] = a - b
for x in keys:
    a = atoms[x]
    for nm, val in (("sq", a * a), ("cu", a ** 3), ("inv", 1 / a if a else np.nan),
                    ("sqrt", np.sqrt(a) if a > 0 else np.nan)):
        exprs[f"{nm}({x})"] = val
# 三原子：a/(b*c)、a*b/c
for x, y, z in itertools.permutations(keys, 3):
    a, b, c = atoms[x], atoms[y], atoms[z]
    if b * c:
        exprs[f"{x}/({y}*{z})"] = a / (b * c)
    if c:
        exprs[f"{x}*{y}/{z}"] = a * b / c
print(f"[W2] 枚舉 {len(exprs)} 條（≤3 原子、預註冊運算）")


def hits(target, tol):
    return sorted([(abs(v / target - 1), k, v) for k, v in exprs.items()
                   if np.isfinite(v) and abs(v / target - 1) <= tol])


def mc_expected(target, lo, hi, tol, n=20000, seed=7):
    rng = np.random.default_rng(seed)
    tg = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    vals = np.array([v for v in exprs.values() if np.isfinite(v) and v > 0])
    vals.sort()
    cnt = 0
    for t in tg:
        lo_i = np.searchsorted(vals, t * (1 - tol), "left")
        hi_i = np.searchsorted(vals, t * (1 + tol), "right")
        cnt += hi_i - lo_i
    return cnt / n


print("\n[W3] 命中（0.1% 容差）")
res = {}
for nm, tgt, rng_ in (("u", U, (1, 100)), ("v", V, (1e-3, 1e-1))):
    h3 = hits(tgt, 1e-3)
    h4 = hits(tgt, 1e-4)
    E3 = mc_expected(tgt, *rng_, 1e-3)
    E4 = mc_expected(tgt, *rng_, 1e-4)
    print(f"\n  {nm} = {tgt:.9f}：0.1% 命中 {len(h3)} 條；1e-4 命中 {len(h4)} 條")
    print(f"     隨機目標期望 E(0.1%) = {E3:.1f} 條、E(1e-4) = {E4:.2f} 條")
    for d, k, val in h3[:6]:
        print(f"       {k:28s} = {val:.9f}  偏差 {d*100:.4f}%")
    if h4:
        print("     1e-4 嚴格命中：")
        for d, k, val in h4[:6]:
            print(f"       {k:28s} = {val:.9f}  偏差 {d*100:.5f}%")
    res[nm] = {"target": tgt, "hits_1e3": len(h3), "hits_1e4": len(h4),
               "E_1e3": E3, "E_1e4": E4,
               "top_1e4": [{"expr": k, "val": val, "dev_pct": d * 100} for d, k, val in h4[:8]],
               "signif_1e3": len(h3) / max(E3, 1e-9), "signif_1e4": len(h4) / max(E4, 1e-9)}

print("\n[W4] 判（W4 門檻：obs/E ≥ 10 才算非偶然）")
for nm, r in res.items():
    print(f"  {nm}：0.1% obs/E = {r['signif_1e3']:.2f}；1e-4 obs/E = {r['signif_1e4']:.2f}")
    ok = r["signif_1e3"] >= 10 or r["signif_1e4"] >= 10
    print(f"     ⟹ {'🟢 有非偶然命中' if ok else '🔴 no-go：命中數與隨機期望同量級 ⟹ 封閉式不可識別'}")

no_go = all(not (r["signif_1e3"] >= 10 or r["signif_1e4"] >= 10) for r in res.values())
print(f"\n[判決] {'🔴 封閉式路線 no-go（u、v 皆不可由預註冊文法識別）' if no_go else '🟢 至少一數有非偶然命中，須逐條審'}")
print("      ⟹ 含義：`(a₀, t)` 唔可以用「搵條式」方式結構化；唯一合法路線 = 由 K_chain⊗K_tree 直接推兩體算子")
print("         （即消除 (a, t) 而非解釋它們），此為 docs/26 所記之『dof=1 之外的第二缺口』。")

with open("data/33_scale_identifiability.json", "w") as f:
    json.dump({"U": U, "V": V, "M0": M0, "n_atoms": len(atoms), "n_exprs": len(exprs),
               "results": res, "no_go": bool(no_go)}, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/33_scale_identifiability.json")