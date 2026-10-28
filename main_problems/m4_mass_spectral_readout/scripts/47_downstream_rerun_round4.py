#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/47 — 第四輪裁定之**下游重推**（2026-10-29）
=========================================================
觸發：調律者指示「有任何下游影響必須重推」。`docs/46` 之兩項新事實對下游有實質影響：

  I1（E4-R2）BMW `Σ` 列印刷值不可加 ⟹ `⟨O⟩_Σ` 由 1.4980 變 1.5469（+0.0488 = 0.19σ）
       下游：`docs/28`（z vs 1）、`docs/29`（χ²/p）、`docs/30`（12 形式 χ²<3 之通過集）、
             `docs/31`（108 模型掃描 best χ²／單調數／`p_null`）、`docs/33/39` 之收據引用
  I2（E4-R1）`σ_X` 部分對 8 通道**完全共線** ⟹ 任何**線性組合** Σc_iE_i 之 X 誤差為
             `(Σ c_i·coeff_i)·σ_X`（rank-1），非各列獨立
       下游：`docs/18` §四、`docs/42` 之 K/N 一致性（`Δ = E_K − E_np`，`coeff` 反號 ⟹ X 項變 +2σ_X）
  I3（E3）機制 `M1/M2` → `K1/K2`：掃全 corpus 是否仍有機制舊稱之引用

判準（事前）：
  A `⟨O⟩_Σ` 改動後 `⟨O⟩≡1` 仍否證（p < 0.01）⟹ `docs/29` 判詞不變
  B `docs/30` 之通過集若改變，須逐個列出（不得只報數目）
  C `docs/31` 之 `p_null` 仍 ≳ 0.05 ⟹ 掃描結論（不可採信）不變
  D K/N 之 `σ_Δ` 用協方差後不得令結論反轉（若 |z| < 2 則仍「一致」，但 σ 須更正）
"""
import itertools
import json
import math

CHK = []


def chk(label, ok, got, want):
    CHK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(want)})
    print(("PASS " if ok else "FAIL ") + f"{label:<62} | {str(got):<26} | {want}")


def chi2_sf(x, k):
    """上尾機率 P(χ²_k > x)。用 mpmath 正則化上不完全 gamma（精確）。
    ⚠ `scripts/30` 內之自寫 `chi2_sf` 在**奇 dof** 分支把 j-sum 寫成 `range(0)`
    ⟹ 第二項消失，k = 3 時返回 1-dof 值（與 M7 指出 docs/29 之同型錯誤相同）。
    本函數不沿用該實作（見 `docs/47` §四 發現 4）。"""
    import mpmath as mp
    with mp.workdps(30):
        return float(mp.gammainc(mp.mpf(k) / 2, mp.mpf(x) / 2, mp.inf) / mp.gamma(mp.mpf(k) / 2))


print("=" * 112)
print("I1 — BMW Σ 列不自洽：⟨O⟩_Σ 重推（1.4980 → 1.5469）")
print("=" * 112)
X = 2.560
DQED_S, DQCD_S, TOT_S = 0.08, -7.67, -7.84          # BMW Σ 列（MeV）
C_S = 2                                              # Σ⁺−Σ⁻：兩次 u↔d
O_OLD = abs(DQCD_S) / (C_S * X)                      # 1.4980
DQCD_FROM_TOTAL = TOT_S - DQED_S                     # −7.92
O_NEW = abs(DQCD_FROM_TOTAL) / (C_S * X)             # 1.5469
S_O = 0.262                                          # docs/28 之保守 σ（不改）
print(f"  ⟨O⟩_Σ(old, ΔQCD 直給) = |{DQCD_S}|/({C_S}·{X}) = {O_OLD:.4f}")
print(f"  ⟨O⟩_Σ(new, total−ΔQED) = |{DQCD_FROM_TOTAL:.2f}|/({C_S}·{X}) = {O_NEW:.4f}   Δ = {O_NEW-O_OLD:+.4f} = {(O_NEW-O_OLD)/S_O:.2f}σ")
chk("I1a ⟨O⟩_Σ_old = 1.4980（docs/28）", abs(O_OLD - 1.498) < 5e-4, f"{O_OLD:.4f}", "1.4980")
chk("I1b ⟨O⟩_Σ_new = 1.5469", abs(O_NEW - 1.5469) < 5e-4, f"{O_NEW:.4f}", "1.5469")

print("\n" + "=" * 112)
print("I1-P2 — docs/29：χ²(⟨O⟩≡1) 三錨 → p、等效 z")
print("=" * 112)
ANCH_OLD = [(0, 0.8906, 0.107), (1, O_OLD, 0.262), (2, 2.2930, 0.351)]
ANCH_NEW = [(0, 0.8906, 0.107), (1, O_NEW, 0.262), (2, 2.2930, 0.351)]
chi_old = sum(((o - 1) / s) ** 2 for _, o, s in ANCH_OLD)
chi_new = sum(((o - 1) / s) ** 2 for _, o, s in ANCH_NEW)
p_old, p_new = chi2_sf(chi_old, 3), chi2_sf(chi_new, 3)   # chi2_sf = 上尾 P(χ²>x)
print(f"  χ²_old = {chi_old:.4f}（docs/29 報 18.2207）⟹ p = {p_old:.4e}")
print(f"  χ²_new = {chi_new:.4f} ⟹ p = {p_new:.4e}；√χ² 等效 z：{math.sqrt(chi_old):.3f} → {math.sqrt(chi_new):.3f}")
chk("I1-P2a χ²_old ≈ 18.2207（本實現 18.2290，⟨O⟩_N 取位 0.891 vs 0.8906）", abs(chi_old - 18.2207) < 0.02, f"{chi_old:.4f}", "18.2207 ± 0.02")
chk("I1-P2a2 p(χ²_old, 3 dof) = 3.96e−4（M7/scripts33 之正確值；1 dof 為 1.97e−5）", abs(p_old - 3.96e-4)/3.96e-4 < 0.03, f"{p_old:.3e}", "3.96e-4 ± 3%")
chk("I1-P2b χ² 上升（更偏離 ⟨O⟩≡1）", chi_new > chi_old, f"{chi_old:.3f} → {chi_new:.3f}", "上升")
chk("I1-P2c 判詞不變：p 仍 < 0.01（p_old 應 ≈ 3.96e−4）", p_new < 0.01 and p_old < 0.01, f"{p_old:.3e} → {p_new:.3e}", "< 0.01")

print("\n" + "=" * 112)
print("I1-P3 — docs/30：12 形式之 χ²<3 通過集（新舊對照）")
print("=" * 112)
B_, NQ, ND, NG = 4, 3, 88, 3
MODELS = {
    "A1: 1 + n_s·2/N_q": lambda n: 1 + 2 * n / NQ,
    "A2: N_q(n_s+1)/B": lambda n: NQ * (n + 1) / B_,
    "A3: 1 + n_s·1/N_q": lambda n: 1 + n / NQ,
    "A4: 1 + n_s·(B−1)/B": lambda n: 1 + n * (B_ - 1) / B_,
    "A5: 1 + n_s/B": lambda n: 1 + n / B_,
    "A6: 1 + n_s/2": lambda n: 1 + n / 2,
    "A7: (1 + n_s/N_q)²": lambda n: (1 + n / NQ) ** 2,
    "A8: B/(B − n_s)": lambda n: B_ / (B_ - n),
    "A9: (N_q + n_s·NG)/N_q": lambda n: (NQ + n * NG) / NQ,
    "A10: 1 + n_s·(B−2)/(B−1)": lambda n: 1 + n * (B_ - 2) / (B_ - 1),
    "A11: 1 + n_s·ND/(ND+NQ·B)": lambda n: 1 + n * ND / (ND + NQ * B_),
    "A12: 1 + n_s·2/B": lambda n: 1 + 2 * n / B_,
}
res = {}
for name, f in MODELS.items():
    co = sum((((o - f(n))) / s) ** 2 for n, o, s in ANCH_OLD)
    cn = sum((((o - f(n))) / s) ** 2 for n, o, s in ANCH_NEW)
    res[name] = {"chi2_old": round(co, 4), "chi2_new": round(cn, 4),
                 "p3_old_correct": round(chi2_sf(co, 3), 6), "p3_new_correct": round(chi2_sf(cn, 3), 6),
                 "pass_old": co < 3.0, "pass_new": cn < 3.0}
    flag = ("✅→✅" if res[name]["pass_old"] and res[name]["pass_new"] else
            ("✅→🔴" if res[name]["pass_old"] else ("🔴→✅" if res[name]["pass_new"] else "🔴→🔴")))
    print(f"  {name:<28} χ² {co:8.3f} → {cn:8.3f}   p(3dof) {chi2_sf(co,3):8.4f} → {chi2_sf(cn,3):8.4f}   {flag}")
p_old_set = [k for k, v in res.items() if v["pass_old"]]
p_new_set = [k for k, v in res.items() if v["pass_new"]]
print(f"\n  通過集：{len(p_old_set)}/12 → {len(p_new_set)}/12")
print(f"  退出：{[k for k in p_old_set if k not in p_new_set] or '無'}")
print(f"  新入：{[k for k in p_new_set if k not in p_old_set] or '無'}")
chk("I1-P3a 通過集 = 7/12（docs/30 所報）", len(p_old_set) == 7, f"{len(p_old_set)}/12", "7/12")
chk("I1-P3b 新錨下通過集不變（逐個比對）", set(p_old_set) == set(p_new_set),
    f"{len(p_new_set)}/12 同集" if set(p_old_set) == set(p_new_set) else str(p_new_set), "同集")

print("\n" + "=" * 112)
print("I1-P4 — docs/31：108 模型掃描（新錨）之 best χ²／單調數／p_null")
print("=" * 112)
B4, LD = 4.0, 12
ENV = {"N": ("u", "d"), "Σ": ("u", "s"), "Ξ": ("s", "s")}
KAP = (0.5, 1, 2, 3, 5, 10, 30, 100, 300)


def rho(dep, kappa, env, res_):
    d = list(range(LD + 1))
    V = [0.0] * (LD + 1)
    for pp in env:
        for i in d:
            V[i] -= kappa * B4 ** (-abs(i - dep[pp]))
    # 三對角矩陣：用 numpy 較穩
    import numpy as np
    H = np.diag(V) + np.diag([-math.sqrt(B4 - 1)] * LD, 1) + np.diag([-math.sqrt(B4 - 1)] * LD, -1)
    _, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    nodes = [1.0] + [B4 * (B4 - 1) ** (k - 1) for k in range(1, LD + 1)]
    v = v / math.sqrt(sum(n * x * x for n, x in zip(nodes, v)))
    return [x * x for x in v] if res_ == "site" else [n * x * x for n, x in zip(nodes, v)]


def scan(anch):
    a = {"N": anch[0][1:], "Σ": anch[1][1:], "Ξ": anch[2][1:]}
    ms = []
    for perm in itertools.permutations(range(3)):
        dep = {"u": perm[0], "d": perm[1], "s": perm[2]}
        for rs in ("site", "shell"):
            for k in KAP:
                tr = tuple(float(3 * rho(dep, k, ENV[name], rs)[dep["d"]]) for name in ("N", "Σ", "Ξ"))
                c = sum((tr[i] - a[kk][0]) ** 2 / a[kk][1] ** 2 for i, kk in enumerate(("N", "Σ", "Ξ")))
                ms.append({"triple": tr, "chi2": c})
    ms.sort(key=lambda m: m["chi2"])
    return ms


def pnull(ms, best_chi, nmc=20000, seed=4729):
    import numpy as np
    rng = np.random.default_rng(seed)
    sig = [ANCH_OLD[0][2], ANCH_OLD[1][2], ANCH_OLD[2][2]]
    grid = np.array([m["triple"] for m in ms])
    mins = []
    for _ in range(nmc):
        fake = np.exp(rng.uniform(math.log(0.2), math.log(4.0), 3))
        mins.append((((grid - fake[None, :]) / np.array(sig)[None, :]) ** 2).sum(axis=1).min())
    return float((np.array(mins) <= best_chi).mean())


ms_old, ms_new = scan(ANCH_OLD), scan(ANCH_NEW)
mono_old = sum(1 for m in ms_old if m["triple"][0] < m["triple"][1] < m["triple"][2])
mono_new = sum(1 for m in ms_new if m["triple"][0] < m["triple"][1] < m["triple"][2])
pn_old, pn_new = pnull(ms_old, ms_old[0]["chi2"]), pnull(ms_new, ms_new[0]["chi2"])
print(f"  best χ²：{ms_old[0]['chi2']:.2f}（docs/31 報 26.68）→ {ms_new[0]['chi2']:.2f}")
print(f"  單調升模型：{mono_old}/108（docs/31 報 12）→ {mono_new}/108")
print(f"  p_null：{pn_old:.3f} → {pn_new:.3f}（2×10⁴ 次）")
chk("I1-P4a best χ²_old ≈ 26.68（本實現 26.67，取位 0.01）", abs(ms_old[0]["chi2"] - 26.68) < 0.02, f"{ms_old[0]['chi2']:.2f}", "26.68 ± 0.02")
chk("I1-P4b 單調數_old = 12/108", mono_old == 12, f"{mono_old}/108", "12/108")
chk("I1-P4c 新錨下 p_null 仍 ≳ 0.05（掃描結論不變）", pn_new > 0.05, f"{pn_new:.3f}", "> 0.05")

print("\n" + "=" * 112)
print("I2 — 共線 σ_X 對線性組合之影響（docs/42 之 K/N 一致性）")
print("=" * 112)
# data/42：E_old/E_new、s_old/s_new、coeff
d42 = json.load(open("data/42_flag2024_x_revision.json", encoding="utf-8"))["data"]["channels"]
E_K, E_N = d42["K+−K0"], d42["n−p"]
SX2 = 0.1159
S_DM = {"K": 0.021, "N": 0.000001}
for tag, key_e, key_s in (("old（σ_X = 0.09434）", "E_old", "s_old"), ("new（σ_X = 0.1159）", "E_new", "s_new")):
    d = E_K[key_e] - E_N[key_e]
    s_indep = math.hypot(E_K[key_s], E_N[key_s])              # docs/42 之算法（當獨立）
    # 協方差：coeff_K = −1, coeff_N = +1 ⟹ ∂Δ/∂X = +1 − (−1) = +2
    var_dm = S_DM["K"] ** 2 + S_DM["N"] ** 2
    s_cov = math.sqrt(var_dm + (2 * (SX2 if "new" in tag else 0.09434)) ** 2)
    print(f"  {tag}: Δ = {d:+.4f} MeV")
    print(f"     獨立算法 σ_Δ = {s_indep:.4f} ⟹ z = {abs(d)/s_indep:.2f}σ（docs/42 現載）")
    print(f"     協方差 σ_Δ = {s_cov:.4f} ⟹ z = {abs(d)/s_cov:.2f}σ  （X 項 = |coeff_K − coeff_N|·σ_X = 2σ_X）")
    if "new" in tag:
        s_new_cov, z_new_cov, s_new_ind, z_new_ind = s_cov, abs(d) / s_cov, s_indep, abs(d) / s_indep
chk("I2a docs/42 之 1.21σ 可重現（獨立算法）", abs(z_new_ind - 1.21) < 0.02, f"{z_new_ind:.2f}σ", "1.21σ")
chk("I2b 協方差修正後 σ_Δ = 0.233（非 0.164）", abs(s_new_cov - 0.2328) < 1e-3, f"{s_new_cov:.4f}", "0.2328")
chk("I2c 結論不反轉（協方差 z 仍 < 2σ）", z_new_cov < 2, f"{z_new_cov:.2f}σ", "< 2σ")

print("\n" + "=" * 112)
print("I3 — 機制 M1/M2 → K1/K2 之殘留引用掃描")
print("=" * 112)
import pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
hits = {}
for f in list((ROOT / "docs").glob("*.md")) + list((ROOT / "scripts").glob("*.py")) + [ROOT / "01_problem_record.md", ROOT / "final_result.md", ROOT / "README.md"]:
    if not f.exists():
        continue
    txt = f.read_text(encoding="utf-8", errors="ignore")
    n = txt.count("機制 M1") + txt.count("機制 M2") + txt.count("機制M1") + txt.count("機制M2")
    if n:
        hits[f.name] = n
print("  " + ("；".join(f"{k}: {v}" for k, v in sorted(hits.items())) or "零命中"))
ALLOWED = {"15_exponent8.md", "15_exponent8.py", "01_problem_record.md",
           "46_handoff_adjudication_round4.md", "46_handoff_adjudication_round4.py",
           "47_downstream_rerun_round4.py"}
chk("I3a 舊稱只餘『原文保留 ＋ 改名記錄』檔（無下游以舊稱為 live 標籤）",
    set(hits) <= ALLOWED, sorted(hits), "⊆ 允許集 6 檔")

DATA = {
    "I1_sigma_anchor": {"O_old": round(O_OLD, 4), "O_new": round(O_NEW, 4), "delta_sigma": round((O_NEW - O_OLD) / S_O, 3)},
    "I1_chi2_docs29": {"chi2_old": round(chi_old, 4), "chi2_new": round(chi_new, 4),
                       "p_old": p_old, "p_new": p_new},
    "I1_docs30": {"n_pass_old": len(p_old_set), "n_pass_new": len(p_new_set),
                  "dropped": [k for k in p_old_set if k not in p_new_set],
                  "added": [k for k in p_new_set if k not in p_old_set], "table": res},
    "I1_docs31": {"chi2_old": round(ms_old[0]["chi2"], 4), "chi2_new": round(ms_new[0]["chi2"], 4),
                  "mono_old": mono_old, "mono_new": mono_new, "pnull_old": pn_old, "pnull_new": pn_new},
    "I2_KN": {"s_indep": round(s_new_ind, 4), "z_indep": round(z_new_ind, 3),
              "s_cov": round(s_new_cov, 4), "z_cov": round(z_new_cov, 3)},
    "I3_rename_hits": hits,
    "n_pass": sum(1 for c in CHK if c["ok"]), "n_total": len(CHK),
}
print("\n" + "=" * 112)
print(f"RESULT: {DATA['n_pass']}/{DATA['n_total']} PASS")
print("=" * 112)
with open("data/47_downstream_rerun_round4.json", "w", encoding="utf-8") as f:
    json.dump({"checks": CHK, "data": DATA}, f, ensure_ascii=False, indent=1)
print("→ data/47_downstream_rerun_round4.json")