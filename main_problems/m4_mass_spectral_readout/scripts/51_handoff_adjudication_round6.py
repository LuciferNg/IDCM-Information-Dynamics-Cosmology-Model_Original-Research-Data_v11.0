#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/51 — 第六輪入站裁定之計算支撐（2026-11-23）
=====================================================
入站：M1 2026-11-23 03:40（`docs/13_residual_push.md`）兩項影響：

  ① `coverage(k)` 延伸至 `k = 8`（uniform-in-t, `ε = 1e−3`）：
     `0.1968 / 5.5353 / 12.269 / 17.0396 / 21.4821 / 25.8069 / 30.0554 / 34.2289 %`
     ★`coverage(8, uniform) = 34.229% ≈ coverage(4, log-uniform 池內) = 35.155%`
       ⟹「加支撐」與「換測度」在數值上等價 ⟹ 引用必須同時報 `(k, 測度)`
  ② `s` 之 `ρ(k=3) = 7` 降為**條件性**（質量 1σ `σ_rel = 8.6%`，`N = 1e4`，74% 抽樣保有編碼；
     `μ` 為 `P = 1.000` 穩定）⟹ 任何「`μ1→s7→τ513→c4285`」梯度須改寫為
     「`μ`（穩定唯一）→ `s`（條件性）→ `τ513` → `c4285`」

判準（事前）：
  A M1 之 k=1..8 序列嚴格單調（算術核；數值屬 M1 產權，只核內部一致性）
  B `34.229 < 35.155` 且比值 ≈ 0.974（等價性陳述之數值前提）
  C M4 docs 內**零**「μ1→s7→τ513→c4285」梯度引用 ⟹ ② 對本 O 下游影響 = 0
  D M4 自身之測度效應方向一致（log-uniform > uniform，沿用 `docs/48` 之 38.11× 量級）
"""
import glob
import json
import os

CHK = []
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)


def chk(label, ok, got, want):
    CHK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(want)})
    print(("PASS " if ok else "FAIL ") + f"{label:<62} | {str(got):<30} | {want}")


print("=" * 104)
print("A — M1 coverage(k) 序列（k = 1..8, uniform-in-t, ε = 1e−3）")
print("=" * 104)
COV = [0.1968, 5.5353, 12.269, 17.0396, 21.4821, 25.8069, 30.0554, 34.2289]
print("  " + "  ".join(f"k={i+1}: {v:.4f}%" for i, v in enumerate(COV)))
mono = all(COV[i] < COV[i + 1] for i in range(len(COV) - 1))
chk("A1 序列嚴格單調（8 值）", mono, f"{len(COV)} 值", "嚴格單調")
chk("A2 支撐閘 coverage ≤ 0.1 ⟹ 只到 k = 2（k=3 為 12.269%）", COV[1] / 100 <= 0.1 < COV[2] / 100,
    f"k2={COV[1]}% / k3={COV[2]}%", "k≤2")

print("\n" + "=" * 104)
print("B — 「加支撐」≈「換測度」之數值前提")
print("=" * 104)
C8_U = COV[7]
C4_LOG = 35.155
print(f"  coverage(k=8, uniform) = {C8_U:.3f}% ｜ coverage(k=4, log-uniform 池內) = {C4_LOG:.3f}%")
print(f"  比值 = {C8_U/C4_LOG:.4f}（差 {abs(C8_U-C4_LOG):.3f}pp）")
chk("B1 34.229 < 35.155（等價性成立之方向）", C8_U < C4_LOG, f"{C8_U} < {C4_LOG}", "<")
chk("B2 兩者差 < 1pp ⟹ 可稱『數值等價』", abs(C8_U - C4_LOG) < 1.0, f"{abs(C8_U-C4_LOG):.3f}pp", "<1pp")

print("\n" + "=" * 104)
print("C — M4 docs 有無「μ1→s7→τ513→c4285」梯度引用")
print("=" * 104)
PAT = ["τ513", "c4285", "μ1→s7", "μ₁→s₇"]
hits = []
for p in sorted(glob.glob("docs/*.md")) + ["01_problem_record.md", "final_result.md"]:
    if not os.path.exists(p):
        continue
    t = open(p, encoding="utf-8", errors="ignore").read()
    for k in PAT:
        if k in t:
            hits.append((p, k))
print(f"  命中：{hits or '零'}")
chk("C1 M4 零引用該梯度 ⟹ ② 對本 O 下游影響 = 0", not hits, hits or "0", "0")
# 反向核：M4 是否用「s 之 ρ=7」作結構證據
rho7 = []
for p in sorted(glob.glob("docs/*.md")):
    t = open(p, encoding="utf-8", errors="ignore").read()
    if ("ρ(k=3)" in t or "ρ(3)" in t) and "7" in t:
        rho7.append(os.path.basename(p))
print(f"  用 ρ(k=3)=7 之 docs：{rho7 or '零'}")
chk("C2 M4 未以 ρ(k=3)=7 作結構證據（若有須改『條件性』）", not rho7, rho7 or "0", "0")

print("\n" + "=" * 104)
print("D — M4 自身測度效應（沿用 docs/48：同一域、同一池）")
print("=" * 104)
d48 = json.load(open("data/48_handoff_adjudication_round5.json", encoding="utf-8"))["data"]["H1_convention"]
print(f"  M4 F1 池：log-uniform {d48['coverage_log_uniform_pct']}% vs uniform-in-t {d48['coverage_uniform_t_pct']}%"
      f"（{d48['ratio']}×）")
chk("D1 M4 側方向一致（log-uniform > uniform）", d48["coverage_log_uniform_pct"] > d48["coverage_uniform_t_pct"],
    f"{d48['ratio']}×", ">1×")

DATA = {
    "coverage_k1_8_uniform": COV,
    "gate": {"k_max_under_0.1": 2, "k3_pct": COV[2]},
    "equivalence": {"c8_uniform": C8_U, "c4_log_inner": C4_LOG, "diff_pp": round(abs(C8_U - C4_LOG), 3)},
    "m4_gradient_citations": hits, "m4_rho3_evidence_docs": rho7,
    "m4_own_measure_factor": d48["ratio"],
    "n_pass": sum(1 for c in CHK if c["ok"]), "n_total": len(CHK),
}
print("\n" + "=" * 104)
print(f"RESULT: {DATA['n_pass']}/{DATA['n_total']} PASS")
print("=" * 104)
with open("data/51_handoff_adjudication_round6.json", "w", encoding="utf-8") as f:
    json.dump({"checks": CHK, "data": DATA}, f, ensure_ascii=False, indent=1)
print("→ data/51_handoff_adjudication_round6.json")