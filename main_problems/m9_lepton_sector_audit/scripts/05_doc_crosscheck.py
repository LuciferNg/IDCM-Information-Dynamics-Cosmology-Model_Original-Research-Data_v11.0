#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""M9 scripts/05 — 文檔 × 數據 交叉核對（data ↔ docs 一致性）

零重算物理：只讀 data/*.json 與 docs/*.md，核對 headline 數字逐項對得上。
任何 FAIL = docs 與腳本輸出不一致 ⟹ 必須修正 docs（唔可以改數字凑）。
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
DOCS = os.path.join(ROOT, "docs")

checks = []


def load(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return json.load(f)


def doc(name):
    with open(os.path.join(DOCS, name), encoding="utf-8") as f:
        return f.read()


def chk(label, cond, detail=""):
    checks.append({"check": label, "pass": bool(cond), "detail": detail})
    print(("PASS " if cond else "FAIL ") + label + ("  " + detail if detail else ""))


d1 = load("01_p1_p2_sigma.json")
d2 = load("02_phi_lock.json")
d3 = load("03_language.json")
d4 = load("04_spectral.json")
m1, m2, m3, m4 = (doc("01_sigma_table.md"), doc("02_phi_lock_gate.md"),
                  doc("03_language_consistency.md"), doc("04_spectral_uniqueness.md"))

print("=" * 70)
print("M9 scripts/05 — 文檔 × 數據 交叉核對")
print("=" * 70)

# ── 01 ───────────────────────────────────────────────────────────
claims = {c["label"]: c for c in d1["claims"]}
c1 = float(claims["C1 phi^(11+7/88)"]["sigma"])
c3 = float(claims["C3 integer 207"]["sigma"])
rel_u = float(d1["sigma_criterion"]["rel_u"])
chk("01-a C1 σ = 428.4（json）", abs(c1 - 428.3718895) < 1e-6, f"{c1:.6f}")
chk("01-b docs/01 記 428.4σ", "428.4σ" in m1 and "428.4" in m1)
chk("01-c C3 σ = 50 373（json）", abs(c3 - 50373.32609) < 1e-4, f"{c3:.4f}")
chk("01-d docs/01 記 50 373σ", "50 373" in m1)
chk("01-e u_rel = 2.225e−8（json）", abs(rel_u - 2.225e-8) < 1e-11, f"{rel_u:.4e}")
chk("01-f docs/01 記 2.225e−8", "2.225e−8" in m1)
chk("01-g C2 exact-def σ = 0（json）",
    float(claims["C2' (g7/g1)(F_mu/F_e) exact-def"]["sigma"]) == 0.0)
chk("01-h docs/01 記 tautology", "恆等" in m1 or "tautology" in m1)

# ── 02 ───────────────────────────────────────────────────────────
chk("02-a 三個指數全部 ∉ K_φ（json）",
    all(not r["in_K_phi"] for r in d2["membership"]))
chk("02-b docs/02 記 ∉ K_φ", "∉ K_φ" in m2)
mc528 = {r["tau"]: r["p_mc"] for r in d2["mc"] if r["D"] == 528}
chk("02-c 528 格 P(1e−5) = 0.0219（json）", abs(mc528[1e-5] - 0.021907) < 5e-6,
    f"{mc528[1e-5]:.6f}")
chk("02-d docs/02 記 0.0219", "0.0219" in m2)
chk("02-e 528 半格 = 0.04558 %（json density）",
    abs([r for r in d2["density"] if r["D"] == 528][0]["half_cell_rel"] - 4.558e-4) < 1e-6)
chk("02-f docs/02 記 0.0456 %", "0.0456" in m2)
chk("02-g 上游 sha256 一致（json ↔ docs/02）",
    d2["upstream_audit"]["sha256"] in m2 and
    d2["upstream_audit"]["sha256"].startswith("585d6a64"))
chk("02-h 記載 0.0003 % 由 m_e=0.511 重現（json）",
    abs(float(d2["upstream_audit"]["repro_k_pct_with_0511"]) - 0.000283917) < 1e-9)
chk("02-i 記載 0.0015 % 由 m_e=0.511 重現（json）",
    abs(float(d2["upstream_audit"]["repro_val_pct_with_0511"]) - 0.00151374) < 1e-8)
chk("02-j docs/02 記 0.0015137", "0.0015137" in m2)
chk("02-k 真值 0.000953 %（json）",
    abs(float(d2["upstream_audit"]["true_val_pct"]) - 0.000953004) < 1e-9)
chk("02-l docs/02 記 0.000953", "0.000953" in m2)
chk("02-m 殘差 1.05 % 格（json）",
    any(abs(float(r["resid_cell"]) - 0.01046) < 1e-4 for r in d2["grid_residual"]
        if r["D"] == 528 and r["candidate"].startswith("C1")))
chk("02-n docs/02 記 1.05 %", "1.05 %" in m2 or "1.05" in m2)

# ── 03 ───────────────────────────────────────────────────────────
chk("03-a 加法恆等 42+457=499（json）",
    d3["additivity_b"]["mu/e"] + d3["additivity_b"]["tau/mu"] == d3["additivity_b"]["tau/e"]
    == 499)
chk("03-b docs/03 記 42/457/499", "42" in m3 and "457" in m3 and "499" in m3)
chk("03-c CAP = 784.9646644（json）", abs(float(d3["cap"]) - 784.9646644) < 1e-6,
    d3["cap"])
chk("03-d docs/03 記 CAP", "784.9647" in m3 or "784.9646644" in m3)
chk("03-e τ/e = 4.43× CAP（json best dev 77.4 %）",
    abs(float(d3["l0_best"]["tau/e"]["dev_pct"]) - 77.4256) < 1e-3)
chk("03-f docs/03 記 4.43× CAP", "4.43" in m3)
chk("03-g 528 半格粗 μ/e 精度 >1e4×", "20 490" in m3 or "2.049e+04" in m3)

# ── 04 ───────────────────────────────────────────────────────────
l0 = {r["target"]: r for r in d4["L0"]}
chk("04-a L0 μ/e dev = −1.1524 %（json）",
    abs(float(l0["mu/e"]["dev_pct"]) + 1.152358) < 1e-5)
chk("04-b docs/04 記 1.1524 %", "1.1524" in m4)
lad = {r["level"]: r for r in d4["dof_ladder"]}
chk("04-c L2 P_MC(1e−5) = 0.524（json）", abs(lad["L2"]["p_mc_1e-05"] - 0.524) < 5e-3,
    f"{lad['L2']['p_mc_1e-05']:.4f}")
chk("04-d docs/04 記 52 %", "52 %" in m4 or "0.524" in m4)
chk("04-e 位元預算全負（json）", all(b["bits"] < 0 for b in d4["bit_budget"]))
chk("04-f docs/04 記負位元", "−14.7" in m4 or "-14.7" in m4)
ho = {r["target"]: r for r in d4["holdout"]}
chk("04-g holdout p/e dev = −57.2495 %（json）",
    abs(float(ho["p/e"]["L0"]["dev_pct"]) + 57.249488) < 1e-5)
chk("04-h docs/04 記 −57.25 %", "57.25" in m4)
chk("04-i O209 讀值差 1314×（json）",
    abs(float(d4["assignment_reading"][0]["factor_off"]) - 1314.0) < 1.0)
chk("04-j docs/04 記 1 314×", "1 314" in m4)

# ── 交叉檔 ───────────────────────────────────────────────────────
chk("X-a 四份 JSON pass == total（無 FAIL）",
    all(d["pass"] == d["total"] for d in (d1, d2, d3, d4)),
    "; ".join(f"{d['pass']}/{d['total']}" for d in (d1, d2, d3, d4)))
fr = open(os.path.join(ROOT, "final_result.md"), encoding="utf-8").read()
chk("X-b final_result 記 29/29",
    bool(re.search(r"29/29", fr)) and (d1["pass"] + d2["pass"] + d3["pass"] + d4["pass"]) == 29)
chk("X-c final_result 記 428.4σ", "428.4σ" in fr)
chk("X-d final_result 記 CAP 784.9646644", "784.9646644" in fr)

npass = sum(1 for c in checks if c["pass"])
print("-" * 70)
print(f"PASS {npass}/{len(checks)}")
out = os.path.join(DATA, "05_doc_crosscheck.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump({"checks": checks, "pass": npass, "total": len(checks)},
              f, ensure_ascii=False, indent=2)
print("->", out)
raise SystemExit(0 if npass == len(checks) else 1)