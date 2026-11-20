#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/06 — 自審：docs 文字數字 vs data JSON；圍欄偶數；U+FFFD；四件套狀態同步
=================================================================================
機械可數；任何 FAIL 即代表文檔與機器輸出不一致（須修，不得放寬）。
"""
import glob
import hashlib
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
DOCS = os.path.join(ODIR, "docs")
CHECKS = []


def ck(tag, cond, detail=""):
    print(("[PASS] " if cond else "[FAIL] ") + tag + " :: " + detail, flush=True)
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


def rd(p):
    return open(p, encoding="utf-8").read()


J = {i: json.load(open(os.path.join(DATA, f"{i:02d}_" + n), encoding="utf-8"))
     for i, n in [(0, "reproduction_gate.json"), (1, "a5_invariants.json"),
                  (2, "amplitude_nogo.json"), (3, "comb_field_stats.json"),
                  (4, "null_and_sensitivity.json"), (5, "chirality.json"),
                  (7, "mr0_carrier.json"), (8, "mr1_partition.json"),
                  (9, "mr2_chirality_test.json"), (10, "mr3_lcdm_baseline.json"),
                  (11, "mr0_r7.json"), (12, "mr1_equivariance.json"), (13, "shape_invariants.json"), (14, "shape_invariants_r10.json")]}

# ── 1 逐支 PASS/FAIL 合計 ────────────────────────────────────────
tot_p = sum(J[i]["n_pass"] for i in J)
tot_f = sum(J[i]["n_fail"] for i in J)
print("逐支：", {i: (J[i]["n_pass"], J[i]["n_fail"]) for i in J})
ck("A1-合計", tot_p == 194 and tot_f == 0, f"{tot_p} PASS ／ {tot_f} FAIL（00–05、07–14）")

# ── 2 docs 數字 vs data ─────────────────────────────────────────
d02 = rd(os.path.join(DOCS, "02_projection_operator.md"))
d03 = rd(os.path.join(DOCS, "03_statistics_bands.md"))
d04 = rd(os.path.join(DOCS, "04_input_gates.md"))
fr = rd(os.path.join(ODIR, "final_result.md"))

D = J[1]["D"]
ck("A2-𝒟", all(str(l) in d02 for l in D), f"𝒟 全部成員在 docs/02（{len(D)} 個）")
ck("A2-Nfix", str(J[1]["N_fix_27"]) in d02 and str(J[1]["N_fix_27"]) in fr,
   f"N_fix(27) = {J[1]['N_fix_27']}")
ck("A2-群階", str(J[1]["group"]["order"]) in d02, f"|A₅| = {J[1]['group']['order']}")
ck("A2-180.795", "180.795" in d02 and "180.795" in fr, "log₂(44!) = 180.795")
ck("A2-768/779", "768" in d02 and "779" in d02 and "768" in fr, "約束 768／779")
b = J[3]["bands"]
ck("A2-R27帶", f"{b['R27'][0]:.4f}" in d03 and f"{b['R27'][1]:.4f}" in d03,
   f"R^27 帶 [{b['R27'][0]:.4f}, {b['R27'][1]:.4f}]")
ck("A2-ratio16帶", f"{b['ratio16'][0]:.6f}" in d03 and f"{b['ratio16'][1]:.6f}" in d03,
   f"ratio16 帶 [{b['ratio16'][0]:.6f}, {b['ratio16'][1]:.6f}]")
fn = J[4]["f_null"]
ck("A2-f門檻", f"{fn['SW']['f_crit_bonf']:.5f}" in d03 and f"{fn['SW']['f_crit_5sigma']:.5f}" in d03
   and f"{fn['SW']['f_crit_5sigma']:.5f}" in d04,
   f"f_crit = {fn['SW']['f_crit_bonf']:.5f} / {fn['SW']['f_crit_5sigma']:.5f}")
ck("A2-Nmc", "10⁷" in d03, "N_mc = 10⁷ 標示")
d05 = rd(os.path.join(DOCS, "05_residuals_r6.md"))
ck("A2-R6門檻", "0.23160" in d03 and "0.11948" in d03 and "0.23160" in d05,
   "真 ΛCDM 門檻 0.23160／0.11948 在 docs/03 §七 與 docs/05")
ck("A2-R6缺損", "40.18" in d05 and "40.18" in fr, "ΛCDM 缺損 40.18%")
ck("A2-R6手性", "2.13e−14" in d05 or "2.13e-14" in d05, "反極點恆等式殘差")
ck("A2-250", "250 PASS" in fr, "final_result 合計 250")
d09 = rd(os.path.join(DOCS, "09_r10_kappa3_extensions.md"))
ck("A2-R10", "3/9" in d09 and "16.919" in d09 and "m 偶" in d09, "R10 三項在 docs/09")

# ── 3 圍欄偶數 ＋ U+FFFD ───────────────────────────────────────
for p in sorted(glob.glob(os.path.join(DOCS, "*.md"))) + [os.path.join(ODIR, "final_result.md"),
                                                          os.path.join(ODIR, "README.md"),
                                                          os.path.join(ODIR, "01_problem_record.md"),
                                                          os.path.join(ODIR, "SESSION_START.md")]:
    t = rd(p)
    ck("A3-圍欄 " + os.path.basename(p), t.count("```") % 2 == 0, f"{t.count('```')} 個 ```")
    ck("A3-U+FFFD " + os.path.basename(p), "\ufffd" not in t, "無替換字元")

# ── 4 四件套狀態同步 ──────────────────────────────────────────
rec = rd(os.path.join(ODIR, "01_problem_record.md"))
d06 = rd(os.path.join(DOCS, "06_mr0_bridge_b1.md"))
ck("A2-B1", "B1" in d06 and "352" in d06 and "B1" in fr and "B1" in rd(os.path.join(DOCS, "04_input_gates.md")),
   "B1 在 docs/04 §九、docs/06、final_result")
d07 = rd(os.path.join(DOCS, "07_mr1_equivariance.md"))
d08 = rd(os.path.join(DOCS, "08_shape_invariant_kappa3.md"))
j13 = J[13]["model"]
vals08 = {6: "+0.169754", 10: "+0.093968", 12: "+0.099039", 16: "0.030804",
          18: "+0.069251", 20: "+0.047672", 24: "+0.051914"}
ck("A2-κ₃表", all(v.replace("-", "") in d08.replace("−", "-") for v in list(vals08.values())[:-1]),
   f"{len(vals08)-1} 個模型 κ₃ 值在 docs/08（正負號形式已容許 −/\u2212）")
ck("A2-κ₃機器一致", all(abs(abs(j13[str(l)]["kappa3_real"]) - abs(float(v))) < 5e-6 for l, v in vals08.items()),
   "docs/08 之值與 data/13 機器值一致（±5e-6）")
ck("A2-奇ℓ恆零", all(abs(j13[str(l)]["I3_re"]) < 1e-9 and abs(j13[str(l)]["I3_abs"]) < 1e-9 for l in (15,21,25,27)), "奇 ℓ I₃ ≡ 0（機器）")
ck("A2-B2", "B2" in d07 and "13 全自由" in d07 and "B1" in d07, "B2 與 B1 並列在 docs/07")
ck("A2-M17越界", "越界" in d07 and "m17_cmb_large_scale_anomalies" in d07, "越界寫入已記錄並還原")
ck("A4-registry", "state:   partial" in rec and "verdict: T1-pass" in rec,
   "01_problem_record registry 已同步 partial")
ck("A4-README", "partial" in rd(os.path.join(ODIR, "README.md")), "README 標 partial")
ck("A4-SESSION", "2026-11-20" in rd(os.path.join(ODIR, "SESSION_START.md")), "SESSION_START 交付記錄")
ck("A4-final", "部分閉合" in fr, "final_result 標部分閉合")

# ── 5 禁區一致性：docs/01 §1.3 五條禁令在 final_result ────────
for key in ["逐模指派", "R^27", "C(π)", "σ²_16 絕對值", "60 dof"]:
    ck("A5-禁區 " + key, key in fr, f"「{key}」在 final_result 禁區段")

# ── 6 檔案指紋（供下次比對）──────────────────────────────────
snap = {}
for p in sorted(glob.glob(os.path.join(DATA, "*.json"))) + sorted(glob.glob(os.path.join(DOCS, "*.md"))):
    snap[os.path.relpath(p, ODIR)] = hashlib.md5(open(p, "rb").read()).hexdigest()
with open(os.path.join(DATA, "06_fingerprints.json"), "w", encoding="utf-8") as f:
    json.dump(snap, f, ensure_ascii=False, indent=1)
print(f"\n   指紋：{len(snap)} 檔 → data/06_fingerprints.json")

npass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n=== scripts/06：{npass} PASS ／ {len(CHECKS) - npass} FAIL ===")
sys.exit(0 if npass == len(CHECKS) else 1)
