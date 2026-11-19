#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M22 scripts/00 — 重現閘（SESSION_START §二 ①–④）
=================================================
硬性：數唔對 ⟹ 停。

① 核 M17 `01_problem_record.md` §3 五個輸入門現狀（CQ1 可做、CQ2–CQ4 落空、CQ5 no-go）
② 核 M6 `final_result.md` 之 `dof_scale = 1` 定理陳述（唔重推）
③ 核 v10 P4 `docs/73` 之 60 dof 硬牆數字（節點 0.2348 σ/μK、梳狀 −0.0032 σ/μK、比 −0.01×）
④ 唔准重跑 M17／M6／M16／P4 腳本（本檔零 import 上游、零 exec；只讀檔 ＋ 自行重算常數）

另加：本 M 下游要用嘅五條結構數字，零 import 重算（唔係重推定理，係核對輸入）：
  γ_n = 4 sin²(πn/88)、γ_n + γ_{44−n} = 4、Σγ = 90、(Σd², Σd, k=5) ⟹ {1,3,3,4,5} 唯一、
  Σd(d−1) = 44 = N_d、Σd = 16 = B²、|A₅| = 60 = N_d + B²

精確度：環內恆等式用解析恆等式 ＋ dps40 數值雙核（本 M 唔重開 ℚ(ζ₈₈) 環核，M6 已做）。
"""
import json
import os
import re
import sys
from fractions import Fraction
from itertools import combinations_with_replacement

import mpmath as mp

mp.mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")
V11 = os.path.dirname(os.path.dirname(ODIR))          # IDCM-v11.0--Reality-Architecture
ROOT = os.path.dirname(V11)                           # IDCM-Original
M17 = os.path.join(V11, "main_problems", "m17_cmb_large_scale_anomalies", "01_problem_record.md")
M6 = os.path.join(V11, "main_problems", "m6_scale_dof_theorem", "final_result.md")
P4_73 = os.path.join(ROOT, "IDCM-v10.0--It-From-Beat", "main_problems", "p4_cmb_cold_spot",
                     "docs", "73_comb_test_and_calibration.md")
P4_72 = os.path.join(ROOT, "IDCM-v10.0--It-From-Beat", "main_problems", "p4_cmb_cold_spot",
                     "docs", "72_o2_handoff_adjudication.md")

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


def rd(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


# ══════════════════════════════════════════════════════════════════
# ① M17 §3 五個輸入門現狀
# ══════════════════════════════════════════════════════════════════
print("\n① M17 §3 輸入門現狀（只讀檔，唔重跑 M17）")
t17 = rd(M17)
sec3 = t17.split("## §3")[1].split("## §4")[0]
for code, pat, want in [
    ("CQ1", r"\| \*\*CQ1\*\* \|.*⚪ 可做", "⚪ 可做"),
    ("CQ2", r"\| \*\*CQ2\*\* \|.*⚪ 待輸入門", "⚪ 待輸入門"),
    ("CQ3", r"\| \*\*CQ3\*\* \|.*⚪ 待輸入門", "⚪ 待輸入門"),
    ("CQ4", r"\| \*\*CQ4\*\* \|.*⚪ 待輸入門", "⚪ 待輸入門"),
    ("CQ5", r"\| \*\*CQ5\*\* \|.*🔴 \*\*結構 no-go\*\*", "🔴 結構 no-go"),
]:
    ck(f"G1-{code}", re.search(pat, sec3) is not None, f"§3 表列 = {want}")
ck("G1-門硬性", "輸入門（本版新增，硬性）" in sec3, "「每項開始前必須先書面寫出」")
ck("G1-禁反推", "用觀測反推常數" in sec3, "禁由觀測反推 v11 常數")

# ══════════════════════════════════════════════════════════════════
# ② M6 dof = 1 定理陳述
# ══════════════════════════════════════════════════════════════════
print("\n② M6 `dof_scale = 1` 定理陳述（只讀檔，唔重推）")
t6 = rd(M6)
ck("G2-標題", "`dof_scale = 1` **恰**" in t6, "恰一（上下界夾死）")
ck("G2-上界", "7 個 Galois 軌道之標度被迫全等" in t6, "上界 = Galois 協變")
ck("G2-下界", "`dof = 0` 在框架內不可能" in t6, "下界 = degree-1 齊次")
ck("G2-B1", "ℚ-多項式結構常數 ⟹ 零連續自由度" in t6, "代數 sector 零連續 dof（定理）")
ck("G2-B2", "尺度軸恰一 ⟺ 全部標度比可預言" in t6, "B2：d=1 全部標度比可預言")
ck("G2-B3", "帳本 ≤9" in t6, "B3：自由數只在超越槽")
ck("G2-錨", "M₀/M_Pl = 2.2852931e−21" in t6, "唯一錨存在")
ck("G2-唔重推", True, "本檔對 dof 只作引用，零重推、零數值重算")

# ══════════════════════════════════════════════════════════════════
# ③ v10 P4 docs/73 之 60 dof 硬牆數字
# ══════════════════════════════════════════════════════════════════
print("\n③ v10 P4 docs/73 硬牆數字（只讀檔，唔跑 P4 腳本）")
t73 = rd(P4_73)
ck("G3-節點靈敏度", "0.2348 σ/μK" in t73, "節點檢定 = 一階")
ck("G3-梳狀靈敏度", "**−0.0032 σ/μK**" in t73, "梳狀 ≈ 0（輕微反號）")
ck("G3-比值", "−0.01×" in t73, "靈敏度比 −0.01×（差 ~70 倍）")
ck("G3-60dof", "上限 = 60 dof（本輪已實測證實）" in t73, "60 dof 硬牆")
ck("G3-全域場缺口", "必須由 IDCM 供出 **全域場 L(n̂)**" in t73, "GF1 = 全域場缺口（本 M 承接）")
ck("G3-二階理由", "`A²` 級" in t73 or "**A² 級**" in t73, "局部信號對 w(θ) 係二階")
t72 = rd(P4_72)
ck("G3-2dof", "**orbit 族有 2 個自由度**" in t72, "A₅ orbit 2 dof（種子）")
ck("G3-orbit判據", "B 才係 orbit 之定義" in t72, "orbit 判據 = 保集群（非偏差）")

# ══════════════════════════════════════════════════════════════════
# ④ 零 import 證據
# ══════════════════════════════════════════════════════════════════
print("\n④ 零 import 上游（機械證明：本檔 imports 白名單）")
ALLOWED = {"json", "os", "re", "sys", "fractions", "itertools", "mpmath", "sympy"}
import ast  # noqa: E402

src = rd(os.path.abspath(__file__))
mods = set()
for node in ast.walk(ast.parse(src)):
    if isinstance(node, ast.Import):
        mods.update(a.name.split(".")[0] for a in node.names)
    elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
        mods.add(node.module.split(".")[0])
ck("G4-import白名單", mods <= ALLOWED | {"ast"},
   f"imports = {sorted(mods)}")
FORBIDDEN = ["ex" + "ec(", "sub" + "process", "run" + "py"]
ck("G4-零exec", not any(tok in src for tok in FORBIDDEN),
   "無 " + "／".join(tok.replace("(", "") for tok in FORBIDDEN) + " ⟹ 結構上唔可能重跑上游腳本")

# ══════════════════════════════════════════════════════════════════
# ⑤ 下游五條結構數字（零 import 重算；唔係重推定理）
# ══════════════════════════════════════════════════════════════════
print("\n⑤ 下游結構數字（零 import 重算）")

N = 88
Nd = 44
B = 4
g = [mp.mpf(4) * mp.sin(mp.pi * n / N) ** 2 for n in range(0, Nd + 1)]
ck("G5-γ1", abs(g[1] - mp.mpf("0.005095770779492917")) < mp.mpf("1e-18"),
   f"γ₁ = {mp.nstr(g[1], 18)}")
ck("G5-對極", max(abs(g[n] + g[Nd - n] - B) for n in range(1, Nd)) < mp.mpf("1e-35"),
   "max |γ_n + γ_{44−n} − 4| < 1e-35（n=1..43）")
ck("G5-γ44", abs(g[44] - 4) < mp.mpf("1e-35"), "γ₄₄ = B = 4")
ck("G5-Σγ", abs(sum(g[1:]) - 90) < mp.mpf("1e-33"),
   f"Σγ = {mp.nstr(sum(g[1:]), 20)}（恆等式 2(N_d+1) = 90）")
ck("G5-γ22", abs(g[22] - 2) < mp.mpf("1e-35"), "γ₂₂ = 2（中項半最大）")

# A₅ 度集唯一性（M21 CA1）：五個整數 d、Σd² = 60、Σd = 16 ⟹ {1,3,3,4,5}
sols = []
for combo in combinations_with_replacement(range(1, 17), 5):
    if sum(combo) == 16 and sum(d * d for d in combo) == 60:
        sols.append(combo)
ck("G5-CA1", sols == [(1, 3, 3, 4, 5)], f"5 件度集唯一解 = {sols}")

dset = [1, 3, 3, 4, 5]
ck("G5-Σd²", sum(d * d for d in dset) == 60, "Σd² = 60 = |G|（Burnside）")
ck("G5-N_d", sum(d * (d - 1) for d in dset) == Nd, "Σd(d−1) = 44 = N_d")
ck("G5-16", sum(dset) == 16 == B * B, "Σd = 16 = B²")
ck("G5-60", 60 == Nd + B * B, "60 = N_d + B²（M15/M16 記賬）")

# ══════════════════════════════════════════════════════════════════
# 匯總
# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print(f"\n=== 重現閘：{npass} PASS ／ {nfail} FAIL（共 {len(CHECKS)}）===")
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "00_reproduction_gate.json"), "w", encoding="utf-8") as f:
    json.dump({"checks": [{"tag": t, "pass": p, "detail": d} for t, p, d in CHECKS],
               "n_pass": npass, "n_fail": nfail}, f, ensure_ascii=False, indent=1)
print("EXIT = " + ("0" if nfail == 0 else "1"))
sys.exit(0 if nfail == 0 else 1)
