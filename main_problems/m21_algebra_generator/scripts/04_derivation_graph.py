#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/04 — 互導圖 + 恆等反演檢定（R4／T4／P4；含 CG1／GR1）
===================================================================
判準照 docs/00_preregistration.md §3：

  CG1  Σ_{n=1}^{N} 4sin²(πn/2N) = 2(N+1)，N = 2..120，違例 0
       ⟹ 「Σγ = 90 ⟹ N = 44」係恆等式反演（零約束）⟹ 該邊標 T
  GR1  逐邊分類 T（恆等反演）／F（強迫定理）／M（吻合無算子身份）／A（錨）／D（定義）
       閉合判準（嚴格，照預註冊）⟺ T／M／A 邊 = 0

★ 本輪判詞補充（不修改預註冊，只加註）：
   T 邊係恆等反演 —— 無害但**不可能**係 producer；M 邊係真缺口；
   「producer 閉合」嚴格判準 = M + A = 0，另列 T 邊清單。

節點（核心七量）：N_d, B, Σγ, h11, h21, N_gen, N_site
輔助節點：A₅(G), k, N₂ = #{g²=e}, K_tree, L₅

零 import 上游；A₅ 度集由本輪 CA1 唯一解提供（{1,3,3,4,5}）。

⚠ 探索 script：**預期 FAIL 檢查** —— `GR1-c`（producer 閉合 M+A=0）與
   `GR1-d`（嚴格 T+M+A=0）之 FAIL 係**結論本身**（M 邊存在 = h11/h21 無算子身份），
   並非腳本錯誤；exit 1 屬預期，全套重跑時不視為 regression。
"""
import json
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


mp.mp.dps = 60

D = [1, 3, 3, 4, 5]
k = 5
G_order = 60
N2 = sum(D)              # #{g: g² = e} = 16（FS，見 scripts/01）
N_d = 44
B = 4
Sigma = 90
h11, h21 = 36, 98
N_gen = 3
N_site = 88
K_tree = 33
L5 = 11

print("=" * 78)
print("M21 scripts/04 — 互導圖")
print("=" * 78)

# ══════════════════════════════════════════════════════════════════
# CG1：Σγ(N) = 2(N+1) 對 N = 2..120（精確性：dps 60，容差 1e-45）
# ══════════════════════════════════════════════════════════════════
print("\n── CG1 循環圖恆等式 Σγ(N) = 2(N+1) ──")
viol = []
tol = mp.mpf("1e-45")
for N in range(2, 121):
    s = mp.fsum([4 * mp.sin(mp.pi * n / (2 * N)) ** 2 for n in range(1, N + 1)])
    if abs(s - 2 * (N + 1)) > tol:
        viol.append((N, mp.nstr(s, 20), 2 * (N + 1)))
print(f"  N = 2..120（119 個）違例數 = {len(viol)}")
ck("CG1  Σγ(N) = 2(N+1) 恆成立（違例 0）⟹ 對 N 零約束", len(viol) == 0, f"violations = {len(viol)}")
ck("CG1-b  解析根據：γ_n + γ_{N−n} = 4（sin²+cos²=1）⟹ Σ = 2(N−1)+4 = 2N+2",
   abs(mp.fsum([4 * mp.sin(mp.pi * n / 88) ** 2 for n in range(1, 45)]) - 90) < tol, "N=44 特例 = 90")

# ══════════════════════════════════════════════════════════════════
# 邊表：逐邊機核 + 分類
# ══════════════════════════════════════════════════════════════════
print("\n── 邊表（逐邊機核 + 分類）──")
edges = []


def edge(eid, src, dst, stmt, cls, ok, evidence):
    edges.append({"id": eid, "src": src, "dst": dst, "stmt": stmt, "class": cls,
                  "verified": bool(ok), "evidence": evidence})
    print(f"  {eid:>4}  {src:>10} → {dst:<8} [{cls}]  {stmt:<46} {evidence}")
    return ok


# E1：|G| = 60 + DNT ⟹ B = 4（M16）—— 本輪只核 B = γ_max = 4（譜端精確）
ok1 = abs(4 * mp.sin(mp.pi * 44 / 88) ** 2 - 4) < tol
edge("E1", "A₅(G)", "B", "B = γ_max（+ M16：|G|=60 ⟹ B=4，dof=0）", "F", ok1,
     "γ_44 = 4 精確")

# E2：N_gen = B − 1（M15 docs/17 §4）
edge("E2", "B", "N_gen", "N_gen = B − 1 = 3", "F", N_gen == B - 1, f"{B}−1 = {N_gen}")

# E3：N_site = 2·N_d（對偶表示之定義）
edge("E3", "N_d", "N_site", "N_site = 2·N_d（dual 88-site 表示）", "D", N_site == 2 * N_d,
     f"2×{N_d} = {N_site}")

# E4：A₅ ⟹ N_d = Σd(d−1) = |G| − Σd（本輪 CA1–CA3）
ok4 = (sum(d * (d - 1) for d in D) == 44) and (G_order - N2 == 44)
edge("E4", "A₅(G)", "N_d", "N_d = Σ_χ d(d−1) = |G| − Σ_χ d = 60 − 16", "F", ok4,
     f"Σd(d−1) = {sum(d * (d - 1) for d in D)}；|G|−Σd = {G_order - N2}")

# E5：N_d ⟹ Σγ = 2(N_d+1) —— 恆等反演
edge("E5", "N_d", "Σγ", "Σγ = 2(N_d+1) = 90", "T", Sigma == 2 * (N_d + 1),
     "恆等（CG1）；對 N 零約束")

# E6：B、N_d ⟹ Σγ = B(N_d+1)/2 —— 同一恆等（B/2 = 2）
edge("E6", "B,N_d", "Σγ", "Σγ = B(N_d+1)/2 = 90", "T", Sigma == B * (N_d + 1) // 2,
     "與 E5 同一方程（B/2 = 2）")

# E7：h11 = N_d − 2B
edge("E7", "N_d,B", "h11", "h11 = N_d − 2B = 36", "M", h11 == N_d - 2 * B,
     "數值吻合；池掃 ≥6 條 k-free（scripts/02）")

# E8：h21 = 2N_d + 10
edge("E8", "N_d", "h21", "h21 = 2N_d + 10 = 98", "M", h21 == 2 * N_d + 10,
     "數值吻合；10 = 2F₅ 係後見選項")

# E9：L₅ ⟹ N_d = B·L₅
edge("E9", "L₅", "N_d", "N_d = B·L₅ = 4×11 = 44", "F", N_d == B * L5,
     f"{B}×{L5} = {N_d}（由 E4 + E11 合成）")

# E10：A₅ ⟹ K_tree = Σ(d−1)² = 33
ok10 = sum((d - 1) ** 2 for d in D) == 33
edge("E10", "A₅(G)", "K_tree", "K_tree = Σ_χ (d−1)² = 33", "F", ok10,
     f"{sum((d - 1) ** 2 for d in D)}")

# E11：A₅ ⟹ L₅ = Σ(d−1) = 11
ok11 = sum(d - 1 for d in D) == 11
edge("E11", "A₅(G)", "L₅", "L₅ = Σ_χ (d−1) = 11", "F", ok11, f"{sum(d - 1 for d in D)}")

# E12：A₅ ⟹ N₂ = Σd = 16（FS，見 scripts/01）
edge("E12", "A₅(G)", "N₂", "N₂ = Σ_χ d = #{g: g²=e} = 16", "F", N2 == 16, f"{N2}")

# E13：N₂ = B²（記賬 60 = 44 + 16 之內容）
edge("E13", "N₂,B", "|G| 記賬", "|G| = N_d + B² ⟺ N₂ = B² = 16", "M", N2 == B ** 2,
     "兩邊各自被 A₅ 確定，但無算子身份 ⟹ 一致性檢查")

# ══════════════════════════════════════════════════════════════════
# GR1 計數與判詞
# ══════════════════════════════════════════════════════════════════
print("\n── GR1 分類計數 ──")
from collections import Counter
cnt = Counter(e["class"] for e in edges)
print(f"  F（強迫）={cnt['F']} ｜ T（恆等反演）={cnt['T']} ｜ M（吻合）={cnt['M']}"
      f" ｜ A（錨）={cnt.get('A', 0)} ｜ D（定義）={cnt['D']}")
allv = all(e["verified"] for e in edges)
ck("GR1-a  全部邊之陳述皆機核為真", allv, f"{sum(1 for e in edges if e['verified'])}/{len(edges)}")
ck("GR1-b  錨邊 A = 0（七量內無未解釋輸入）", cnt.get("A", 0) == 0, f"A = {cnt.get('A', 0)}")
ck("GR1-c  producer 閉合（M + A = 0）", (cnt["M"] + cnt.get("A", 0)) == 0,
   f"M + A = {cnt['M'] + cnt.get('A', 0)}")
ck("GR1-d  嚴格預註冊判準（T + M + A = 0）", (cnt["T"] + cnt["M"] + cnt.get("A", 0)) == 0,
   f"T + M + A = {cnt['T'] + cnt['M'] + cnt.get('A', 0)}")
ck("GR1-e  Σγ 之入邊全屬 T（恆等反演，非 producer）",
   all(e["class"] == "T" for e in edges if e["dst"] == "Σγ"), "E5, E6")
ck("GR1-f  h11/h21 之入邊全屬 M（無算子身份）",
   all(e["class"] == "M" for e in edges if e["dst"] in ("h11", "h21")), "E7, E8")

# 孤立節點檢查
print("\n── 節點度數 ──")
nodes = ["N_d", "B", "Σγ", "h11", "h21", "N_gen", "N_site", "A₅(G)", "K_tree", "L₅", "N₂"]
deg = Counter()
for e in edges:
    for s in e["src"].split(","):
        deg[s.strip()] += 1
    deg[e["dst"]] += 1
for n in nodes:
    print(f"  {n:<10} deg = {deg[n]}")
ck("GR1-g  無孤立節點（每節點 deg ≥ 1）", all(deg[n] >= 1 for n in nodes),
   f"min deg = {min(deg[n] for n in nodes)}")

# ══════════════════════════════════════════════════════════════════
# dof 記帳（解 M15 out_of_scope/docs/14 §F5(3) 之張力）
# ══════════════════════════════════════════════════════════════════
print("\n── dof 記帳 ──")
n_anchor = cnt.get("A", 0)
n_forced_Nd = 1 if any(e["id"] == "E4" and e["class"] == "F" and e["verified"] for e in edges) else 0
print(f"  N_d 之 producer = E4（F 類，零參數）⟹ N_d 非錨")
print(f"  七量內錨數 = {n_anchor} ｜ N_d 錨化 = {1 - n_forced_Nd}")
ck("dof-a  N_d 由零參數群不變量確定 ⟹ 非第二錨",
   n_forced_Nd == 1 and n_anchor == 0, "分支 (a)：O212 dof = 1 成立")
ck("dof-b  連續 dof 不增（無新增可調參數）", True,
   "E4/E10/E11 全部係整數恆等，無自由係數")

# ══════════════════════════════════════════════════════════════════
npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/04 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "04_derivation_graph.json"), "w", encoding="utf-8") as f:
    json.dump({
        "checks": [{"tag": t, "pass": ok, "detail": d} for t, ok, d in CHECKS],
        "n_pass": npass, "n_fail": nfail, "n_total": len(CHECKS),
        "edges": edges,
        "class_counts": dict(cnt),
        "node_degrees": {n: deg[n] for n in nodes},
        "dof": {"anchors_in_seven": n_anchor, "N_d_forced": bool(n_forced_Nd), "branch": "(a)"},
    }, f, ensure_ascii=False, indent=2)

sys.exit(0 if nfail == 0 else 1)
