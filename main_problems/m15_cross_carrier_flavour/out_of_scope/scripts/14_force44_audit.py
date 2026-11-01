#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/14 — 審計上游 `FORCE_44`（v10 p9）：`N_d = 44` 係定理抑或錨？
================================================================================
背景：M15 `docs/13 §Q3` 把 `K_tree` 之 k-free 問題歸約至「為何 N_d = L₃L₅ = 44」，
      並誤記為「上游命題、非本 O 可再推」。實際上 v10 `p9_kchain_number_field/
      05_structural_proof_FORCE_44.{md,py}` 宣稱「✅ 定理證明完成（7/7）」⟹ 必須審。

被審對象（**只讀**）：v10 p9 `05_structural_proof_force_44.{md,py}`。

判準（事前寫死）：
  F1  引理 1「Σγ = 2(N+1)」若對所有 N 成立 ⟹ 提供零約束（只是把 Σγ 寫成 N 之函數）。
  F2  引理 2「90 = L₁₁ − L₉ − L₇ − L₃」之**選擇性**：若同型 4 項交錯和可解釋大量整數 ⟹ 非強迫。
  F3  檢查 3「2(N+1) = 90 ⟹ N = 44」：若 90 係字面輸入 ⟹ 以答案求答案（違 A1 出處閘）。
  F4  腳本內 `N_d = 44` 若硬編碼 ⟹ 全部檢查皆為「在答案上驗恆等」。
  F5  結論：若 F1–F4 全中 ⟹ FORCE_44 非定理；`N_d = 44` 降級為**錨（dof +1）**。
"""
import itertools
import json
import re
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 30
P9 = Path("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat/"
          "main_problems/p9_kchain_number_field")
LOG, RES = [], []


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail=""):
    st = "PASS" if cond else "FAIL"
    say(f"[{st}] {tag}" + (f"  |  {detail}" if detail else ""))
    RES.append({"tag": tag, "status": st, "detail": detail})
    return cond


say("=" * 78)
say("M15 scripts/14 — 上游 FORCE_44 審計：定理抑或錨")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§F1  引理 1 是否恆等（Σ_{n=1}^{N} 4 sin²(πn/2N) = 2N + 2，∀N）")
say("-" * 78)
bad = []
for N in range(2, 121):
    s = sum(4 * mp.sin(mp.pi * n / (2 * N)) ** 2 for n in range(1, N + 1))
    if abs(s - (2 * N + 2)) > mp.mpf("1e-25"):
        bad.append(N)
say(f"  N = 2..120 逐一驗證：違例 {bad if bad else '無'}")
say(f"  N = 44 ⟹ 2(N+1) = {2 * 44 + 2}")
chk("F1  引理 1 對所有 N 恆等 ⟹ 零約束", not bad, f"違例 {bad}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§F2  引理 2 之選擇性（odd Lucas 4 項相異、任意符號）")
say("-" * 78)
L = {1: 1, 3: 4, 5: 11, 7: 29, 9: 76, 11: 199, 13: 521}
vals, idxs = list(L.values()), list(L.keys())
expl = {}
for comb in itertools.combinations(range(len(vals)), 4):
    base = [vals[i] for i in comb]
    for signs in itertools.product([1, -1], repeat=4):
        S = sum(s * b for s, b in zip(signs, base))
        if S > 0 and S % 2 == 0:
            N = S // 2 - 1
            if 3 <= N <= 200:
                expl.setdefault(N, []).append((tuple(idxs[i] for i in comb), signs, S))
say(f"  可解釋之 N（3..200）數目 = {len(expl)} / 198 = {len(expl) / 198 * 100:.1f}%")
say("  最小 8 個被解釋之 N：" + str(sorted(expl)[:8]))
say(f"  N = 44（即 Σγ = 90）之表示數 = {len(expl.get(44, []))}")
for r in expl.get(44, []):
    terms = " ".join(f"{'+' if s > 0 else '−'}L_{i}" for i, s in zip(r[0], r[1]))
    say(f"    {terms} = {r[2]}")
chk("F2  級聯非唯一（90 有 ≥ 2 條同型表示）且解釋率 > 20% ⟹ 非強迫",
    len(expl.get(44, [])) >= 2 and len(expl) / 198 > 0.2,
    f"90 之表示 {len(expl.get(44, []))} 條；解釋率 {len(expl) / 198 * 100:.1f}%")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§F3/F4  被審腳本之靜態審計（字面值當結構）")
say("-" * 78)
src = (P9 / "05_structural_proof_force_44.py").read_text()
lits = {}
m1 = re.search(r"^\s*N_d\s*=\s*(\d+)", src, re.M)
lits = {"N_d =": m1.group(1).strip() if m1 else None,
        "90//2-1（檢查 3）": "found" if re.search(r"90\s*//\s*2\s*-\s*1", src) else None,
        "Σγ(目標)": "found" if re.search(r"Σγ\(目標\)", src) else None}
say(f"  `N_d = {lits['N_d =']}` 硬編碼：{'是' if lits['N_d ='] == '44' else '否'}")
say(f"  檢查 3 之式：`C3 = (90 // 2 - 1 == 44)` ⟹ 對兩個字面值做算術")
say(f"  腳本自印 `Σγ(目標) = 90` ⟹ **90 係目標值，非由 φ 級聯『逼出』**")
say("")
say("  即：")
say("    引理 1：Σγ = 2(N+1)                      [∀N 恆等]")
say("    引理 2：90 = L₁₁ − L₉ − L₇ − L₃            [算術恆等，與 N 無關]")
say("    檢查 3：2(N+1) = 90 ⟹ N = 44              [把 90 當已知目標 ⟹ N 已知]")
say("  ⟹ 三者合起來只證明「若 Σγ = 90 則 N = 44」，而 Σγ = 90 ⟺ N = 44（因 90 = 2·44+2）")
say("  ⟹ **循環**：以答案（44）求得之數（90）反求答案。")
chk("F3  檢查 3 以字面 90 求解 ⟹ 以答案求答案（違 A1）",
    lits["90//2-1（檢查 3）"] is not None, "`90 // 2 - 1 == 44`")
chk("F4  腳本硬編碼 N_d = 44 ⟹ 全部 7 項檢查皆在答案上驗恆等",
    lits["N_d ="] == "44", "N_d = " + str(lits["N_d ="]))

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§F5  判詞 ＋ dof 張力")
say("-" * 78)
say("  (1) FORCE_44 **不成立為定理**：F1 恆等、F2 非選擇（31.8% 之 N 可被同型級聯解釋，")
say("      且 90 本身有 3 條相異表示）、F3 以字面 90 求解、F4 硬編碼 44。")
say("  (2) ⟹ `N_d = 44`（⟺ 88 site）之地位 = **錨（anchor）**，dof +1，不可由 corpus 推出。")
say("  (3) ★ dof 張力：O212（v10 `o212_anchor_theorem`）宣稱「錨 = G ≡ ω_ring ≡ M_Pl（唯一測量）」")
say("      ⟹ dof = 1。若 `N_d = 44` 係第二個獨立錨 ⟹ **dof = 2**。二擇一：")
say("        (a) 88 另有（未註冊之）結構強迫 ⟹ dof 仍 1；")
say("        (b) IDCM 實際 dof = 2 ⟹ O212 之計數需修訂。")
say("      ⚠ 此屬**跨版本**（O212 在 v10）⟹ M15 只可記錄 ＋ 同日版本（v11）M6 通報；禁跨版本推送。")
say("  (4) ⟹ M15 之最終鏈在其末端係**錨**，非推導：")
say("        3 = rank(P_κ) = 60·Ω_b = γ_max − 1 ，其上游 = N_d = 44（錨）")
chk("F5  FORCE_44 判非定理；N_d = 44 降為錨；dof 張力已定位", True,
    "跨版本（O212）⟹ 只記錄；v11 內通報 M6")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "lemma1_violations": bad,
               "n_explainable": len(expl), "explain_rate_pct": round(len(expl) / 198 * 100, 2),
               "representations_for_90": [[list(r[0]), list(r[1]), r[2]] for r in expl.get(44, [])],
               "hardcoded_N_d": lits["N_d ="], "check3_literal": lits["90//2-1（檢查 3）"]},
              f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)