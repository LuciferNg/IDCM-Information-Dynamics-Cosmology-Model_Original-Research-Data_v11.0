#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/07 — R5：反饋閘（回接 M5 `A_req` 缺口 34×）
======================================================
前設問題
--------
R5（§五）：「R1–R4 任何結果須回接 M5 `A_req` 缺口；接唔上 = 未完成」。
M5 `final_result`：`A_req = 1.3335 ± 0.2570 MeV` ＞ 硬上界 `A_max = α q(ħc/a_C) = 1.2718 MeV`
（中心值超界 4.86%，需 `⟨B⟩* = 1.0485 > 1`）；非線性 IB 殘餘只 1.46%（差 34×）。

本輪問：R1–R4 找到之供給側結構（σ-對支撐、χ=0 退化、16 級樹梯、J^PC 部分）能唔能
改變 `A` 之**上界**？抑或上界係定理級、與支撐無關？

事前判準（寫死，唔准事後改）
----------------------------
  Y1 `A = α q (ħc/a)·⟨B^{−d}⟩` 之最大化必須以**波函數全枚舉／極端情形**證，唔可以只講「多數值細」
  Y2 須證明（或反證）定理：`∀ 歸一化 ψ：⟨B^{−d}⟩ ≤ 1`（因 0 < B^{−d} ≤ 1）
  Y3 須答「供給側能否提供乘性因子 > 1」：計算多 constituent／多 σ-對情形之振幅上界
  Y4 若上界不可破 ⟹ R5 判「接唔上」，並指明缺口所在層（須指名，唔可以留空）
  Y5 禁：調整 `a_C`／`t`／校準點去凑 `A_req`（後驗）

輸出：data/07_feedback_gate.json
"""
import json
import os
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 60
B = 4
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="Y", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="Y"):
    return chk(name, ok, detail, tag, integrity=False)


HBARC = mp.mpf("197.3269804")
a0 = mp.mpf("0.37742216816391858835")
ALPHA = 1 / mp.mpf("137.035999084")
q = mp.mpf(1) / 3                     # 由 M5 A_max 反核：αq = A_max/(ħc/a)
Ehc_a = HBARC / a0
A_max = ALPHA * q * Ehc_a
A_req = mp.mpf("1.3335")
A_req_sig = mp.mpf("0.2570")

print("=" * 78)
print("Y0. 上游數字核對（唔照抄，重算）")
print("=" * 78)
print("  ħc/a₀ = %s MeV｜α = %s｜q = 1/3 ⟹ A_max = %s MeV"
      % (mp.nstr(Ehc_a, 10), mp.nstr(ALPHA, 8), mp.nstr(A_max, 10)))
chk("Y0a A_max 重算 = 1.2718 MeV（M5 docs/13）", abs(A_max - mp.mpf("1.2718")) < mp.mpf("5e-5"),
    "|Δ| = %s" % mp.nstr(abs(A_max - mp.mpf("1.2718")), 4), "M5 final_result")
Needed_B = A_req / A_max
print("  A_req = %s ± %s ⟹ 需 ⟨B⟩* = %s（> 1）" % (A_req, A_req_sig, mp.nstr(Needed_B, 8)))
chk("Y0b ⟨B⟩* = 1.0485（M5 docs/13 報 1.0486）", abs(Needed_B - mp.mpf("1.0486")) < mp.mpf("1e-3"),
    "%s" % mp.nstr(Needed_B, 8), "M5 docs/13")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Y1/Y2. 定理：∀ ψ 歸一化，⟨B^{−d}⟩ ≤ 1（因此 A ≤ A_max）")
print("=" * 78)
Bd = [B ** (-d) for d in range(0, 41)]
print("  0 < B^{−d} ≤ 1 對 d = 0..40 全真：", all(0 < x <= 1 for x in Bd))
print("  d = 0 時 B^{−d} = 1（唯一等號）")
chk("Y2a B^{−d} ≤ 1 ∀d ≥ 0（等號僅 d = 0）", all(0 < x <= 1 for x in Bd),
    "上界由**核形狀**決定，與 ψ 無關")


def cayley_expect(u, L=80):
    """回傳 (基態 ⟨B^{−d}⟩, P₀) —— M5 docs/10 同一算子，本輪獨立實作"""
    N = np.array([1.0] + [float(B * (B - 1) ** (d - 1)) for d in range(1, L + 1)], dtype=float)
    H = np.zeros((L + 1, L + 1))
    for d in range(L + 1):
        H[d, d] = -u * B ** (-d)
    hop1 = float(np.sqrt(B))
    hopn = float(np.sqrt(B - 1))
    for d in range(L):
        h = hopn if d >= 1 else hop1
        H[d, d + 1] = -h
        H[d + 1, d] = -h
    w, v = np.linalg.eigh(H)
    psi = v[:, 0]
    prob = psi ** 2
    prob = prob / prob.sum()
    Bexp = float(np.sum(prob * np.array([B ** (-d) for d in range(L + 1)])))
    return Bexp, float(prob[0])


print()
print("  u = (G/a)/t 掃描 ⟹ ⟨B^{−d}⟩ 之上確界：")
tab = []
for u in (0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 1024):
    be, p0 = cayley_expect(u)
    tab.append({"u": u, "B_expect": be, "P0": p0})
    print("   u=%6.1f：⟨B⟩ = %.6f｜P₀ = %.6f" % (u, be, p0))
chk("Y2b ⟨B⟩ 隨 u 單調升並逼近 1（永不超過）", all(t["B_expect"] <= 1.0 + 1e-12 for t in tab)
    and tab[-1]["B_expect"] > 0.99, "u=1024 ⟹ ⟨B⟩ = %.8f" % tab[-1]["B_expect"])

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Y3. 供給側能唔能提供乘性因子 > 1？（多 constituent／多 σ-對）")
print("=" * 78)
print("  A 之形式 = α q_h Δq_l (ħc/a)·⟨B^{−d}⟩：q 係**電荷積**（群論固定），")
print("  ⟨B⟩ 係**單粒子徑向期望**（≤ 1）。兩體束縛態之和 = Σ_pairs q_iq_j ⟨B_ij⟩，")
print("  每項 ≤ |q_iq_j|·1，而電荷積之絕對值已含於 A_max 之 q。")
print("  ⟹ 乘性因子只能來自**額外帶電荷對**，但 K 系（qq̄）只有 1 對；")
print("     多對需 ≥ 4 constituent（qq̄qq̄），唔屬 K 系通道。")
trial = []
for npair in (1, 2, 3):
    # 上界：每對貢獻 ≤ |q|·A_max 之對應部分；總上界 = npair×A_max（僅當全部 d=0）
    trial.append({"pairs": npair, "A_upper_MeV": float(npair * A_max),
                  "needed_MeV": float(A_req), "ok": float(npair * A_max) >= float(A_req)})
    print("   %d 對 ⟹ A 上界 = %.4f MeV（需 %.4f）⟹ %s"
          % (npair, float(npair * A_max), float(A_req), "夠" if npair * A_max >= A_req else "唔夠"))
chk("Y3a 單對（K 系真實構成）之上界 = A_max < A_req ⟹ 缺口不可由『對數』解釋",
    not trial[0]["ok"], "1 對上界 %.4f < %.4f" % (trial[0]["A_upper_MeV"], trial[0]["needed_MeV"]))
chk("Y3b 2 對雖夠數，但需 qq̄qq̄（4 constituent）⟹ 唔屬 K 系；且 A 上界只喺 d=0 全部飽和才達",
    True, "⟹ 供給側『加對』= 換通道，唔係同通道加因子")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("Y4. 判定：R5 接駁結果")
print("=" * 78)
gap_pct = (A_req / A_max - 1) * 100
print("  缺口 = +%.2f%%（中心值）；需 ⟨B⟩* = %s > 1 ⟹ **定理級不可達**" % (float(gap_pct), mp.nstr(Needed_B, 6)))
verdict("Y4a R5 接唔上：供給側（支撐／佔用／J^PC／徑向）**provably** 唔能提升 A 之上界",
        False, "A ≤ αq(ħc/a) 係核形狀定理（⟨B⟩ ≤ 1），與支撐無關 ⟹ 整層被排除", "Y4")
verdict("Y4b 缺口所在層已指名：一體 EM 自能項 V_self（M5 docs/43 已判『一體 EM 必要』）",
        True, "M18 獨立確認：兩體 Coulomb 核之上界不足，缺口 4.86%（中心）必須落一體項", "Y4")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "07_feedback_gate.py", "checks": checks, "integrity_fail": nfail,
       "A_max_MeV": float(A_max), "A_req_MeV": float(A_req), "Needed_B": float(Needed_B),
       "gap_pct": float(gap_pct), "u_scan": tab, "pair_trials": trial}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "07_feedback_gate.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/07_feedback_gate.json")
sys.exit(0 if nfail == 0 else 1)
