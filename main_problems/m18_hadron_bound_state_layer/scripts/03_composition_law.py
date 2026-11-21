#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/03 — R2：mode 佔用選擇律（T2）
==========================================
前設問題
--------
M4 Q2/Q3/Q8 阻塞點：邊個 mode 組合構成一個強子？M5 判「缺口必在供給側（mode 佔用）」。
R2 問：`mode 佔用` 之選擇律能否由 γ_n + γ_{44−n} = 4 之對極結構導出（零觀測輸入）？

本輪採用**譜閉合判準**（零新假設之最簡形式）：
  一個 N_q 體複合態嘅內容 S（multiset，每個 constituent 佔一個鏈模）為「在殼」
  ⟺ 其能量量子 E(S) = Σ_{n∈S} γ_n **本身係譜值** γ_m（m ∈ I₄₄）
理由（結構，非擬合）：𝒯 之譜係離散封閉集；複合態若唔落在譜上，則無本徵模可承載 ⟹
「唔存在」。此判準**唔含任何可調參數**，亦唔用任何觀測質量。

已註冊對照（M14 三波共振定律，N-無關）：
  γ_a + γ_b = γ_m 有解；sum 通道 ⟺ a+b = N_c（落 m = N_c）；
  γ_b − γ_a = γ_m 有解；beat 通道 ⟺ b = N_c（落 m = N_c − a）
  本輪須**獨立重算**（唔 import M14 scripts），並核數字。

事前判準（寫死，唔准事後改）
----------------------------
  S1 全部命中必須以 mpmath dps=60（tol 1e−45）掃，再以 ℚ(ζ₈₈) 環元素**精確確認**
     （避免 float 預篩漏解 —— M14 已兩次踩過）
  S2 兩體：全部 a ≤ b（946 對）× 全部 m ∈ I₄₄ ⟹ 完全枚舉
  S3 三體：全部 combinations_with_replacement(1..44, 3) = 15180 條，完全枚舉（本輪新算）
  S4 判定律必須係**閉式**（可對任意 N_c 陳述），並做 N_c 剛性掃描
  S5 選擇律只可用譜結構；任何用到觀測質量之步驟 ⟹ 該支結果即時作廢
  S6 必須報「選擇性比率」= 在殼內容 / 全部內容；若比率 ≈ 1 ⟹ 無鑑別力

輸出：data/03_composition_law.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
TOL = mp.mpf("1e-45")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
NC = 44
checks = []


def chk(name, ok, detail, tag="S", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="S"):
    return chk(name, ok, detail, tag, integrity=False)


x = sp.symbols("x")
PHI88 = sp.Poly(sp.cyclotomic_poly(88, x), x, domain=sp.QQ)


def ring(p):
    return sp.Poly(p, x, domain=sp.QQ).rem(PHI88)


def gam_ring(n):
    return ring(sp.Poly(2, x, domain=sp.QQ) - sp.Poly(x**n + x**(88 - n), x, domain=sp.QQ))


def gam(n):
    return mp.mpf(4) * mp.sin(mp.pi * n / 88) ** 2


G = {n: gam(n) for n in range(1, 45)}
GR = {n: gam_ring(n) for n in range(1, 45)}
spec_idx = {}          # 數值 → m（唯一化，tol 內）
for m in range(1, 45):
    spec_idx.setdefault(mp.nint(G[m] / TOL * mp.mpf("1e-45")) if False else m, m)
vals = sorted((G[m], m) for m in range(1, 45))


def on_shell(E):
    """E 係否譜值；返回所有 m（tol 內），再精確確認另行做"""
    out = []
    for v, m in vals:
        if abs(E - v) < TOL:
            out.append(m)
    return out


def on_shell_exact(E_ring):
    return [m for m in range(1, 45) if (E_ring - GR[m]).is_zero]


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("S2. 兩體 + 通道（sum，複合）：全部 946 對，完全枚舉")
print("=" * 78)
sum_hits = []
for a in range(1, 45):
    for b in range(a, 45):
        E = G[a] + G[b]
        ms = on_shell(E)
        if ms:
            sum_hits.append((a, b, tuple(ms)))
print("  命中數 =", len(sum_hits))
print("  樣本 =", sum_hits[:6], "…", sum_hits[-3:])
mset = sorted({m for h in sum_hits for m in h[2]})
chk("S2a 兩體 + 通道全部落 γ₄₄ = 4（單一 m）", mset == [44], "m 值集 = %s" % mset)
pred_sum = [(a, 44 - a, (44,)) for a in range(1, 23)]
chk("S2b 解集 = {a + b = 44}（22 條，含 (22,22)）", sorted(sum_hits) == sorted(pred_sum),
    "差集 = %s" % (set(sum_hits) ^ set(pred_sum)))
# 精確確認
exact_all = all(on_shell_exact(GR[a] + GR[b]) == [44] for a, b, _ in sum_hits)
chk("S2c 22 條全部 ℚ(ζ₈₈) 精確確認", exact_all, "環確認 22/22")
# 非命中樣本的反證（精確：無其他 (a,b,m) 解）
exact_extra = []
for a in range(1, 45):
    for b in range(a, 45):
        if (a, 44 - a) == (a, b) or b == 44 - a:
            continue
        if a + b == 44:
            continue
        if on_shell_exact(GR[a] + GR[b]):
            exact_extra.append((a, b))
chk("S2d 非 (a+b=44) 之 924 對：環內 0 其他解", len(exact_extra) == 0, "額外解 %d" % len(exact_extra))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("S2'. 對照：M14 已註冊 law 之獨立重算（+ 通道 21/43、− 通道）")
print("=" * 78)
beat_hits = []
for a in range(1, 45):
    for b in range(a + 1, 45):
        E = G[b] - G[a]
        if E <= 0:
            continue
        ms = on_shell(E)
        if ms:
            beat_hits.append((a, b, tuple(ms)))
print("  − 通道命中數 =", len(beat_hits), "｜樣本 =", beat_hits[:5])
beat_ok = all(b == 44 and set(ms) == {44 - a} for a, b, ms in beat_hits)
chk("S2'a − 通道全部 = {b = 44，m = 44−a}（43 條）", beat_ok and len(beat_hits) == 43,
    "命中 %d；b=44 全真 = %s" % (len(beat_hits), beat_ok))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("S3. 三體 + 通道（baryon 供給）：15180 條 multiset，完全枚舉")
print("=" * 78)
tri_hits = []
for S in itertools.combinations_with_replacement(range(1, 45), 3):
    E = G[S[0]] + G[S[1]] + G[S[2]]
    ms = on_shell(E)
    if ms:
        tri_hits.append((S, tuple(ms)))
print("  命中數 =", len(tri_hits), "｜全部內容 =", 15180)
ms3 = {}
for S, ms in tri_hits:
    for m in ms:
        ms3[m] = ms3.get(m, 0) + 1
print("  m 分佈 =", dict(sorted(ms3.items())))
print("  樣本 =", tri_hits[:8])
# 獨立方法核掃描正確性（ordered tuples 44³ = 85184，tol 1e−50）＋ 最小間隙量度
mind = mp.mpf(10)**9
mind_S = None
n_hit_ordered = 0
for a in range(1, 45):
    for b in range(1, 45):
        for c in range(1, 45):
            E = G[a] + G[b] + G[c]
            d = min(abs(E - v) for v, _ in vals)
            if d < mind:
                mind, mind_S = d, (a, b, c)
            if d < mp.mpf("1e-50"):
                n_hit_ordered += 1
print("  獨立核（ordered 44³=85184）：命中 %d｜最小間隙 = %s @ %s"
      % (n_hit_ordered, mp.nstr(mind, 6), mind_S))
chk("S3a0 獨立方法一致（ordered scan 亦 0 命中；最小間隙 ≫ tol）",
    n_hit_ordered == 0 and mind > mp.mpf("1e-12"),
    "最小間隙 %s ⟹ 零命中穩健（唔係容差 artifact）" % mp.nstr(mind, 6))
verdict("S3a 三體譜閉合 = 空集（0/15180 在殼內容）", len(tri_hits) == 0,
        "★ 結構 no-go：鏈譜對三體加法**唔封閉** ⟹ baryon 唔可能以『譜閉合』承載")
# 精確確認全部命中
ex3 = sum(1 for S, ms in tri_hits if on_shell_exact(GR[S[0]] + GR[S[1]] + GR[S[2]]))
chk("S3b 全部命中 ℚ(ζ₈₈) 精確確認", ex3 == len(tri_hits), "環確認 %d/%d" % (ex3, len(tri_hits)))
chk("S3c 選擇性比率 = 命中/全部", True,
    "%.3f%%（%d/%d）" % (100.0 * len(tri_hits) / 15180, len(tri_hits), 15180))

# 結構刻畫：命中集是否 = {σ-pair + 22} 或含 44？
def chi(S):
    return sum(G[n] - 2 for n in S)


with44 = sum(1 for S, _ in tri_hits if 44 in S)
has22 = sum(1 for S, _ in tri_hits if 22 in S)
chi0 = sum(1 for S, _ in tri_hits if abs(chi(S)) < TOL)
print("  含 44 之命中 = %d｜含 22 之命中 = %d｜χ=0 之命中 = %d" % (with44, has22, chi0))
for S, ms in tri_hits[:12]:
    print("    %-16s χ = %+8.4f  m = %s" % (str(S), float(chi(S)), ms))

print()
print("=" * 78)
print("S4. 結構剛性：兩體 law 對 N_c 掃描（{3,5,6,11,22,44}）")
print("=" * 78)
rig = []
for Nc in (3, 5, 6, 11, 22, 44):
    gg = [None] + [mp.mpf(4) * mp.sin(mp.pi * n / (2 * Nc)) ** 2 for n in range(1, Nc + 1)]
    hits = [(a, b) for a in range(1, Nc + 1) for b in range(a, Nc + 1)
            if min(abs(gg[a] + gg[b] - gg[m]) for m in range(1, Nc + 1)) < TOL]
    pred = [(a, Nc - a) for a in range(1, Nc // 2 + 1)]
    rat = [n for n in range(1, Nc + 1) if abs(gg[n] - round(gg[n])) < mp.mpf("1e-40")]
    rig.append({"N_c": Nc, "hits": len(hits), "pred": len(pred),
                "extra": len(hits) - len(pred),
                "rational_interior": [n for n in rat if n != Nc],
                "law_holds": set(pred) <= set(hits)})
    print("   N_c=%2d：命中 %3d｜a+b=N_c %3d｜額外 %d｜內部有理譜值 %s"
          % (Nc, len(hits), len(pred), len(hits) - len(pred),
             [n for n in rat if n != Nc]))
chk("S4a a+b=N_c 之解必屬命中集（必要部分，6/6）", all(r["law_holds"] for r in rig), "6/6")
verdict("S4b 命中集 ≡ {a+b=N_c} ⟺ 無內部有理譜值",
        all(r["extra"] == 0 for r in rig if not r["rational_interior"]),
        "N_c=6 例外 +%d 條（γ₂=1、γ₃=2 有理 ⟹ 2γ₂=γ₃、γ₂+γ₃=γ₄）⟹ 定律**精確性**需 Niven 級條件"
        % [r["extra"] for r in rig if r["N_c"] == 6][0])
chk("S4c N_c=44 唯一內部有理譜值 = γ₂₂ = μ² = 2，其加法巧合 2μ² = 4 = γ₄₄ 已含於 {a+b=N_c}",
    [r["rational_interior"] for r in rig if r["N_c"] == 44][0] == [22],
    "γ_n 有理 ⟺ n ∈ {0,22,44}（γ = 0,2,4）；22 係唯一內部有理值；"
    "2γ₂₂ = γ₄₄ 恰好就係 (22,22) ∈ {a+b=44} ⟹ 零額外巧合（與 N_c=6 不同）")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "03_composition_law.py", "checks": checks, "integrity_fail": nfail,
       "sum_hits": len(sum_hits), "sum_all_m44": mset == [44],
       "beat_hits": len(beat_hits), "tri_hits": len(tri_hits),
       "tri_total": 15180, "tri_selectivity_pct": 100.0 * len(tri_hits) / 15180,
       "tri_m_dist": {str(k): v for k, v in sorted(ms3.items())},
       "rigidity": rig, "tri_min_gap": str(mp.nstr(mind, 10)), "tri_contain_44": with44, "tri_contain_22": has22, "tri_chi0": chi0}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "03_composition_law.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/03_composition_law.json")
sys.exit(0 if nfail == 0 else 1)
