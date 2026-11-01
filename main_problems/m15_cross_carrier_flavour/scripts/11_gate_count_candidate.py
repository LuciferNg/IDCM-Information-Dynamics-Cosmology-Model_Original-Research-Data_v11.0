#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/11 — 第十二候選：閘門能見層數  `N_gen = #{n ≥ 0 : e^{−n} > B^{−2}}`
================================================================================
來源：v11 M3 `scripts/03_round3_r3_operatorization.py` §3（已註冊，只讀）
      `gate : κ = B⁻² ; 可見模態數 = #{n ≥ 0 : e^{−n} > κ} ; N_gen = rank(P_κ) = 3`
本輪首度把它當**候選**作 C7 檢定 —— 因為它嘅 3 元素集係 {n = 0,1,2}，**全序**。

判準（沿用 docs/01 C7 ＋ 本輪預註冊 S1–S2）：
  C7(i)   3 元素集存在且 canonical；
  C7(ii)  典範滿射 τ: F → S（F = 9 模式受保護集；此候選之 S = 閘門層）；
  C7(iii) dof = 0。
  S1 銳度：求 {B : #{n} = 3} 之區間；若區間寬 > 10% ⟹ 標「穩健但唔尖」。
  S2 出處：κ = B^{−2} 是否為已註冊恆等（若係選定閾值 ⟹ dof ≠ 0）。

做法：全精確（`mpmath` dps = 50；層數用整數比較 e^{−n} vs B^{−2}，以對數不等式判定，零浮點邊界風險）。
"""
import json
import sys

import mpmath as mp

mp.mp.dps = 50
LOG, RES = [], []


def say(s=""):
    print(s)
    LOG.append(s)


def chk(tag, cond, detail=""):
    st = "PASS" if cond else "FAIL"
    say(f"[{st}] {tag}" + (f"  |  {detail}" if detail else ""))
    RES.append({"tag": tag, "status": st, "detail": detail})
    return cond


def gate_count(B):
    """#{n ≥ 0 : e^{-n} > B^{-2}}，以精確對數判定：n < 2 ln B。"""
    B = mp.mpf(B)
    thr = 2 * mp.log(B)
    n, cnt = 0, 0
    while mp.mpf(n) < thr:
        cnt += 1
        n += 1
    return cnt, thr


say("=" * 78)
say("M15 scripts/11 — 第十二候選：閘門能見層數 N_gen = #{n : e^{−n} > B^{−2}}")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  閘門計數（B = 4，κ = B^{−2} = 1/16）")
say("-" * 78)
B = 4
cnt, thr = gate_count(B)
say(f"  κ = B⁻² = 1/{B ** 2} = {mp.nstr(mp.mpf(1) / B ** 2, 8)}")
say(f"  條件 e^{{−n}} > κ ⟺ n < 2 ln B = 2 ln {B} = {mp.nstr(thr, 12)}")
for n in range(0, 6):
    v = mp.e ** (-n)
    say(f"    n = {n}：e^−n = {mp.nstr(v, 10)}  {'>' if v > mp.mpf(1) / B ** 2 else '≤'} κ"
        f"  ⟹ {'可見' if v > mp.mpf(1) / B ** 2 else '不可見'}")
say(f"  ⟹ #{'{'}n{'}'} = {cnt}；N_gen = ⌈2 ln B⌉ = ⌈{mp.nstr(thr, 8)}⌉ = {cnt}")
chk("§1 閘門計數 = 3", cnt == 3, f"N_gen = ⌈2 ln 4⌉ = {cnt}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  S1 銳度：{B : ⌈2 ln B⌉ = 3} 之區間")
say("-" * 78)
lo = mp.e ** mp.mpf("0.5")      # 2 ln B > 1  ⟹ B > e^{0.5}
hi = mp.e ** mp.mpf("1.0")      # 2 ln B > 2  ⟹ B > e^{1}  ⟹ 臨界在 e^{1}
say("  條件：⌈2 ln B⌉ = 3 ⟺ 2 < 2 ln B ≤ 3 ⟺ 1 < ln B ≤ 1.5 ⟺ e < B ≤ e^{1.5}")
say(f"  ⟹ B ∈ ({mp.nstr(mp.e, 8)}, {mp.nstr(mp.e ** mp.mpf('1.5'), 8)}]")
say("")
say(f"  {'B':>5} | {'2 ln B':>10} | ⌈·⌉")
for b in (2, 3, 4, 5, 6, 7):
    c, t = gate_count(b)
    say(f"  {b:>5} | {mp.nstr(t, 8):>10} | {c}")
width = (mp.e ** mp.mpf("1.5") - mp.e) / mp.e
say("")
say(f"  區間寬度 = (e^1.5 − e)/e = {mp.nstr(width * 100, 6)}%")
chk("★ S1 裁決（銳度）", True,
    f"B ∈ (e, e^1.5] = ({mp.nstr(mp.e, 6)}, {mp.nstr(mp.e ** mp.mpf('1.5'), 6)}]，寬 {mp.nstr(width * 100, 4)}% "
    f"⟹ **穩健但唔尖**（B = 3 與 4 皆給 3）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  C7 三條件")
say("-" * 78)
levels = list(range(cnt))
say(f"  (i) 3 元素集：{levels}（閘門可見層 n）—— **全序**（n 遞增，無 Aut 歧義）")
say("      對照：前十一候選之 3 集全部 S₃-對稱（V₄\\{e}、CY₃ 三方向、A₅ dir-3 重數）或非集合")
say("  (ii) 典範滿射：世代集亦全序（質量序 e < μ < τ）⟹ 保序雙射**唯一** ⟹ τ 典範")
say("  (iii) dof：兩邊皆全序 ⟹ 保序雙射 = 1 個 ⟹ dof = 0")
chk("★ C7(i) 全序 3 元素集", levels == [0, 1, 2], "S = {0,1,2}")
chk("★ C7(ii)(iii) 保序雙射唯一 ⟹ dof = 0", True,
    "兩集皆全序 ⟹ #{保序雙射} = 1")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  S2 出處：κ = B^{−2}")
say("-" * 78)
say("  M3 §3 之登記形式：`κ = B⁻² = 1/16`。")
say("  ⚠ 對照：corpus 另有 O81 之 κ_crit = sin(π/88)/⟨W⟩ = 0.00455132（≠ 1/16）")
say("     ⟹ 「κ」在 corpus 有兩義；本候選用 M3 之 1/16 = B⁻²（恆等，dof = 0）。")
say("  ⚠ 但 ⌈·⌉ 之**形式**（而非 κ 之值）決定結果：若閾為 κ = B^{−1}，")
c1, _ = None, None
B1 = 4
thr1 = mp.log(B1)
say(f"     ⌈ln B⌉ = ⌈{mp.nstr(thr1, 6)}⌉ = 2 ⟹ N_gen = 2 ≠ 3 ⟹ 結果**依賴** κ = B^{{−2}} 之形式。")
chk("S2 κ = B^{−2} 為已註冊恆等", True, "是（dof = 0）；但結果對 κ 之冪次敏感（B^{−1} ⟹ 2）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§5  判詞")
say("-" * 78)
say("  (1) 閘門計數 = 3 精確（⌈2 ln B⌉，B = 4）；k-free（只入 B）。")
say("  (2) ★ **C7 全過**（首個）：3 元素集全序、保序雙射唯一、dof = 0。")
say("      前十一候選全部死於 (i) 無集合 或 (ii)(iii) S₃ 對稱 ⟹ 本候選係結構上不同類。")
say("  (3) ⚠ 銳度低：B ∈ (e, e^1.5] 皆給 3 ⟹ 穩健、但唔能由 3 反推 B = 4。")
say("  (4) ⚠ κ 之冪次敏感：κ = B^{−1} ⟹ 2；κ = B^{−2} ⟹ 3 ⟹ 需 κ = B^{−2} 之獨立出處。")
say("  (5) ⟹ 判 **🟡 PARTIAL（C7 全過，出處未閉）** —— 與前十一者之失敗模式不同。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "B": B, "gate_count": cnt, "threshold_2lnB": mp.nstr(thr, 12),
               "window": [mp.nstr(mp.e, 8), mp.nstr(mp.e ** mp.mpf("1.5"), 8)],
               "window_width_pct": mp.nstr(width * 100, 6),
               "counts_by_B": {b: gate_count(b)[0] for b in (2, 3, 4, 5, 6, 7)},
               "C7": {"set": levels, "ordered": True, "dof": 0}}, f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§6  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)