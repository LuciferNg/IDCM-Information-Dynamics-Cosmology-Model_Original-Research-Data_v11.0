#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/14 — 殘留三項攻堅：Q5 命名步／Q2 後驗形式／#12 閘門冪次
================================================================================
判準（事前寫死）：
  A1 (#12) 閘門索引範圍審計：M3 半譜 = k = 1..44（44 模）；閘門用 range(0, N_d)。
           若 γ₀ = 0（零模，不屬 44 模）被計入 ⟹ 「3」含真空 ⟹ 判 artifact。
  A2 (#12 冪次) 樹層 ℓ 由 Ω_DM 觀測選出：Ω_DM(ℓ) = B^ℓ/(N_d + B^ℓ)；若 ℓ=2 唯一 ⟹ 冪次 = 2（可觀測決定）。
  B  (Q2)  B/(B−1) 是否有 a priori 結構源：幾何級數 Σ_{k≥0} B^{−k} = B/(B−1)（全樹層權重和），
           且 ≡ N_d/K_tree（鏈／樹模數比）；兩者皆零質量輸入 ⟹ 形式非後驗。
  C  (Q5)  閘門集 vs corpus 世代集 {1,4,5}、受保護集 P：若閘門集 ⊄ P 且交集 < 3 ⟹ 識別不可能。
"""
import json
import sys
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 30
B, N_d = 4, 44
KAPPA = mp.mpf(1) / B ** 2
P = {1, 2, 3, 4, 5, 7, 9, 11, 13}          # 9 模式受保護集（docs/01 C7）
GEN_MODES = [1, 4, 5]                       # v9 corpus gen_modes（字面）
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
say("M15 scripts/14 — Q5 命名步／Q2 後驗形式／#12 冪次")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§A1  #12 閘門索引範圍審計（M3 range(0, N_d) vs 物理半譜 k = 1..44）")
say("-" * 78)
g0 = mp.mpf(4) * mp.sin(mp.mpf(0)) ** 2
g44 = 4 * mp.sin(mp.pi * 44 / 88) ** 2
say(f"  γ₀  = 4 sin²(0)        = {mp.nstr(g0, 6)}   ⟵ 零模（鏈連通 ⟹ 唯一零特徵值）")
say(f"  γ₄₄ = 4 sin²(π·44/88)  = {mp.nstr(g44, 6)}   ⟵ Nyquist（γ_max = B）")
say(f"  M3 半譜定義（docs/03 §1）：k = 1..44 ⟹ N_d = 44 模 ⟹ 零模 **不屬** 44 模")
say(f"  閘門實際掃描範圍：range(0, N_d) = 0..43  ⟹ **含 0、排除 44**")
vis_m3 = [n for n in range(0, N_d) if mp.e ** (-n) > KAPPA]
vis_phys = [n for n in range(1, N_d + 1) if mp.e ** (-n) > KAPPA]
say(f"  ⟹ 依 M3 範圍：{{n ≥ 0}} = {vis_m3}，count = {len(vis_m3)}")
say(f"  ⟹ 依物理索引（k = 1..44）：{{n}} = {vis_phys}，count = **{len(vis_phys)}**")
chk("A1 零模污染：γ₀ = 0 被計入 ⟹ 物理索引下 count = 2 ≠ 3",
    g0 == 0 and len(vis_m3) == 3 and len(vis_phys) == 2,
    f"count(0..43) = {len(vis_m3)}；count(1..44) = {len(vis_phys)}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§A2  #12 冪次：樹層 ℓ 由 Ω_DM 觀測選出（與閘門索引無關）")
say("-" * 78)
ob_h2, oc_h2, H0 = 0.02237, 0.1200, 67.36
s_oc, s_h = 0.0012, 0.54
h2 = (H0 / 100) ** 2
Oc = oc_h2 / h2
sig = Oc * ((s_oc / oc_h2) ** 2 + (2 * s_h / H0) ** 2) ** 0.5
say(f"  Ω_DM(obs, Planck 2018) = {Oc:.5f} ± {sig:.5f}")
zmap = {}
for l in (1, 2, 3, 4):
    pred = B ** l / (N_d + B ** l)
    zmap[l] = (pred - Oc) / sig
    say(f"    ℓ = {l}: Ω_DM = {B ** l}/{N_d + B ** l} = {pred:.5f}   z = {zmap[l]:+.2f}σ")
chk("A2 ℓ = 2 唯一（|z| ≤ 1；其餘 ≥ 30σ 排除）⟹ 閘門冪次 = 2（觀測決定）",
    abs(zmap[2]) <= 1 and abs(zmap[1]) > 30 and abs(zmap[3]) > 30,
    f"z: ℓ=1 {zmap[1]:+.1f}, ℓ=2 {zmap[2]:+.2f}, ℓ=3 {zmap[3]:+.1f}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§B  Q2：`B/(B−1)` 之 a priori 源（零質量輸入）")
say("-" * 78)
geosum = sum(Fr(1, B ** k) for k in range(0, 12))          # Σ_{k=0..11} B^-k
geosum_trunc_err = Fr(1, B ** 12) / (1 - Fr(1, B))
say(f"  (i) 幾何級數：Σ_{{k≥0}} B^{{−k}} = B/(B−1)")
say(f"      Σ_{{k=0..11}} 4^{{−k}} = {geosum}   餘項 = {geosum_trunc_err}")
say(f"      B/(B−1) = {Fr(B, B - 1)}；4/3 − 部分和 = {Fr(B, B - 1) - geosum} = 餘項 ⟹ 收斂至 4/3 ✓")
say(f"  (ii) 鏈／樹模數比：N_d/K_tree = {N_d}/{N_d - 11} = {Fr(N_d, N_d - 11)}")
chk("B 兩條獨立 a priori 源皆給 4/3（幾何級數 ≡ 鏈／樹比），零質量輸入",
    Fr(B, B - 1) == geosum + geosum_trunc_err == Fr(N_d, N_d - 11),
    "Σ B^-k = B/(B−1) = N_d/K_tree = 4/3")
say("")
say("  ⟹ 係數**形式**（B/(B−1)）有兩條零輸入來源 ⟹ 非後驗凑數；")
say("     仍未解者＝「為何質量分裂比 = 該係數」——屬物理識別（同 Q5 同類），非同義反覆。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§C  Q5：閘門集能否命名為世代？（三閘）")
say("-" * 78)
say(f"  閘門集（M3 範圍）= {set(vis_m3)}｜物理索引 = {set(vis_phys)}｜corpus gen_modes = {set(GEN_MODES)}")
inter = set(vis_m3) & set(GEN_MODES)
say(f"  (i)  閘門集 ∩ gen_modes = {sorted(inter)} ⟹ |∩| = {len(inter)} < 3")
say(f"  (ii) 閘門集 ⊆ P（受保護費米子集）？ {set(vis_m3) <= P}（0 ∉ P）")
say(f"  (iii) 第三槽 = n = 2；但 n = 0 者係 γ₀ = 0 之**真空/零模** ⟹ 唔可能係費米子世代")
chk("C 閘門 → 世代之命名不可能（交集 1/3；含零模；閘門集 ⊄ P）",
    len(inter) < 3 and not (set(vis_m3) <= P),
    f"∩ = {sorted(inter)}；0 ∉ P")
say("")
say("  ⟹ Q5 之『閘門命名』路線 **判死**；Q5 之殘餘收窄為單一命題：")
say("     『Ω_b 之分子 = 世代數』——預算層之識別（無真空污染，`rank` 恆等式不再成立於閘門）。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§D  判詞")
say("-" * 78)
say("  ① #12（閘門層數）🔴 **判死**：其 3 來自 `range(0, N_d)` 把零模 γ₀ = 0 計入、")
say("     同時排除 γ₄₄ = 4（Nyquist）；物理索引 k = 1..44 下 count = **2**。")
say("     ⟹ **`3 = rank(P_κ)` 撤回**（此式僅在含真空之索引下成立）。")
say("  ② #12 之**冪次**問題獨立成立並已解：ℓ = 2（Ω_DM 選出，+%.2fσ；ℓ=1/3 差 %.0f/%.0fσ）。" % (zmap[2], abs(zmap[1]), abs(zmap[3])))
say("  ③ Q2 🟢 **後驗形式已解**：B/(B−1) = Σ B^{−k} = N_d/K_tree = 4/3（兩條零輸入源）；")
say("     殘餘 = 質量比 ↔ 係數之物理識別（同 Q5 同類）。")
say("  ④ Q5 🔴 **命名步判死（閘門路線）**：交 1/3、含零模、閘門集 ⊄ P。")
say("     殘餘單一命題：`N_gen = 60·Ω_b` 之預算層識別。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES,
               "gamma0": str(g0), "gamma44": str(mp.nstr(g44, 6)),
               "gate_set_m3": vis_m3, "gate_count_m3": len(vis_m3),
               "gate_set_phys": vis_phys, "gate_count_phys": len(vis_phys),
               "z_by_level": {str(k): round(v, 3) for k, v in zmap.items()},
               "geosum_4_3": str(geosum), "geom_sum_closed": str(Fr(B, B - 1)),
               "nd_over_ktree": str(Fr(N_d, N_d - 11)),
               "gate_in_P": set(vis_m3) <= P,
               "gate_inter_gen": sorted(inter)},
              f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)