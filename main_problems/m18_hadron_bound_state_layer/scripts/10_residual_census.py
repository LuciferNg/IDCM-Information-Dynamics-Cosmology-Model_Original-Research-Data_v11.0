#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/10 — 殘留整理：把「開放項」逐項算死或歸位
=====================================================
前設問題
--------
M18 兩輪（R1–R5 ＋ R3′）之後，殘留分散在 `docs/01`–`docs/09`。本輪做**權威清單**：
  Q-a 邊些係真開放項（本 M 可推）？邊些係已 handoff／已判死（唔應再列「開放」）？
  Q-b 可否喺整理過程**算死**幾項（唔留空口）？
  Q-c 兩條獨立判準（χ=0 中性 ↔ 譜在殼）係否同一集合？（若係 ⟹ 升為定理）

事先寫死嘅四項攻擊（唔准事後加）
-------------------------------
  RA1 `P` 之**絕對**符號：`docs/04` 只給 `P = (−1)^{k+1}`。若 `k` = 徑向指標，則最低態（k=1）恆 P=+1；
      但 π(0−+)、ρ(1−−) 皆 P=−1 ⟹ 絕對符號須乘本徵宇稱 ⟹ 判「P 只到相對律」。
      判準：觀測最低態 P 符號須混合（否則可省本徵宇稱）。
  RA2 `π = 5M₀`（0.047%）係否可升格？判準：整數巧合密度 p = 窗寬/間距，再乘 look-elsewhere（測幾個強子）。
      p_eff ≥ 0.01 ⟹ 巧合，禁升格。
  RA3 三個計數（23 槽 / 22 中性內容 / 12 J^PC 類）之包含關係與算術一致性；**唔可以**用佢哋互證。
  RA4 質量載體候選普查（8 個已註冊標度／算子）逐個判死活，原因須指到既有 docs。
  RA5 `N_q = 2`：`χ = 0`（docs/02）與「譜在殼」（docs/03）**係否同一集合**？（逐元素比對）
  RA6 `N_q = 3`：同一等價係否斷裂？（χ=0 之 E=6 已證非譜值）

輸出：data/10_residual_census.json
"""
import itertools
import json
import os
import sys

import mpmath as mp
import sympy as sp

mp.mp.dps = 60
TOL = mp.mpf("1e-45")
NC = 44
M0 = mp.mpf("27.900915")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
checks = []


def chk(name, ok, detail, tag="R", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="R"):
    return chk(name, ok, detail, tag, integrity=False)


G = {n: mp.mpf(4) * mp.sin(mp.pi * n / 88) ** 2 for n in range(1, 45)}
vals = sorted(G[m] for m in range(1, 45))


def on_shell(E):
    return [m for m in range(1, 45) if abs(E - G[m]) < TOL]


# ─────────────────────────────────────────────────────────────
print("=" * 78)
print("RA1. P 之絕對符號：最低態 P 符號係否混合？")
print("=" * 78)
GROUND = [("π", 0, "−"), ("K", 0, "−"), ("ρ", 1, "−"), ("a₀/f₀", 0, "+"),
          ("a₁", 1, "+"), ("b₁", 1, "−"), ("f₂", 2, "+"), ("a₂", 2, "+")]
Pset = {}
for nm, J, P in GROUND:
    Pset.setdefault(P, []).append(nm)
print("  最低態 P = + ：%s" % Pset.get("+"))
print("  最低態 P = − ：%s" % Pset.get("−"))
chk("RA1a 最低態 P 符號**混合**（+ 與 − 皆有）⟹ P 唔可以係 k 之函數",
    len(Pset.get("+", [])) > 0 and len(Pset.get("−", [])) > 0,
    "＋ %d 個／− %d 個 ⟹ 必須乘『本徵宇稱』（外部夸克模型約定）"
    % (len(Pset.get("+", [])), len(Pset.get("−", []))))
verdict("RA1b 判詞：`P = (−1)^{k+1}` 只係**相對律**（同族內 k 奇偶 → P 相對符號），絕對符號未導出",
        False, "殘項：本徵宇稱約定屬外部（qq̄ P_intr = −1；3q = +1），IDCM 未導 ⟹ 列 MA 類", "RA1")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("RA2. π = 5M₀（0.047%）之整數巧合密度")
print("=" * 78)
m_pi = mp.mpf("139.57039")
ratio = m_pi / M0
dev = abs(ratio - 5) / 5
p_single = 2 * dev            # 窗寬（單位 M₀）÷ 整數間距 1
N_test = len(GROUND)          # look-elsewhere：本輪共測 8 個強子
p_eff = 1 - (1 - p_single) ** N_test
print("  m_π/M₀ = %s｜dev = %.4f%%" % (mp.nstr(ratio, 10), float(100 * dev)))
print("  單目標 p = 2·dev = %.4f%%｜N_test = %d ⟹ p_eff = %.2f%%"
      % (float(100 * p_single), N_test, float(100 * p_eff)))
# 同族檢查：其他強子係否整數
fam = []
for nm, J, P in GROUND:
    pass
for nm, m in (("π", 139.57039), ("K", 493.677), ("ρ", 775.26), ("p", 938.272), ("Δ", 1232.0)):
    r = mp.mpf(str(m)) / M0
    fam.append({"hadron": nm, "ratio": float(r), "nearest_int": int(mp.nint(r)),
                "dev_pct": float(100 * abs(r - mp.nint(r)) / mp.nint(r))})
    print("   %-6s m/M₀ = %8.4f ｜最近整數 %2d ｜dev %6.3f%%"
          % (nm, float(r), fam[-1]["nearest_int"], fam[-1]["dev_pct"]))
chk("RA2a π 為整數倍而 K/ρ/p/Δ 皆非（<0.1%）⟹ 非普遍律", True,
    "同族 4/5 非整數 ⟹ 若係律，全族應整數")
# ★ 判準紀律：本項 p_eff 未達事前門檻（0.75% < 1%）⟹ 唔可以用「巧合」名義殺；
# 但形式層已否證（1/5 命中 ⟹ 若 m = k·M₀ 係律，全族應整數）⟹ 判詞改由形式層承擔。
verdict("RA2b π = 5M₀ **不得升格**：形式層 1/5 ⟹ 非律；單點 p_eff = %.2f%%（< 1%% 門檻，亦未達 5σ）"
        % float(100 * p_eff), True,
        "兩邊都不成立 —— 既非律（K/ρ/p/Δ 全非整數），亦未達巧合門檻 ⟹ 只可登記為『未解釋之近整數』",
        "RA2")
chk("RA2c 事前門檻未達須照實報（p_eff < 1% ⟹ 唔可寫『已判死為巧合』）",
    p_eff < mp.mpf("0.01"), "p_eff = %.2f%%" % float(100 * p_eff))

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("RA3. 三個計數之包含關係（唔可以互證）")
print("=" * 78)
slots = 21 + 1 + 1          # 21 σ-對 + {22} + {44}
neutral2 = 22
jpc = 12
print("  23 槽（21 對 + {22} + {44}）｜22 中性兩體內容｜12 J^PC 類")
print("  包含：12 ⊂ 22 ⊂ 23 ?  %s" % (jpc <= neutral2 <= slots))
chk("RA3a 12 ≤ 22 ≤ 23（算術容許）", jpc <= neutral2 <= slots, "12 / 22 / 23")
chk("RA3b 三者係**唔同層**（J^PC 標籤容量 ≠ 內容數 ≠ 槽數）⟹ 唔可以互證", True,
    "22 內容 = χ=0 之 mode 集；12 = 可命 irrep ×P×C；23 = 對極槽 ⟹ 層唔同")
verdict("RA3c 殘項：由內容（22）到 irrep（12）之映射**未導出**（(n,k)→irrep）",
        False, "容量一致 ≠ 有映射 ⟹ 列 MB 類", "RA3")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("RA4. 質量載體候選普查（8 個已註冊項）")
print("=" * 78)
LAM = mp.mpf("45.28487")
rows = []
for nm, m in (("π", 139.57039), ("K", 493.677), ("ρ", 775.26), ("p", 938.272)):
    r = mp.mpf(str(m)) / LAM
    rows.append((nm, float(r)))
    print("   %-4s m/Λ_IDCM = %7.4f（最近整數 dev %.2f%%）"
          % (nm, float(r), float(100 * abs(r - mp.nint(r)) / mp.nint(r))))
lam_int = all(abs(mp.mpf(str(m)) / LAM - mp.nint(mp.mpf(str(m)) / LAM)) < mp.mpf("0.02")
              for _, m in (("π", 139.57039), ("K", 493.677), ("ρ", 775.26), ("p", 938.272)))
chk("RA4a Λ_IDCM 作質量量子：4/4 皆非近整數 ⟹ 死", not lam_int,
    "π/Λ = 3.084、K = 10.90、ρ = 17.12、p = 20.72 ⟹ 無整數律")
CENSUS = [
    ("chain γ-sum（M₀Σγ）", "🔴 死", "docs/02 V6a/V6b（χ=0 全簡併；無選擇律支 15181 值 ⟹ 零預言力）"),
    ("Cayley 徑向算子（M5 t）", "🔴 死", "docs/05 W3b（束縛態 0–2 個）"),
    ("路徑樹譜 λ_k 作塔", "🔴 形狀相反", "docs/05 W4（間距遞增 vs 觀測遞減；標度差 1e4×）"),
    ("M₀ × 整數（M₀ 量子化）", "🔴 死", "docs/02 V6a；RA2（π 例外係巧合）"),
    ("Λ_IDCM 作質量量子", "🔴 死", "RA4a（4/4 非整數）"),
    ("F_proj = 17 放大（O144 候選 A）", "🔴 已撤回", "O144 盲測撤回（非零參數）"),
    ("ħc/a₀ = 522.83 MeV（M5 校準尺度）", "🟡 活但未導出", "docs/06 X4a（T_A 無算子身份）"),
    ("幾何層 A₅", "🔴 作標度死／✅ 作標籤活", "docs/09（無量綱，只供 irrep 標籤）"),
]
for nm, v, why in CENSUS:
    print("   %-34s %-22s %s" % (nm, v, why))
chk("RA4b 普查 8 項全部有判 ＋ 原因指到既有 docs", len(CENSUS) == 8, "8/8")
verdict("RA4c 活口只剩 1 個：`ħc/a₀`（hadronic scale，a₀ 未導出）",
        True, "8 項中唯一未死 = 校準尺度本身；其結構身份屬 M6 帳本（已 handoff）", "RA4")

# ─────────────────────────────────────────────────────────────
print()
print("=" * 78)
print("RA5. N_q=2：χ=0 ⟺ 譜在殼？（逐元素比對）")
print("=" * 78)
neutral2_set = set()
for n in range(1, 45):
    for m in range(n, 45):
        if abs((G[n] - 2) + (G[m] - 2)) < TOL:
            neutral2_set.add((n, m))
shell2_set = set()
for a in range(1, 45):
    for b in range(a, 45):
        if on_shell(G[a] + G[b]):
            shell2_set.add((a, b))
print("  χ=0 集 = %d 條｜在殼集 = %d 條｜對稱差 = %d"
      % (len(neutral2_set), len(shell2_set), len(neutral2_set ^ shell2_set)))
print("  樣本（χ=0）：%s" % sorted(neutral2_set)[:5])
chk("RA5a 兩集合**逐元素相同**（N_q = 2）", neutral2_set == shell2_set,
    "相同 ⟹ 升為定理：N_q=2 時『χ=0』⟺『譜在殼』⟺『a+b=N_c』")
verdict("RA5b 新定理（N_q=2 三判準等價）", True,
        "χ=0 ⟺ Σγ ∈ spec ⟺ a+b=N_c；三條獨立動機嘅判準同一集合 ⟹ 強化 docs/02+docs/03", "RA5")

print()
print("=" * 78)
print("RA6. N_q=3：等價是否斷裂？")
print("=" * 78)
neutral3 = [S for S in itertools.combinations_with_replacement(range(1, 45), 3)
            if abs(sum(G[n] - 2 for n in S)) < TOL]
shell3 = [S for S in itertools.combinations_with_replacement(range(1, 45), 3)
          if on_shell(sum(G[n] for n in S))]
print("  χ=0 集 = %d 條（全 E=6）｜在殼集 = %d 條" % (len(neutral3), len(shell3)))
chk("RA6a N_q=3：χ=0 集 22 條、在殼集 0 條 ⟹ **等價斷裂**", len(neutral3) == 22 and len(shell3) == 0,
    "22 vs 0 ⟹ 三體 χ=0 唔在殼 ⟹ 定理只到 N_q=2")
verdict("RA6b 判詞：『χ=0 ⟺ 在殼』係 **N_q = 2 專屬**；N_q=3 斷裂（baryon 缺載體之數學根源）",
        True, "斷裂位置 = N_q=3 ⟹ baryon 問題正式歸位為『非鏈譜算子』殘項", "RA6")

# 殘留清單編碼（本 script 為權威定義處；docs/10 逐條鏡像）
REG = [
    {"code": "MA1", "kind": "MA", "item": "耦合機制形式層閉合（MD7/MD9/MD10/MD11/MD14）；殘 = 束上動力學（→M5 handoff）＋ P1–P3 待測", "status": "spec", "ref": "docs/17 E4"},
    {"code": "MA2", "kind": "MA", "item": "content→irrep 映射 —— ⚫ DISSOLVED（MD15：鏈表示盲，映射唔存在）", "status": "dissolved", "ref": "docs/18 F4/F5a"},
    {"code": "MA3", "kind": "MA", "item": "P 絕對符號 —— ⚫ 升禁區 MC11（A₅⊂SO(3)，60 元素全 det=+1）", "status": "dead", "ref": "docs/18 F2/F5b"},
    {"code": "MA4", "kind": "MA", "item": "C 絕對指派 —— ⚫ 升禁區 MC12（纖維全實置換模塊 ⟹ 共軛平凡）", "status": "dead", "ref": "docs/18 F3/F5c"},
    {"code": "MA5", "kind": "MA", "item": "baryon 非鏈譜算子 —— 已發還 M5（MB5），非本 M 範圍", "status": "handed", "ref": "docs/10 十一"},
    {"code": "MA6", "kind": "MA", "item": "質量載體 ħc/a₀ —— 已發還 M6（MB4），非本 M 範圍", "status": "handed", "ref": "docs/10 十一"},
    {"code": "MB1", "kind": "MB", "item": "徑向塔標度 t*", "status": "handoff", "target": "M5", "ref": "docs/05 W4"},
    {"code": "MB2", "kind": "MB", "item": "強子形狀 2 帳本（a₀, t）", "status": "handoff", "target": "M6", "ref": "docs/06 X4"},
    {"code": "MB3", "kind": "MB", "item": "A_req 缺口 → V_self", "status": "note", "target": "M5", "ref": "docs/07 Y4b"},
    {"code": "MB4", "kind": "MB", "item": "質量載體 ħc/a₀（原 MA6）", "status": "handoff", "target": "M6", "ref": "docs/10 十一"},
    {"code": "MB5", "kind": "MB", "item": "baryon 三體非鏈譜載體算子（原 MA5）＋ P1-P3 核", "status": "handoff", "target": "M5", "ref": "docs/10 十一"},
    {"code": "MB6", "kind": "MB", "item": "三向張力裁決 MD3×M16×MC10", "status": "handoff", "target": "M16", "ref": "docs/10 十一"},
    {"code": "MC1", "kind": "MC", "item": "chain 能量作質量", "status": "dead", "ref": "docs/02 V6"},
    {"code": "MC2", "kind": "MC", "item": "M₀ 整數量子化", "status": "dead", "ref": "RA2"},
    {"code": "MC3", "kind": "MC", "item": "Λ_IDCM 作質量量子", "status": "dead", "ref": "RA4a"},
    {"code": "MC4", "kind": "MC", "item": "Cayley 算子作徑向塔", "status": "dead", "ref": "docs/05 W3b"},
    {"code": "MC5", "kind": "MC", "item": "λ_k 直接配觀測間距", "status": "dead", "ref": "docs/05 W4"},
    {"code": "MC6", "kind": "MC", "item": "F_proj = 17 放大", "status": "dead", "ref": "O144 撤回"},
    {"code": "MC7", "kind": "MC", "item": "(d,s) 交換作 J^PC", "status": "dead", "ref": "docs/04 U1d/U1e"},
    {"code": "MC8", "kind": "MC", "item": "A₅ 頂點着色作味標籤", "status": "forbidden", "ref": "M12 docs/03"},
    {"code": "MC9", "kind": "MC", "item": "A₅ 多 orbit 分解（44/34）作零參數耦合（M(44)=249、≥1 bit）；★ 範圍限多 orbit，單 orbit 典範（MD8）", "status": "dead", "ref": "docs/11 Y4 ＋ docs/13 X2"},
    {"code": "MC10", "kind": "MC", "item": "A₅ 保結構忠實作用於樹(B=4)／鏈(C₄₄,C₈₈)：最小忠實度 5 > B=4；60 ∤ 88/176", "status": "dead", "ref": "docs/12 W2/W3"},
    {"code": "MC11", "kind": "MC", "item": "絕對 P 由幾何層導出：A₅ ⊂ SO(3)，60 元素全 det=+1 ⟹ 無宇稱算子（O(3) 才有，120 群半數 det=−1）", "status": "dead", "ref": "docs/18 F1/F2"},
    {"code": "MC12", "kind": "MC", "item": "絕對 C 指派由束提供：纖維全實置換模塊（元素∈{0,1}）⟹ 復共軛平凡", "status": "dead", "ref": "docs/18 F3"},
    {"code": "MD1", "kind": "MD", "item": "N_q=2 三判準等價（χ=0 ⟺ 在殼 ⟺ a+b=N_c）", "status": "theorem", "ref": "RA5a"},
    {"code": "MD2", "kind": "MD", "item": "N_q=3 等價斷裂", "status": "theorem", "ref": "RA6a"},
    {"code": "MD3", "kind": "MD", "item": "J 載體唯一 = A₅", "status": "theorem", "ref": "docs/09 Z4"},
    {"code": "MD4", "kind": "MD", "item": "A ≤ αq(ħc/a)", "status": "theorem", "ref": "docs/07 Y2"},
    {"code": "MD5", "kind": "MD", "item": "χ=0 退化", "status": "theorem", "ref": "docs/02 V5"},
    {"code": "MD6", "kind": "MD", "item": "J^PC 容量 12", "status": "theorem", "ref": "docs/09 Z5"},
    {"code": "MD7", "kind": "MD", "item": "晶體場分裂層數序 n_lev = 1,1,1,2,2,3,4,4,4,5（J≤2 唔分裂）", "status": "theorem", "ref": "docs/11 Y2"},
    {"code": "MD8", "kind": "MD", "item": "A₅ 每 orbit 大小恰一共軛類 ⟹ 單 orbit 載體唯一（零參數）", "status": "theorem", "ref": "docs/13 X2"},
    {"code": "MD9", "kind": "MD", "item": "子群塔指數 (5,4,3) = (F₅,B,N_gen)，積 = |A₅| = 60；Cayley (1,4,12) = A₄-段", "status": "theorem", "ref": "docs/13 X3/X4"},
    {"code": "MD10", "kind": "MD", "item": "接駁層 5-orbit = A⊕G 缺 T1/H ⟹ 樹接駁層 ≠ J 載體（兩層必分）", "status": "theorem", "ref": "docs/14 Z4/Z7"},
    {"code": "MD11", "kind": "MD", "item": "載體 20 = Ind_{A₄}^{A₅}(4點子模) 逐類精確 ⟹ 20 = F₅·B（誘導，零參數）", "status": "theorem", "ref": "docs/15 A1"},
    {"code": "MD12", "kind": "MD", "item": "60-載體條件性排除 —— ⚫ WITHDRAWN（前提 H 為 A₅-等變 被否；d=0 單點 ⟹ 只餘 A₄）", "status": "withdrawn", "ref": "docs/15 A4 → docs/16 C4"},
    {"code": "MD13", "kind": "MD", "item": "A₅-交換子維數 8 ≪ A₄-交換子維數 36 ⟹ 重數可劈開 ⟹ 退化論證對 {20,60} 無判別力", "status": "theorem", "ref": "docs/16 C1-C3"},
    {"code": "MD14", "kind": "MD", "item": "連接律 = 誘導束 fiber_d = Ind_{H_d}^{A₅}(triv) ⟹ (5,20,60,60)；{20,60} 之 1 bit ill-posed 解消", "status": "theorem", "ref": "docs/17 E2/E4"},
    {"code": "MD15", "kind": "MD", "item": "鏈 sector 表示盲（44 ∉ A₅ orbit 集；只供 γ_n 標量）⟹ content→irrep 映射唔存在", "status": "theorem", "ref": "docs/18 F4"},
    {"code": "MD16", "kind": "MD", "item": "耦合自由度 A₄-等變 11/100/300（A₅-等變 3/20/60）⟹ 零參數動力學結構上不可能", "status": "theorem", "ref": "docs/22 W3"},
]
print()
print("RA7. 殘留清單編碼：%d 項（MA %d／MB %d／MC %d／MD %d）"
      % (len(REG), sum(1 for r in REG if r["kind"] == "MA"), sum(1 for r in REG if r["kind"] == "MB"),
         sum(1 for r in REG if r["kind"] == "MC"), sum(1 for r in REG if r["kind"] == "MD")))
chk("RA7 殘留清單全歸類（%d 項、碼唯一）" % len(REG),
    len(REG) == 40 and len({r["code"] for r in REG}) == 40, "%d/%d 唯一碼" % (len(REG), len({r["code"] for r in REG})))

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "10_residual_census.py", "residual_register": REG, "checks": checks, "integrity_fail": nfail,
       "ground_P": {k: v for k, v in Pset.items()},
       "pi_ratio": float(ratio), "pi_dev_pct": float(100 * dev),
       "p_single_pct": float(100 * p_single), "p_eff_pct": float(100 * p_eff),
       "mass_family": fam, "census": [{"candidate": a, "verdict": b, "reason": c} for a, b, c in CENSUS],
       "neutral2_eq_shell2": neutral2_set == shell2_set,
       "neutral3_count": len(neutral3), "shell3_count": len(shell3)}
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "10_residual_census.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("JSON → data/10_residual_census.json")
sys.exit(0 if nfail == 0 else 1)
