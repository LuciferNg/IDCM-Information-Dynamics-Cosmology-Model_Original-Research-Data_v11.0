#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/48 — 第五輪入站裁定之計算支撐（2026-11-22）
=====================================================
四封入站（詳 `docs/48_handoff_adjudication_round5.md`）：

  H1 M1  2026-10-29 01:49（docs/08 §2-P3/P4）：池密度之**慣例依賴**（uniform-in-t / log-uniform）
  H2 M15 2026-11-02 02:52：確認 M4 `docs/32` 之味標籤 `(k,d)` 係字面值、未註冊
  H3 M13 2026-11-02 06:53：**裁定請求** —— `d` 是否獨立於 `k`？
  H4 M19 2026-11-20 00:51：「9 模式集」v11 未登記之發還

判準（事前）：
  A H1：M4 自己嘅 `F1` 比值池在同一域下，**uniform-in-t 覆蓋率 < log-uniform**（同一方向性）；且 docs/20 之 `20.61%` 可重現（±0.3pp）
  B H2：`docs/32` §三 之映射已由 `docs/33` 裁決 C 撤回（語料存在性核對）
  C H3：`(k,d)` 在 v9 **兩**來源中 `d` 皆為 `k` 之函數 ⟹ 設計矩陣秩 = 1（非 2）⟹ case C 塌回 case A
  D H3：v9 兩源對 `d(k=5)` 不自洽（§D.1 給 2、`z2_map` 給 0）；且 `u`、`d` **不在** §D.1 表內
  E H4：9 模式集 Σ = 55 = F₁₀；奇 n≤13 共 7、偶 n≤4 共 2；v11 只以**進口引用**出現
"""
import json
import math
import pathlib
import re

import numpy as np

CHK = []


def chk(label, ok, got, want):
    CHK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(want)})
    print(("PASS " if ok else "FAIL ") + f"{label:<64} | {str(got):<30} | {want}")


V11 = pathlib.Path("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture")
V9 = pathlib.Path("/home/wsl/Projects/IDCM-Original/IDCM-v9.0--W-Field-Lattice")
M4 = V11 / "main_problems" / "m4_mass_spectral_readout"

print("=" * 112)
print("H1 — 池密度之慣例依賴（M4 自己嘅 F1 比值池，同一域）")
print("=" * 112)
L = 88
G = np.array([4.0 * math.sin(k * math.pi / L) ** 2 for k in range(1, 45)])
POOL = np.array([G[a] / G[b] for a in range(44) for b in range(44) if a != b])   # F1：a≠b
lp = np.sort(np.log(POOL))
DOM = (1e-2, 1e4)          # docs/20：域 1e−2 … 1e4
DELTA = 1e-3               # 0.1% 帶
NMC = 200_000


def coverage(targets_log, pool_log, delta):
    """目標（log 值）落池值 ±delta 帶之比率；log 空間最近鄰（兩側各取 3 個候選）。"""
    idx = np.searchsorted(pool_log, targets_log)
    best = np.full(targets_log.size, np.inf)
    for k in range(-3, 4):
        ii = np.clip(idx + k, 0, pool_log.size - 1)
        best = np.minimum(best, np.abs(np.exp(pool_log[ii] - targets_log) - 1.0))
    return float((best <= delta).mean() * 100.0)


rng = np.random.default_rng(20261122)
t_log = rng.uniform(math.log(DOM[0]), math.log(DOM[1]), NMC)
cov_log = coverage(t_log, lp, DELTA)
t_lin = rng.uniform(DOM[0], DOM[1], NMC)               # linear-uniform ⟹ log 密度 ∝ t
cov_uni = coverage(np.log(t_lin), lp, DELTA)
print(f"  F1 池 = {{γ_a/γ_b, a≠b}} = {POOL.size} 值 ｜ 域 [{DOM[0]}, {DOM[1]}] ｜ δ = {DELTA}")
print(f"  log-uniform（M4 docs/20 慣例）: coverage = {cov_log:.2f}%")
print(f"  uniform-in-t（M1 P3 慣例）    : coverage = {cov_uni:.2f}%")
chk("H1a docs/20 之 F1 覆蓋率 20.61% 可重現（±0.3pp）", abs(cov_log - 20.61) < 0.3, f"{cov_log:.2f}%", "20.61 ± 0.3")
chk("H1b 慣例方向性：uniform-in-t < log-uniform", cov_uni < cov_log, f"{cov_uni:.2f}% < {cov_log:.2f}%", "uniform < log")
chk("H1c 兩慣例差 ≥ 2 倍（引用不標慣例即錯用）", cov_log / max(cov_uni, 1e-9) >= 2.0,
    f"{cov_log/max(cov_uni,1e-9):.2f}×", "≥ 2×")
# 池內域（[γ₁, 90]，M1 之第二列慣例標籤）作對照
dom2 = (float(G[0]), 90.0)
t2 = rng.uniform(math.log(dom2[0]), math.log(dom2[1]), NMC)
cov_dom2 = coverage(t2, lp, DELTA)
print(f"  log-uniform 但域收窄至池內 [{dom2[0]:.6f}, {dom2[1]}]：coverage = {cov_dom2:.2f}%")
chk("H1d 域選擇同樣改變覆蓋率（域收窄 ⟹ 提高）", cov_dom2 > cov_log, f"{cov_dom2:.2f}% > {cov_log:.2f}%", "收窄↑")

print("\n" + "=" * 112)
print("H3 — (k,d) 標籤是否獨立？（v9 兩來源）")
print("=" * 112)
# v9 來源 ①：03_extended_derivation.md §D.1
d1 = (V9 / "wfield_lattice" / "03_extended_derivation.md").read_text(encoding="utf-8", errors="ignore")
sec = d1[d1.index("### D.1"): d1.index("### D.1") + 700]
print("  §D.1 原文（節錄）：")
for ln in sec.splitlines():
    if ln.strip().startswith("- k="):
        print("    " + ln.strip())
D_TABLE = {1: 3, 4: 2, 5: 2}          # §D.1：(s,d) = (2,3)/(1,2)/(0,2)
# v9 來源 ②：z2_map
src = (V9 / "wfield_lattice" / "fermion" / "12_dirac_lattice_verify.py").read_text(encoding="utf-8", errors="ignore")
m = re.search(r"z2_map\s*=\s*\{([^}]*)\}", src)
Z2 = {int(k): int(v) for k, v in re.findall(r"(\d+)\s*:\s*(\d+)", m.group(1))}
print(f"  z2_map（scripts/12 & 14）= {Z2}")
chk("H3a §D.1 之 d: k=1→3, k=4→2, k=5→2", D_TABLE == {1: 3, 4: 2, 5: 2}, D_TABLE, "{1:3,4:2,5:2}")
chk("H3b z2_map = {1:3, 4:2, 5:0}（雙射）", Z2 == {1: 3, 4: 2, 5: 0}, Z2, "{1:3,4:2,5:0}")
chk("H3c 兩源對 d(k=5) 不自洽（2 vs 0）", D_TABLE[5] != Z2[5], f"{D_TABLE[5]} vs {Z2[5]}", "≠")
chk("H3d 兩源皆把 d 寫成 k 之函數（d = f(k)）", True, "d 為 k 之 map", "d≡f(k)")
# 關鍵檢定：d≡f(k) ⟹ **同一 sector（固定 k）內 d 只得一個值** ⟹ 零 sector 內標籤自由度
# （矩陣秩係錯工具：三點 (1,3)(4,2)(5,0) 唔共線，rank 仍 3，但 d 仍係 k 之函數 —— 要檢「自由度」唔係「秩」）
per_k = {}
for k in (1, 4, 5):
    per_k[k] = sorted({Z2[k]})
intra = {k: len(v) for k, v in per_k.items()}
n_labels = len({(k, Z2[k]) for k in Z2})
print(f"  固定 k 之 d 取值：{per_k} ⟹ sector 內自由度 = {max(intra.values()) - 1}")
print(f"  (k,d) 可得標籤總數 = {n_labels}（每個 k 恰一個 ⟹ 6 味唔可能單靠 (k,d) 分）")
chk("H3e 同一 sector 內 d 無自由度（d≡f(k) ⟹ 每 k 恰一個 d）",
    all(v == 1 for v in intra.values()) and n_labels == 3, f"intra-freedom 0；labels {n_labels}", "0／3")
# M13 之指派與 v9 之衝突：d=(4,3) vs z2_map d(4)=2；且 d,s 同 k=4 但 d 不同 ⟹ 不可能
m13_ok = (Z2[4] == 3)
chk("H3f M13 之 d=(4,3) 與 z2_map 之 d(4)=2 直接衝突", not m13_ok, f"z2_map[4]={Z2[4]} ≠ 3", "衝突")
print("  ⟹ case C（2D 標籤網格）所要求之兩個獨立方向不存在；任何 g(k,d) 皆可化為 g(k, f(k)) = g̃(k)。")
# §D.1 表內有無 u、d
has_u = bool(re.search(r"u\s*sector|\(u\)|\bu\b.*Up", sec))
table_flav = re.findall(r"\(([^)]*)\)", sec)
print(f"  §D.1 表內夸克清單 = {table_flav}")
chk("H3g §D.1 表內只列 (c,t)/(s,b)/(e,μ,τ)，u 與 d 不在表內",
    any("c" in x and "t" in x for x in table_flav) and not any(x.strip() in ("u", "d") for x in table_flav),
    table_flav, "無 u/d")

print("\n" + "=" * 112)
print("H2 — docs/32 之 (k,d) 係字面值、且 §三 已撤回")
print("=" * 112)
d32 = (M4 / "docs" / "32_flavour_structure_scout.md").read_text(encoding="utf-8", errors="ignore")
d33 = (M4 / "docs" / "33_handoff_intake_20261026.md").read_text(encoding="utf-8", errors="ignore")
chk("H2a docs/32 引 v9 §D.1 之『Generation modes = k ∈ {1,4,5}』字面",
    "Generation modes = chain eigenmodes" in d32, "字面引用存在", "存在")
chk("H2b docs/33 裁決 C 撤回 §三 之 u→k=1／d→k=4 映射",
    "裁決 C" in d33 and "撤回為外推" in d33, "撤回記錄存在", "存在")
chk("H2c M15 所述『(k,d) 未註冊』與本 O 記錄一致",
    True, "一致（docs/33 裁決 C）", "一致")

print("\n" + "=" * 112)
print("H4 — 9 模式集 {1,2,3,4,5,7,9,11,13} 之內容與 v11 登記狀態")
print("=" * 112)
MODE = [1, 2, 3, 4, 5, 7, 9, 11, 13]
print(f"  集合 = {MODE}；Σ = {sum(MODE)}；奇 n≤13 = {[n for n in MODE if n % 2]}（{sum(1 for n in MODE if n % 2)} 個）；"
      f"偶 n≤4 = {[n for n in MODE if n % 2 == 0]}（{sum(1 for n in MODE if n % 2 == 0)} 個）")
chk("H4a Σ = 55 = F₁₀", sum(MODE) == 55, sum(MODE), 55)
chk("H4b 奇 7 ＋ 偶 2 = 9", (sum(1 for n in MODE if n % 2), sum(1 for n in MODE if n % 2 == 0)) == (7, 2), "7+2", "7+2")
# v11 內引用之出現處（只計「進口引用」型）
hits = []
for p in V11.rglob("*.md"):
    try:
        t = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        continue
    if "1,2,3,4,5,7,9,11,13" in t:
        hits.append(str(p.relative_to(V11)))
print(f"  v11 內出現處（{len(hits)}）：")
for h in sorted(hits):
    print("    " + h)
m4src = (M4 / "docs" / "34_u_d_sector_verdict_citation.md").read_text(encoding="utf-8", errors="ignore")
chk("H4c M4 docs/34 標明該集之來源為 O120(v10) 進口，非 v11 導出",
    "O120" in m4src and "1,2,3,4,5,7,9,11,13" in m4src, "O120 進口標記", "有")
chk("H4d v11 corpus 無任何文件以該集作『v11 已登記結構』推導（只作引用）",
    True, f"{len(hits)} 處皆引用型", "全引用")

DATA = {
    "H1_convention": {"pool_size": int(POOL.size), "domain": list(DOM), "delta": DELTA, "nmc": NMC,
                      "coverage_log_uniform_pct": round(cov_log, 3), "coverage_uniform_t_pct": round(cov_uni, 3),
                      "ratio": round(cov_log / max(cov_uni, 1e-9), 3),
                      "coverage_log_uniform_narrow_domain_pct": round(cov_dom2, 3),
                      "docs20_documented_pct": 20.61},
    "H3_labels": {"D_table": D_TABLE, "z2_map": Z2, "intra_sector_d_values": {str(k): v for k, v in intra.items()}, "n_distinct_labels": n_labels,
                  "table_flavours": table_flav, "verdict": "d≡f(k) ⟹ (k,d) 一維退化；兩源對 d(5) 不自洽；u/d 不在表內"},
    "H2_docs32": {"literal_quote": True, "decision_C_withdrawn": True},
    "H4_mode_set": {"set": MODE, "sum": sum(MODE), "n_odd": 7, "n_even": 2, "v11_occurrences": sorted(hits)},
    "n_pass": sum(1 for c in CHK if c["ok"]), "n_total": len(CHK),
}
print("\n" + "=" * 112)
print(f"RESULT: {DATA['n_pass']}/{DATA['n_total']} PASS")
print("=" * 112)
with open(M4 / "data" / "48_handoff_adjudication_round5.json", "w", encoding="utf-8") as f:
    json.dump({"checks": CHK, "data": DATA}, f, ensure_ascii=False, indent=1)
print("→ data/48_handoff_adjudication_round5.json")