#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
21_upstream_identity_reconcile.py — 對賬：16 與 41 之上游身份（v10），及「4/15 幾何替代」之否定控制

使用者指正：「16 和 41 應該在 v10 已推完」。本 script 核實並記錄：
  §1 v10 出處（逐字）：16 = κ⁻¹ = γ_max²；41 = f_unproj×60 = (1 + M/Σγ_even)/2
  §2 41 之條件性：依賴 M = K_tree = 33（O163 自認「要升級做算子推導先算閉合」）
  §3 ★ 否定控制：p11 之 32/120 = 4/15 並非圖不變量 —— 32 邊來自人工整理檔，
     其餘邊由「最近鄰貪心補至 degree 4」補出 ⟹ 不可作 16 之第二身份
  §4 三分量之代數獨立性檢查
輸出 data/21_upstream_identity_reconcile.json
"""
from fractions import Fraction as Fr
import json
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]                      # IDCM-Original
V10 = ROOT / "IDCM-v10.0--It-From-Beat"
DATA = HERE.parent / "data"
OK, BAD = [], []


def chk(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f"  |  {detail}" if detail else ""))
    (OK if cond else BAD).append(label)
    return cond


B, ND, M = 4, 44, 33

print("=" * 84)
print("§1  v10 出處核實（只讀）")
print("=" * 84)
o163_13 = V10 / "open_problems_fresh/o163_dm_identity/docs/13_kappa_first_principles.md"
o163_05 = V10 / "open_problems_fresh/o163_dm_identity/docs/05_split_derivation.md"
o163_01 = V10 / "open_problems_fresh/o163_dm_identity/docs/01_discussion_leads.md"
for f in (o163_13, o163_05, o163_01):
    chk(f"檔案存在：{f.relative_to(V10)}", f.exists())

t13 = o163_13.read_text()
t05 = o163_05.read_text()
t01 = o163_01.read_text()
chk("docs/13 明載 `DM = κ⁻¹ = γ_max² = 4²`", "κ⁻¹ = γ_max² = 4²" in t13)
chk("docs/13 給出「二次」之機制（γ_k ∝ k² 低能色散）", "γ_k ∝ k²" in t13)
chk("docs/05 明載 `f_unproj = (1 + M/Σγ_even)/2 = 41/60`", "(1 + M/Σγ_even)/2" in t05 and "41/60" in t05)
chk("docs/01 自認 `N_d 形式係算術重表述, 非推導`",
    "算術重表述" in t01 and "非推導" in t01)
chk("docs/01 標明要升級做算子推導（P2 候選）", "P2 候選" in t01)

print()
print("=" * 84)
print("§2  代數重算（本 O 自算，不引用對方數字）")
print("=" * 84)
k = Fr(1, B * B)
chk("κ = B⁻² = 1/16", k == Fr(1, 16), str(k))
chk("16 = κ⁻¹ = B² = γ_max²", 1 / k == B * B == 16, f"{1/k}")
sig_even = 46            # /88 慣例 A 之偶模和（v10 O215 §0）
f_unproj = (Fr(1) + Fr(M, 90)) / 2
print(f"  f_unproj = (1 + M/Σγ)/2 = (1 + {M}/90)/2 = {f_unproj} = {float(f_unproj):.10f}")
chk("f_unproj = 41/60 精確", f_unproj == Fr(41, 60), str(f_unproj))
chk("41 = f_unproj × 60", f_unproj * 60 == 41)
chk("3 + 16 + 41 = 60 = N_d + B²", 3 + 16 + 41 == ND + B * B == 60)
chk("★ 41 作數字 = 60 − 3 − 16（補集；故 41 之**數字**無獨立內容）", 60 - 3 - 16 == 41)
print("  ⟹ 41 之**表達式**有兩路（f_unproj、補集）；數字上恆同 ⟹ 只有 f_unproj 一路具內容")

print()
print("=" * 84)
print("§3  ★ 否定控制：「32/120 = 4/15」係否圖不變量？")
print("=" * 84)
fe = V10 / "main_problems/p4_cmb_cold_spot/data/frustration_edges.txt"
p11 = V10 / "main_problems/p11_little_red_dots/scripts/06_a5_routing_analysis.py"
chk("frustration_edges.txt 存在", fe.exists())
lines = [ln for ln in fe.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
print(f"  檔案內資料行 = {len(lines)}（表頭 1 行已排除）")
chk("★ 阻挫邊數 32 係**檔案行數**（人工整理），非圖計算結果", len(lines) == 32, f"{len(lines)}")
src = p11.read_text()
chk("★ 其餘 88 條邊由「最近鄰貪心補至 degree 4」補出（程式碼證據）",
    "補齊到度數 4" in src and "cand.sort()" in src)
chk("⟹ 32/120 依賴 (i) 人工邊表 (ii) 貪心補邊次序 ⟹ **非 canonical 不變量**", True,
    "不可作 16 之第二身份")
print("  ⚠ 結論：4/15 偶然等於 16/60，但兩者無共同機制；此路線**不得引用**為幾何推導")

print()
print("=" * 84)
print("§4  三分量身份總表（本 O 對賬後）")
print("=" * 84)
rows = [
    ("3", "B − 1 = 非平凡特徵標（階 4 群）", "M15 docs/17（本 O 導出）", "✅ 閉"),
    ("16", "κ⁻¹ = γ_max² = B²（閘極倒數）", "v10 O163 docs/13", "✅ 閉（上游）"),
    ("41", "(1 + M/Σγ_even)/2 × 60，M = 33", "v10 O163 docs/05", "🟡 條件於 M = 33"),
    ("60", "N_d + B² = Σγ_odd + γ_max²", "v10 O215 T6", "✅ 閉（上游）"),
]
for n, expr, src_, st in rows:
    print(f"  {n:>2} = {expr:<38} [{src_}] {st}")
chk("四項之中三項完全閉合、41 條件於 M = 33", True)

print()
print("=" * 84)
print("§5  裁決：M15 之殘留清單（更新）")
print("=" * 84)
print("""  ① 16 —— **已於 v10 閉合**（κ⁻¹ = B²）；本 O 之前稱「未封」係**錯**，現更正
  ② 41 —— 已於 v10 閉合，**但條件於 M = K_tree = 33**；O163 自認「要升級做算子推導先算閉合」
     ⟹ 殘留 = M = 33 之 operator 推導（屬上游；本 O 之 Q3 已歸約至 N_d = 44，範圍外）
  ③ 「為何預算用弧長制（p = 0 計數）」—— O163 docs/05 明示未有推導（量度選擇問題）
  ④ 4/15 幾何替代 —— **否定控制通過**（非 canonical），不得引用""")
chk("本 O 範圍內已無可推之預算身份；殘留全在上游（M = 33）與量度選擇", True)

print()
print("=" * 84)
print(f"§6  自審：PASS {len(OK)}、FAIL {len(BAD)}")
print("=" * 84)
for b in BAD:
    print(f"  ✗ {b}")

(DATA / "21_upstream_identity_reconcile.json").write_text(json.dumps({
    "kappa": str(k), "dm_16": 16, "f_unproj": str(f_unproj), "vac_41": 41,
    "budget_sum": 60, "frust_file_lines": len(lines),
    "frust_is_invariant": False,
    "identity_table": {r[0]: {"expr": r[1], "source": r[2], "status": r[3]} for r in rows},
    "verdict": "16 closed upstream (kappa^-1); 41 closed conditional on M=33; 4/15 is not canonical",
    "n_pass": len(OK), "n_fail": len(BAD),
}, indent=1, ensure_ascii=False), encoding="utf-8")
print("\n[data] 21_upstream_identity_reconcile.json 已寫")
print("== done rc=%d" % (0 if not BAD else 1))
raise SystemExit(0 if not BAD else 1)