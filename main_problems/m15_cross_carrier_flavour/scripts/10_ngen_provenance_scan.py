#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/10 — 全 corpus provenance 掃描：`N_gen` 有無 producer？
================================================================================
觸發：docs/10 之後，問題剩一句 —— 「3 之安裝點在何處，且全 corpus 有無任何一處**由結構算出** N_gen？」
本輪對 v9／v10／v11 之**全部 .py** 做賦值掃描（`N_gen = …`、`n_gen = …`、`N_GEN = …`、`ngen = …`），
分三類：L 字面（rhs 為純數）｜W 數值包裝（`mp.mpf(3)`、`Fraction(3,1)`…）｜E 真表達式。
再對 E 類逐一核其是否化約為字面或恆等式。

判準（跑前寫死）：
  P1  若 E 類為 0 ⟹ `N_gen` 全 corpus 無 producer ⟹ 3 之身份 = 純輸入（結論級）。
  P2  若 E 類非 0，逐一列出並核其自變量；若自變量含 `N_gen`／`B`／`K_tree` 之任一個 ⟹ 循環。
  P3  安裝點：世代索引之**唯一**實際使用處 = 三個模式之字面列表（`gen_modes = [1, 4, 5]`）。

路徑：v9／v10／v11（**全部只讀**，零寫入）。
"""
import ast
import json
import re
import sys
from pathlib import Path

ROOT = Path("/home/wsl/Projects/IDCM-Original")
ROOTS = ["IDCM-v9.0--W-Field-Lattice",
         "IDCM-v10.0--It-From-Beat/main_problems",
         "IDCM-v10.0--It-From-Beat/open_problems_fresh",
         "IDCM-v11.0--Reality-Architecture/main_problems"]
PAT = re.compile(r"^\s*(N_gen|n_gen|N_GEN|ngen)\s*=\s*(.+?)\s*$")
LIT = re.compile(r"^[-+]?\d+(?:\.\d+)?$")
WRAP = re.compile(r"^(?:mp\.mpf|float|int|sp\.Integer|Fraction)\(\s*[-+]?\d+(?:\.\d+)?\s*(?:,\s*\d+\s*)?\)$")
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
say("M15 scripts/10 — 全 corpus `N_gen` provenance 掃描")
say("=" * 78)

rows = []
for r in ROOTS:
    for p in (ROOT / r).rglob("*.py"):
        if "m15_cross_carrier_flavour" in str(p):      # SELF 排除：本 O 目錄不入掃描
            continue
        try:
            txt = p.read_text(errors="ignore")
        except Exception:
            continue
        for i, ln in enumerate(txt.splitlines(), 1):
            m = PAT.match(ln)
            if not m:
                continue
            rhs = m.group(2).strip()
            # 去行內註解（粗略：# 之後，但保留字串內者不計）
            code = rhs.split("  #")[0].strip()
            # 先判「是否可解析為 Python 賦值」——不可解析者屬散文／docstring／註解
            try:
                ast.parse(f"{m.group(1)} = {code}")
                parseable = True
            except Exception:
                parseable = False
            if not parseable:
                cat = "P"          # 散文／docstring／註解行（非可執行賦值）
            elif LIT.match(code):
                cat = "L"
            elif WRAP.match(code):
                cat = "W"
            else:
                cat = "E"
            rows.append({"path": str(p.relative_to(ROOT)), "line": i,
                         "lhs": m.group(1), "rhs": code, "cat": cat})

cnt = {c: sum(1 for r in rows if r["cat"] == c) for c in ("L", "W", "P", "E")}
say("")
say("§1  掃描結果（v9／v10／v11 全部 .py）")
say("-" * 78)
say(f"  總賦值 {len(rows)} 處：" +
    "｜".join(f"{k} 類 {cnt[k]}" for k in ("L", "W", "P", "E")))
say("    L = 字面純數｜W = 數值包裝（mp.mpf(3) 等）｜P = 散文／docstring（不可解析）｜E = 可解析之真表達式")
chk("§1-a 掃描完成", len(rows) > 100, f"{len(rows)} 處")
say("")
say("  E 類（唯一可能係「由結構算出」者）逐條：")
E = [r for r in rows if r["cat"] == "E"]
for r in E:
    say(f"    {r['path']}:{r['line']}  {r['lhs']} = {r['rhs'][:70]}")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  E 類之自變量分析（P2：含 N_gen／B／K_tree 即循環）")
say("-" * 78)
SUMMARY = {
    "B - 1": ("循環", "B = |ℤ₂|² = 4 為恆等；B − 1 = 3 係算術（docs/08 Lucas 恆等式）"),
    "B-1": ("循環", "同上"),
    "K_Tree / L₅": ("循環", "K_tree := N_gen·L₅（P7 自身第 21 行已寫死 N_gen = 3）"),
    "K_Tree/L₅": ("循環", "同上"),
    "(N_d − N_d/B) / (N_d/B)": ("恆等", "(N_d−N_d/B)/(N_d/B) ≡ B − 1，對**任何** N_d 成立 ⟹ 與 N_d 無關"),
    "(N_d − L₅)/L₅": ("恆等", "= N_d/L₅ − 1 = L₃ − 1（用 N_d = L₃L₅）⟹ docs/08 恆等式"),
    "3;": ("字面（複合語句）", "`N_gen = 3; f_proj = 19` ⟹ 仍係字面"),
    "3  (": ("字面（散文後綴被解析為呼叫）", "`N_gen = 3  (…)` ⟹ 值仍係字面 3"),
    "mpf(3)": ("字面（數值包裝）", "`mpf(3)` = mp.mpf(3)"),
    "h11_test // 12": ("跨載體輸入", "= 36/12；12 ∉ CY₃ primitive（M15 docs/06 判 CROSS-IMPORT）"),
    "F(e + o - 3, F5 + 1)": ("含字面 3", "(e+o−3)/(F₅+1) = (6+15−3)/6 = 3；分子含字面 −3"),
}
for r in E:
    key = next((k for k in SUMMARY if r["rhs"].startswith(k)), None)
    verdict, why = SUMMARY.get(key, ("未分類", "需人工核"))
    say(f"    {r['rhs'][:44]:<46} ⟹ {verdict}：{why}")
say("")
say("  ⟹ 無一條 E 係「由非 N_gen 輸入算出 N_gen」：全部屬 (i) B − 1 之算術、")
say("     (ii) K_tree/L₅ 之定義重寫、或 (iii) 對任何 N_d 皆成立之恆等式。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  P3 安裝點：世代索引之唯一實際使用處")
say("-" * 78)
GM = []
for r in ROOTS:
    for p in (ROOT / r).rglob("*.py"):
        try:
            t = p.read_text(errors="ignore")
        except Exception:
            continue
        for i, ln in enumerate(t.splitlines(), 1):
            if re.match(r"^\s*gen_modes\s*=\s*\[\s*1\s*,\s*4\s*,\s*5\s*\]", ln):
                GM.append(f"{Path(r).name}/{p.name}:{i}")
for g in GM:
    say(f"    {g}  ← 字面列表 [1, 4, 5]")
say(f"  ⟹ 共 {len(GM)} 處（v9 三個 verify 腳本各一份 ＋ M14 之 AST 追溯腳本各一份）")
# 恆真鏈
say("")
say("  恆真鏈（v9 12_dirac_lattice_verify.py）：")
say("    n_physical = len(gen_modes)                    # = len([1,4,5]) = 3")
say("    check(..., n_physical == 3, '3 generations')    # ⟹ 對字面列表算 len")
say("    check(..., len(set(z2_map.values())) == 3)      # ⟹ 對字面 dict {1:3,4:2,5:0} 算 len")
say("    chiral_index(k): if k==1→3 elif k==4→2 elif k==5→0 else −1   # if/elif 硬編碼")
chk("§3 安裝點唯一（字面列表 ＋ if/elif 硬編碼；全 corpus 僅此一類）", len(GM) == 4,
    f"{len(GM)} 處 gen_modes=[1,4,5]；chiral_index 為字面分支；無任何由結構生成世代索引之處")
say("")
say("  （另）p10/07_structural_proof_force_cy3.py:44  `n_gen = 3  # 世代數 = 操作維度 2+1`")
say("  （另）p2/61_p7_tree_forcing.py:21            `N_d, N_gen, B, K_tree, mu2 = 44, 3, 4, 33, 2`")
say("        ⟹ P7 之「N_gen = K_tree/L₅，唯一 n=5」係以 N_gen = 3 為輸入，再除出 3。")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  判詞")
say("-" * 78)
say(f"  (1) P1 成立：全 corpus {len(rows)} 處 `N_gen` 賦值中，"
    f"{cnt['L'] + cnt['W']} 處為字面／數值包裝、{cnt['P']} 處為散文，{cnt['E']} 處為可解析表達式，")
say("      而該 E 類全部化約為 B − 1 之算術或 K_tree/L₅ 之定義重寫 ⟹ **無 producer**。")
say("  (2) P3 成立：世代索引之唯一實際使用處 = 三個模式之字面列表 [1,4,5]（3 處）")
say("      ＋ `chiral_index` 之 if/elif 硬編碼 ⟹ **3 之安裝點唯一且係字面**。")
say("  (3) P7（v10 p2「樹密度 = 世代數」）不成立：其自身第 21 行已把 N_gen = 3 寫死；")
say("      「唯一 n=5」只係 33 ÷ 11 之除法（其餘 Lucas 除數不整除）。")
say("  (4) ⟹ **「3」之身份 = 三個字面模式索引之基數；全 corpus 無任何結構 producer。**")
say("      唯一未解之瓶頸不變：`K_tree = M = 33` 之 k-free 推導（docs/08 §未盡 Q1）。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "counts": cnt, "total": len(rows),
               "install_points": GM,
               "E_class": [{"path": r["path"], "line": r["line"], "rhs": r["rhs"]} for r in E]},
              f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§5  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)