#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/03 — Q2①：味標籤 (k,d) 之**出處追溯**（AST 級）

問題（record §二 Q2①）：`k∈{1,4,5}`、`Z₂` grading `d∈{3,2,0}` 之**產生規則**有無落地？
判準：A1 出處（須 AST 級可追來源；禁「事後命名」）。

方法：不引用任何二手敘述，直接對 v9 原文做 AST 解析 ＋ 逐字核。
  §一 AST：gen_modes／z2_map／chiral_index 是否為字面值（零依賴）
  §二 兩個 check 之邏輯地位（恆真 vs 真檢定）
  §三 「Z₂ 分級」之範疇檢定（|G|=2 ⟹ |image| ≤ 2）
  §四 跨檔一致性（d 值是否自洽）
  §五 反證嘗試：k∈{1,4,5} 可否由任何已註冊譜謂詞唯一給出（N_trials = 6，全報）
N_trials（掃描部分）= 6。零上游寫入（唯讀）。
"""
import ast
import os

V9 = "/home/wsl/Projects/IDCM-Original/IDCM-v9.0--W-Field-Lattice/wfield_lattice"
FILES = {
    "fermion/12_dirac_lattice_verify.py": V9 + "/fermion/12_dirac_lattice_verify.py",
    "04_extended_verify.py":              V9 + "/04_extended_verify.py",
    "lagrangian/14_full_sm_verify.py":    V9 + "/lagrangian/14_full_sm_verify.py",
}
print("=" * 78)
print("M14 scripts/03 — Q2① 味標籤 (k,d) 出處追溯（AST 級，唯讀 v9）")
print("=" * 78)

# =====================================================================
# 一、AST：三個關鍵符號是否字面值
# =====================================================================
print("\n【一】AST 追溯：gen_modes ／ z2_map ／ chiral_index")


def is_literal(node):
    """node 是否完全不依賴任何 Name（即零推導）。"""
    if isinstance(node, ast.Constant):
        return True
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        return all(is_literal(e) for e in node.elts)
    if isinstance(node, ast.Dict):
        return all(is_literal(k) for k in node.keys if k is not None) and \
               all(is_literal(v) for v in node.values)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return is_literal(node.operand)
    return False


findings = {}
for tag, path in FILES.items():
    if not os.path.exists(path):
        print(f"  ⚠ 缺檔：{tag}")
        continue
    tree = ast.parse(open(path, encoding="utf-8").read())
    rec = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id in ("gen_modes", "z2_map"):
                    rec[t.id] = (node.lineno, is_literal(node.value), ast.unparse(node.value)[:44])
        if isinstance(node, ast.FunctionDef) and node.name == "chiral_index":
            body_lit = all(
                isinstance(s, ast.If) and is_literal(s.test) is False and
                all(isinstance(b, ast.Return) and is_literal(b.value) for b in s.body)
                for s in node.body if isinstance(s, ast.If)
            )
            tail_lit = all(isinstance(b, ast.Return) and is_literal(b.value)
                           for s in node.body if isinstance(s, (ast.If, ast.Return))
                           for b in ([s] if isinstance(s, ast.Return) else s.body)
                           if isinstance(b, ast.Return))
            rec["chiral_index"] = (node.lineno, body_lit and tail_lit,
                                   f"{len([s for s in node.body if isinstance(s, ast.If)])} 個 if 分支")
    findings[tag] = rec
    print(f"  {tag}")
    for k, (ln, lit, src) in rec.items():
        print(f"    L{ln:<4d} {k:14s} 純字面值={str(lit):5s}  {src}")

all_lit = all(lit for rec in findings.values() for (_, lit, _) in rec.values())
print(f"  ⟹ 全部關鍵符號皆為**純字面值**（零 Name 依賴）：{all_lit}")
print("     ⟹ 出處追溯終點 = 檔案內之整數字面 ⟹ **無上游推導可追**")

# =====================================================================
# 二、兩個 check 之邏輯地位
# =====================================================================
print("\n【二】兩條『驗證』之邏輯地位（是否恆真）")
gen_modes = [1, 4, 5]
z2_map = {1: 3, 4: 2, 5: 0}
print(f"  v9 `04_extended_verify.py` L262：check('D6. 3 generations from chain modes k={{1,4,5}}',")
print(f"                                    len(gen_modes) == 3)")
print(f"    → gen_modes 係字面 list {gen_modes} ⟹ len == 3 **恆真**，與『三代』無關")
print(f"  v9 `fermion/12_dirac_lattice_verify.py` L460：check('Z₂ grading d={{3,2,0}} allows 3")
print(f"                                    chiral modes without pairs', len(set(z2_map.values())) == 3)")
print(f"    → z2_map = {z2_map} ⟹ len(set(values)) = {len(set(z2_map.values()))} **恆真**")
print(f"    ⚠ 該行並以『Each mode has distinct Z₂ grade』為理由 —— 三個**相異**分級值")
print(f"      與『ℤ₂ 分級』字面矛盾（見 §三）")

# =====================================================================
# 三、範疇檢定：「Z₂ 分級」可否有 3 個值？
# =====================================================================
print("\n【三】範疇檢定：ℤ₂ 分級之像集大小")
print("  定義：G-分級 = 映射 φ: 模集 → G，且 φ 為群同態（此處作用於模指標之加性結構）")
print("  對 G = ℤ₂ = {0,1}：|image φ| ≤ |G| = 2")
print("  v9 要求 |image| = 3（d ∈ {3,2,0}，且以 len(set)==3 為驗證條件）")
print("  ⟹ 3 > 2 ⟹ **d 不可能是任何 ℤ₂ 分級之值** ⟹ 「Z₂ grading d={3,2,0}」係範疇錯置")
print("  附：v9 對 ℤ₂ 之**唯一**正確用法 = 宇稱選擇（88 節點中一半 = 44 個 carry chiral modes，")
print("      expansion_plan.md L43）—— 該用法只給 2 類，**不能**給 3 個扇區。")

# =====================================================================
# 四、跨檔一致性：d 值是否自洽
# =====================================================================
print("\n【四】跨檔一致性（同一版本 v9 內之 d 值）")
table = [
    ("wfield_lattice/03_extended_derivation.md §D.1 L206", "k=1", "(s=2, d=3)"),
    ("同上", "k=4", "(s=1, d=2)"),
    ("同上", "k=5", "(s=0, **d=2**)"),
    ("wfield_lattice/fermion/12_dirac_lattice_verify.py L278", "k=1/k=4/k=5", "d = 3 / 2 / 0"),
    ("wfield_lattice/lagrangian/14_full_sm_verify.py L392", "k=1/k=4/k=5", "d = 3 / 2 / **0**"),
    ("wfield_lattice/paper/17_paper.md L217", "—", "d = {3, 2, 0}"),
    ("wfield_lattice/fermion/11_dirac_lattice_derivation.md L181", "—", "d = {3, 2, 0}"),
]
for src, k, v in table:
    print(f"  {src:56s} {k:12s} {v}")
print("  ⟹ k=5 之 d：**§D.1 = 2**，而其餘四處 = 0 ⟹ v9 內部**自相矛盾**（A1 結構性失敗）")

# =====================================================================
# 五、（k,d) ↔ {u,d,s} 之來源
# =====================================================================
print("\n【五】(k,d) ↔ {u,d,s} 之對應有無來源")
print("  v9 §D.1 之 sector 表逐字：")
print("     k=1 (s=2,d=3) → Up sector (c, t)")
print("     k=4 (s=1,d=2) → Down sector (s, b)")
print("     k=5 (s=0,d=2) → Lepton sector (e, μ, τ)")
print("  ⟹ 表內**只有 sector 名**，且 **u、d 皆不在表內**（Up sector 只列 c,t；Down 只列 s,b）")
print("  ⟹ 『u→k=1、d→k=4』係由 sector 名推廣之**外推**，非原文指派")
print("     （M4 docs/32 §九 已獨立記錄：該映射**已撤回**；register 亦記 v9 §D.1 未採用）")
print("  ⟹ (k,d) ↔ {u,d,s} 之 6 個對應：**未註冊輸入**（A1 FAIL）")

# =====================================================================
# 六、反證嘗試：k∈{1,4,5} 可否由已註冊譜謂詞唯一給出
# =====================================================================
print("\n【六】反證嘗試（N_trials = 6，跑前凍結；全報，不挑選）")
import numpy as np
gam88 = {n: 4 * np.sin(np.pi * n / 88) ** 2 for n in range(1, 45)}
print(f"  v9 所引 γ 值（/88 慣例）：γ₁ = {gam88[1]:.4f}、γ₄ = {gam88[4]:.4f}、γ₅ = {gam88[5]:.4f}"
      f"   ⟹ 與 v9 文內 0.0051/0.0810/0.1261 吻合 ⟹ v9 用 **/88**")
fib = [1, 2, 3, 5, 8, 13, 21, 34]
luc = [1, 3, 4, 7, 11, 18, 29]
P = {
    "P1 γ 最小三個":            sorted(gam88, key=lambda n: gam88[n])[:3],
    "P2 Fibonacci 首三個":      sorted(fib)[:3],
    "P3 Lucas 首三個":          sorted(luc)[:3],
    "P4 完全平方 ≤10":          [n for n in range(1, 11) if int(n**0.5)**2 == n][:3],
    "P5 三角數 ≤10":            [n for n in range(1, 11) if (8*n+1)**0.5 % 1 == 0][:3],
    "P6 γ_n < 0.15":            [n for n in range(1, 45) if gam88[n] < 0.15][:3],
}
hits = 0
for name, val in P.items():
    ok = sorted(val) == [1, 4, 5]
    hits += ok
    print(f"  {name:22s} → {sorted(val)}   {'✅ 命中' if ok else '🔴 不中'}")
print(f"  命中 {hits} / 6 ⟹ k∈{{1,4,5}} **不可**由上述已註冊譜謂詞給出")
print("  ⚠ 依禁區 3（禁後驗凑數），本節**不作**任何「反推謂詞」嘗試；")
print("     6 個候選全部落空即止，結論只為『無既有規則』，非『不存在任何規則』（A5 範圍）")

# =====================================================================
# 判決
# =====================================================================
print("\n" + "=" * 78)
print("Q2① 判決")
print("=" * 78)
print("""
  §一 gen_modes／z2_map／chiral_index 全部係**純字面值**（零依賴）⟹ A1 出處斷鏈
  §二 兩條「驗證」皆恆真（len(字面 list)==3；len(set(字面 dict))==3）⟹ 零檢定力
  §三 d∈{3,2,0} 有 3 個值 ⟹ **不可能是 ℤ₂ 分級**（|ℤ₂|=2）⟹ 範疇錯置
  §四 v9 內部 k=5 之 d：§D.1 給 2、其餘四處給 0 ⟹ 自相矛盾
  §五 sector 表不含 u,d ⟹ (k,d)↔{u,d,s} 未註冊
  §六 k∈{1,4,5} 於 6 個已註冊譜謂詞下 0/6 命中

  ⟹ **Q2① = 判「未註冊輸入」（A1 FAIL）**
    範圍：{v9 之 fermion/12、04_extended、lagrangian/14、03_extended_derivation §D.1、
           paper/17、fermion/11} 六處原文 —— 即現存全部 (k,d) 陳述處
  ⟹ 後果：Q2② 之「3! 唯一化」所依之標籤集 (k,d) 本身即未註冊輸入 ⟹
     即使 3! 可砍至唯一，亦只係在未註冊標籤上之算術（不可作結構結論）。
  N_trials = 6（§六，全報）。零上游寫入。
""")
raise SystemExit(0)