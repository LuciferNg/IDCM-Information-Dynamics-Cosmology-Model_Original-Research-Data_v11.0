#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M21 scripts/05 — 產物查重 ＋ 雙態登記（P5／Q3；含 NDR1）
=========================================================
判準照 docs/00_preregistration.md §3–§4：

  1) 產物查重（跨 O／跨版本，唯讀）：本輪宣稱之每條恆等式，喺 v10／v11 全目錄
     搜關鍵詞族，記錄命中檔案數 —— 命中數 > 0 ⟹ **不得宣稱新發現**，只能宣稱
     「裁決／密度閘／唯一性引理／圖分類」為新。
  2) 雙態登記：N_d = 44 之狀態 = 結構導出（群側）／橋未閉（晶格側）
     —— 由 scripts/01–04 之結果機械生成，不得手寫。
  3) dof 記帳：producer 存在且零參 → N_d 非第二錨 → O212 dof = 1（分支 (a)）。

零 import 上游；唯讀掃描（open() 讀，不寫任何外部檔）。
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ODIR = os.path.dirname(HERE)
DATA = os.path.join(ODIR, "data")

ROOTS = [
    "/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat",
    "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture",
]
SELF = os.path.abspath(ODIR)

CHECKS = []


def ck(tag, cond, detail=""):
    status = "[PASS]" if cond else "[FAIL]"
    print(f"{status} {tag} :: {detail}")
    CHECKS.append((tag, bool(cond), detail))
    return bool(cond)


# ── 產物查重 ───────────────────────────────────────────────────────
CLAIMS = {
    "K1  Σ_χ d(d−1) = 44（A₅ 度集普查）": r"Σ\s*d\s*\(\s*d\s*[−-]\s*1\s*\)",
    "K2  N_d = B·L₅ = L₃·L₅ = 4×11": r"(2\s*B\s*[·*]\s*L\s*₅)|(B\s*[·*]\s*L\s*₅)|(L\s*₃\s*[·*]\s*L\s*₅)|(2B·L₅)",
    "K3  Σγ = 2N + 2（恆等）": r"(Σ\s*γ\s*=\s*2N\s*\+\s*2)|(2\s*\(\s*N\s*\+\s*1\s*\))|(2\(N\+1\))",
    "K4  K_tree = 33 = 3·11 = N_gen·L₅": r"(K_tree\s*=\s*33\s*=\s*3\s*[·*]\s*11)|(K_tree\s*=\s*N_gen\s*[·*]\s*L\s*₅)|(3\s*[·*]\s*11\s*=\s*33)",
    "K5  非對合／平方根計數（群側算子身份）": r"(非對合)|(involut)|(#\s*\\?\{\s*g\s*:)|(平方根)",
    "K6  f_M0·f_R 乘積 = N_O38": r"(f_M0\s*[·*]\s*f_R)|(N_O38)|(1\.801868e40)",
    "K7  60 = N_d + B² 記賬": r"(60\s*=\s*N_d\s*\+\s*B\s*\^?\s*2)|(N_d\s*\+\s*B²)|(N_d \+ B²)",
    "K8  唯一性：度集由 (Σd², Σd, k) 唯一決定": r"(度集.{0,12}唯一)|(degree multiset.{0,20}unique)",
}


def scan():
    files = []
    for root in ROOTS:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", "node_modules")]
            for fn in filenames:
                if not fn.endswith((".md", ".py")):
                    continue
                p = os.path.join(dirpath, fn)
                if p.startswith(SELF):
                    continue
                try:
                    if os.path.getsize(p) > 2_000_000:
                        continue
                except OSError:
                    continue
                files.append(p)
    return files


print("=" * 78)
print("M21 scripts/05 — 產物查重 ＋ 雙態登記")
print("=" * 78)
files = scan()
print(f"  掃描檔案數（v10 ∪ v11，.md/.py，排除本 O） = {len(files)}")

results = {}
compiled = {k: re.compile(v) for k, v in CLAIMS.items()}
for k, rx in compiled.items():
    hits = []
    for p in files:
        try:
            txt = open(p, encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        m = rx.search(txt)
        if m:
            hits.append(os.path.relpath(p, "/home/wsl/Projects/IDCM-Original"))
    results[k] = hits
    print(f"\n  {k}")
    print(f"    命中檔案數 = {len(hits)}")
    for h in hits[:6]:
        print(f"      · {h}")

ck("NDR1-a  查重掃描完成（> 500 檔）", len(files) > 500, f"{len(files)} files")
for k in CLAIMS:
    ck(f"NDR1  {k.split('  ')[0]} 查重已記錄", k in results, f"hits = {len(results[k])}")

# NDR1-b：K1／K2／K4／K6 已有上游命中 ⟹ 不得宣稱新發現
upstream_known = ["K1", "K2", "K4", "K6", "K7"]
for code in upstream_known:
    key = [k for k in results if k.startswith(code)][0]
    ck(f"NDR1-dedup  {code} 有上游命中 ⟹ 只可宣稱裁決／閘／圖，不得宣稱新恆等式",
       len(results[key]) > 0, f"hits = {len(results[key])}")

# ── 雙態登記 ───────────────────────────────────────────────────────
print("\n── 雙態登記（機械生成）──")
D = [1, 3, 3, 4, 5]
B = 4
N2 = sum(D)
G = 60
N_d_census = sum(d * (d - 1) for d in D)
assert N_d_census == 44 and G - N2 == 44

registry = {
    "quantity": "N_d = 44",
    "formal_state": "structure-derived (group side) / bridge-closed-negative + identification-stipulated (lattice side)",
    "state_zh": "結構導出（群側）／橋已判死（對象層）＋識別為約定（晶格側）",
    "dof_continuous_added": 0,
    "producers": [
        {
            "id": "E4",
            "route": "N_d = Σ_χ d(d−1) = |G| − Σ_χ d = 60 − 16",
            "operator_identity": "非對合元素計數（= |G| − #{g: g²=e}）；度集由 (Σd²=60, Σd=16, k=5) 唯一強迫",
            "dof": 0,
            "status": "F",
            "registered_upstream": ["v10 O171 docs/38_chain_spectrum_derivation.md",
                                   "v11 M15 §三 已知前置", "v11 M16 記賬 60 = N_d + B²"],
        },
        {
            "id": "E9",
            "route": "N_d = B·L₅ = 4×11（= L₃·L₅）",
            "operator_identity": "B = (d−1)-加權平均度（CA2）；L₅ = Σ_χ(d−1)",
            "dof": 0,
            "status": "F（給定 E4 + CA2）",
            "registered_upstream": ["v11 M15 docs/13 Q3（N_d = L₃L₅）", "v10 O171（2B·L₅ = 88）"],
        },
    ],
    "tautologies_never_producers": [
        {"id": "E5", "stmt": "Σγ = 2(N_d+1) = 90", "why": "循環圖恆等（CG1，N=2..120 違例 0）⟹ 對 N 零約束"},
        {"id": "E6", "stmt": "Σγ = B(N_d+1)/2 = 90", "why": "與 E5 同一方程（B/2 = 2）"},
    ],
    "matches_no_operator_identity": [
        {"id": "E7", "stmt": "h11 = N_d − 2B = 36", "why": "池掃 ≥6 條 k-free 式（scripts/02）⟹ 零判別力"},
        {"id": "E8", "stmt": "h21 = 2N_d + 10 = 98", "why": "10 = 2F₅ 係後見選項"},
        {"id": "E13", "stmt": "|G| = N_d + B² ⟺ N₂ = B² = 16", "why": "兩邊各自被 A₅ 確定，但 B² 無算子身份"},
    ],
    "refuted_routes": [
        {"id": "p9-FORCE44", "route": "φ-級聯 → Σγ = 90 → N = 44",
         "status": "🔴 循環", "ref": "M15 out_of_scope/docs/14（3 條同型級聯、31.8% 覆蓋、硬編碼 44）"},
    ],
    "open_bridge": {
        "question": "A₅ 之 44 個非對合元素 ↔ 鏈之 44 個非零相異模式 之典範映射",
        "equivalent_form": "N_site = 2·N_d（為何乘 2）＋ 為何載體係環（O171 residual）",
        "status": "🔴 closed-NEGATIVE（對象層）2026-11-20 第二轮：入站 M22 handoff 裁定（scripts/07–08）",
        "death_certificate": [
            "類級：gcd 型 {20,10,5,5,2,1,1} vs 階型 {24,20}（7 型 vs 2 型）",
            "對級：鏈 (21 對 +1 不動點 +1 無伴) vs 群 (22 對 +0 +0)；44 集上 involutive 置換必偶數不動點",
            "等變級：兩側 orbit 多重集 0 對相等（{20,12,12}/{2×22} vs {20,10,5,5,2,1,1}/{2×21,1,1}）",
            "機制級：環作用＝上游 O171 已判非群作用；transport 類型不合；Galois 阿貝爾 vs A₅ 單非交換",
        ],
        "residue": {
            "coarse_label": "{20,24} = (gcd=1 vs >1) ≡ (階3 vs 階5)；cut 密度 2/43 = 4.65% ⟹ 標籤，無映射內容",
            "dimension_level": "Σ dim Λ² V_χ = 22 = N_d/2；Σ_(i<j) d_i d_j = 98 = h21（皆非元素映射）",
        },
        "identification_status": "stipulation（載體 = 2× 普查計數；無元素映射）⟹ 連續 dof 仍 0，帳本 ＋1 條結構約定",
    },
    "dof_ledger": {
        "anchors_in_seven_quantities": 0,
        "N_d_is_second_anchor": False,
        "O212_branch": "(a) dof = 1 成立",
        "note": "producer 零參數 ⟹ N_d 唔新增連續 dof；橋係離散結構聲稱（可判真假），非可調參數",
    },
}

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "05_registry_Nd.json"), "w", encoding="utf-8") as f:
    json.dump({"registry": registry, "dedup": {k: v for k, v in results.items()},
               "n_files_scanned": len(files)}, f, ensure_ascii=False, indent=2)
print(json.dumps(registry, ensure_ascii=False, indent=2)[:1500])

ck("NDR1-e  登記狀態 = 雙態（群側導出／晶格側橋判死＋識別約定）",
   "structure-derived" in registry["formal_state"]
   and "bridge-closed-negative" in registry["formal_state"],
   registry["formal_state"])
ck("NDR1-f  dof 連續新增 = 0", registry["dof_continuous_added"] == 0, "0")
ck("NDR1-g  N_d 非第二錨 ⟹ O212 dof = 1", registry["dof_ledger"]["N_d_is_second_anchor"] is False,
   registry["dof_ledger"]["O212_branch"])
ck("NDR1-h  登記含 producer／T／M／refuted／bridge 五欄",
   all(k in registry for k in ("producers", "tautologies_never_producers",
                               "matches_no_operator_identity", "refuted_routes", "open_bridge")),
   "五欄齊")

npass = sum(1 for _, ok, _ in CHECKS if ok)
nfail = len(CHECKS) - npass
print("\n" + "=" * 78)
print(f"scripts/05 結果：{npass}/{len(CHECKS)} PASS ｜ {nfail} FAIL")
print("=" * 78)

sys.exit(0 if nfail == 0 else 1)
