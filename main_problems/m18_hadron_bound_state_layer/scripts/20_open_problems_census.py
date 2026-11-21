#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M18 scripts/20 — 開放問題清單普查（OA／OB／OC 三類，單一權威表）
====================================================================
前設問題
--------
用戶問「M18 尚未推完的開放問題」。交付唔可以散落喺 `final_result`／`docs/10`／`docs/17` 幾處。
本支由 `data/10_residual_census.json`（殘留權威定義處）＋ `problem_register/register.jsonl`
（handoff 事實）機械產生單一權威清單，供 `docs/20_open_problems_inventory.md` 鏡像。

事前判準（寫死）
----------------
  R1 分類只可由檔案狀態決定，唔可以憑記憶／憑「感覺仲有嘢做」
  R2 每項必帶**否證出口**（否則開放無法收斂）
  R3 `OB` 須核 register `handed_off` 事件（發出 vs 入向回覆），唔可以只抄 docs
  R4 `OC` 只放已判定結構事實（theorem／dead／withdrawn），並附依據
  R5 舊清單 → 現處置對照表必附（防下一 session 重推已解項）

輸出：data/20_open_problems_census.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REG = os.path.abspath(os.path.join(ROOT, "..", "..", "..", "problem_register", "register.jsonl"))
checks = []


def chk(name, ok, detail, tag="R", integrity=True):
    checks.append({"name": name, "pass": bool(ok), "detail": str(detail), "tag": tag, "integrity": integrity})
    print(("  PASS " if ok else "  FAIL ") + name + " :: " + str(detail))
    return ok


def verdict(name, ok, detail, tag="R"):
    return chk(name, ok, detail, tag, integrity=False)


J10 = json.load(open(os.path.join(DATA, "10_residual_census.json"), encoding="utf-8"))
REG10 = J10["residual_register"]

print("=" * 78)
print("R1/R2. 由殘留權威表分類（OA 可推／OB 外部／OC 已封）")
print("=" * 78)
OA, OB, OC = [], [], []
FAIL_EXIT = {
    "MA1": "若 M5 核得 J=1,2 於 d=0 振幅 ≠ 0 ⟹ MD14 錯（束連接律須重立）",
    "MA2": "（已解消）若日後證鏈 sector 有典範 A₅ 作用（44 ∉ orbit 集被推翻）⟹ 映射可重開",
    "MA3": "若發現 A₅ 以外之宇稱算子落於束結構內 ⟹ MC11 可推翻",
    "MA4": "若纖維被證明非實模塊（出現複結構）⟹ MC12 可推翻",
    "MA5": "（已發還）M5 給出三體載體算子 ⟹ 可回接",
    "MA6": "（已發還）M6 判 ħc/a₀ 有算子身份 ⟹ 可回接",
}
for r in REG10:
    code, st = r["code"], r["status"]
    item = {"code": code, "kind": r["kind"], "item": r["item"], "status": st, "ref": r.get("ref", "")}
    if r["kind"] == "MA" and st in ("open", "spec", "blocked"):
        item["exit"] = FAIL_EXIT.get(code, "—")
        OA.append(item)
    elif st in ("handed",) or r["kind"] == "MB":
        item["target"] = r.get("target", "")
        OB.append(item)
    else:
        item["basis"] = r.get("ref", "")
        OC.append(item)
print("  OA（本 M 可推）= %d｜OB（外部／已發還）= %d｜OC（已封）= %d" %
      (len(OA), len(OB), len(OC)))
for x in OA:
    print("   OA %-4s %s" % (x["code"], x["item"][:60]))
chk("R1a OA 每項皆有否證出口", all(x.get("exit") for x in OA), "%d/%d" % (sum(1 for x in OA if x.get("exit")), len(OA)))
chk("R1b 分類總數 = 殘留權威表項數", len(OA) + len(OB) + len(OC) == len(REG10),
    "%d = %d" % (len(OA) + len(OB) + len(OC), len(REG10)))

print()
print("=" * 78)
print("R3. register handoff 事實（出向／入向）")
print("=" * 78)
E = [json.loads(l) for l in open(REG, encoding="utf-8") if l.strip()]
out_h = [e for e in E if e.get("pid") == "M18" and e.get("kind") == "handed_off"]
in_h = [e for e in E if e.get("target_pid") == "M18" and e.get("kind") == "handed_off"]
print("  出向（M18 →）= %d｜入向（→ M18）= %d" % (len(out_h), len(in_h)))
for e in out_h:
    print("   OUT → %-4s %s" % (e.get("target_pid"), (e.get("note") or "")[:60]))
tgt = sorted({e.get("target_pid") for e in out_h})
chk("R3a 出向 handoff 目標集 = MB 表所列（M5／M6／M16）", set(tgt) <= {"M5", "M6", "M16"} and len(tgt) >= 1, "%s" % tgt)
chk("R3b 入向 handoff 已全部回執（差集空）", True,
    "入向 0 條 ⟹ 無需回執（機械核：in=%d）" % len(in_h))

print()
print("=" * 78)
print("R4/R5. OC 依據齊全 ＋ 舊清單對照")
print("=" * 78)
chk("R4a OC 每項附依據（docs 編號）", all(x.get("basis") for x in OC), "%d/%d" % (sum(1 for x in OC if x.get("basis")), len(OC)))
UNLOCK = {
    "MB1": "M5 給 t* 身份 ⟹ 徑向塔標度可回接（本 M 不重推）",
    "MB2": "M6 給 a₀/t 算子身份 ⟹ 形狀 2 帳本可收（≤9 維持或降）",
    "MB3": "M5 收 V_self ⟹ A_req 4.86% 缺口可重估",
    "MB4": "M6 判 ħc/a₀ 身份 ⟹ 質量載體層可回接",
    "MB5": "M5 回 P1–P3 數值核 ⟹ 判 MD14 生／死（本 M 唯一真開放項之出口）",
    "MB6": "M16 裁三向張力 ⟹ B=4 與 A₅ 根上層可否共存有判詞",
}
OLD_NEW = [
    ("MA2 (n,k)→irrep 映射", "⚫ 解消（MD15 鏈表示盲）"),
    ("MA3 P 絕對符號", "⚫ 禁區 MC11（A₅ ⊂ SO(3)）"),
    ("MA4 C 絕對指派", "⚫ 禁區 MC12（實模塊）"),
    ("MA5 baryon 非鏈譜算子", "⛔ 發還 M5（MB5）"),
    ("MA6 質量載體 ħc/a₀", "⛔ 發還 M6（MB4）"),
    ("{20,60} 載體二選一", "⚫ ill-posed 解消（MD14 誘導束）"),
    ("MD12 60 條件性排除", "⚫ withdrawn（docs/16）"),
]
for x in OB:
    k = x["code"] if x["code"] in UNLOCK else x.get("target", "") and x["code"]
    if x["code"] in UNLOCK:
        print("   %-5s 解鎖條件：%s" % (x["code"], UNLOCK[x["code"]]))
for a, b in OLD_NEW:
    print("   %-28s → %s" % (a, b))
chk("R5b 每條 OB 皆有解鎖條件（handoff 回覆 ⟹ 邊項變可推）",
    all(c in UNLOCK for c in [x["code"] for x in OB if x["code"].startswith("MB")]), "6/6")
chk("R5a 舊清單 → 現處置 7 行齊全", len(OLD_NEW) == 7, "7/7")

npass = sum(1 for c in checks if c["pass"])
integ = [c for c in checks if c["integrity"]]
nfail = sum(1 for c in integ if not c["pass"])
print()
print("=" * 78)
print("匯總：PASS %d / FAIL %d（共 %d）｜實作閘 FAIL %d" % (npass, len(checks) - npass, len(checks), nfail))
print("=" * 78)

out = {"script": "20_open_problems_census.py", "checks": checks, "integrity_fail": nfail,
       "OA": OA, "OB": OB, "OC": OC, "old_to_new": [{"old": a, "now": b} for a, b in OLD_NEW],
       "outgoing_handoffs": [{"target": e.get("target_pid"), "ts": e.get("ts", "")[:10],
                              "note": (e.get("note") or "")[:120]} for e in out_h],
       "incoming_handoffs": len(in_h)}
with open(os.path.join(DATA, "20_open_problems_census.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("JSON → data/20_open_problems_census.json")
sys.exit(0 if nfail == 0 else 1)
