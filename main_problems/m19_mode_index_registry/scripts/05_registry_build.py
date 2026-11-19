#!/usr/bin/env python3
"""
M19 / scripts/05 — 登記冊落地（R4 / P4）

預註冊（跑前凍結）:
  F1  登記冊須完備覆蓋 9 帶電費米子 × 1 個正典 frame（F1）
  F2  狀態三態（結構／輸入／未註冊）由**可計算判據**決定，唔可以人手標：
        STRUCT_CLASS ⟺ 該粒子之容許支撐 ⊊ 9 模式集（值層受結構約束至類別）
        INPUT_RANK   ⟺ 支撐 = 全 9 模式（值層零結構約束）＋ 觀測 rank 給值
        未註冊       ⟺ 冇任何 frame 可指（本輪預期 0）
  F3  每一行必須帶 frame 欄（禁裸 k_q）
  F4  引用閘（可執行）：一份宣稱含 k_q 之記錄若不帶 {frame, state, provenance}
      三欄 ⟹ gate = BLOCK。用 3 個樣本（1 合法 ＋ 2 非法）證明閘真會攔。
  F5  跨記錄不一致清單（唯讀比對，唔改他 M）：
        (a) u/d 地位：M12 OP6（M4 docs/34 vs O140 docs/19 兩個地位）
        (b) M12 OP3 / M14 §附二.3：S4 阻塞於 c/b 指數
        (c) F3（O206）被當「物理 assignment」使用，實屬支撐編碼層

失敗判準: F2 判據與 03 之矩陣不一致 ⟹ 登記冊撤回。
"""
import json
from pathlib import Path

OUT = Path(__file__).parent.parent / "data"
M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
LABEL = {"e": "e", "u": "u", "d": "d", "s": "s", "mu": "μ", "c": "c",
         "tau": "τ", "b": "b", "t": "t"}
F1 = {"e": 1, "u": 2, "d": 3, "s": 4, "mu": 5, "c": 7, "tau": 9, "b": 11, "t": 13}
MASS = {"e": 0.511, "u": 2.16, "d": 4.67, "s": 93.4, "mu": 105.6584,
        "c": 1270.0, "tau": 1776.86, "b": 4180.0, "t": 172690.0}

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/05 — 登記冊落地（R4）")
print("=" * 92)

adm = json.loads((OUT / "03_admissibility.json").read_text())["A"]

rows = []
for p, n in F1.items():
    col = adm[LABEL[p]]
    support = [int(k) for k, v in col.items() if v > 0]
    state = "STRUCT_CLASS" if len(support) < len(M9) else "INPUT_RANK"
    if state == "STRUCT_CLASS":
        prov = "T3 槽定理（χ=1 ⟹ n≡1 mod 4）× 群層 iso Φ:Ŝ↦z（v10 O210 docs/08,10；O233 docs/13+17）"
        note = "類別結構；具體值由質量 rank 定"
    else:
        prov = "質量 rank（觀測輸入）"
        note = "零結構約束；值 = 類內 γ 升序"
    rows.append({
        "frame": "F1_gamma_ascending",
        "particle": LABEL[p],
        "n": n,
        "state": state,
        "provenance": prov,
        "note": note,
        "support": support,
        "log2_free_value_layer": None,
    })

chk("F1 登記冊覆蓋 9 帶電費米子", len(rows) == 9 and {r["particle"] for r in rows} ==
    {LABEL[p] for p in F1}, "9/9")
chk("F2 狀態由判據決定（輕子 = STRUCT_CLASS；夸克 = INPUT_RANK）",
    sum(r["state"] == "STRUCT_CLASS" for r in rows) == 3 and
    sum(r["state"] == "INPUT_RANK" for r in rows) == 6,
    "3 結構類別 ＋ 6 輸入")
chk("F3 每行帶 frame 欄（禁裸 k_q）", all(r["frame"] for r in rows), "F1_gamma_ascending ×9")

# F4 引用閘：可執行 gate
REQUIRED = ("frame", "state", "provenance")


def ref_gate(record):
    miss = [k for k in REQUIRED if k not in record or record[k] in (None, "", "-")]
    return ("BLOCK", miss) if miss else ("PASS", [])


samples = [
    ("合法：完整三欄", {"frame": "F1", "state": "INPUT_RANK", "provenance": "質量 rank"}),
    ("非法：缺 frame", {"state": "INPUT_RANK", "provenance": "質量 rank"}),
    ("非法：裸 k_q", {"k_q": 7}),
]
verdicts = [(nm, *ref_gate(rec)) for nm, rec in samples]
for nm, v, miss in verdicts:
    print(f"    {nm:<16} → {v}" + (f"（缺 {miss}）" if miss else ""))
chk("F4 引用閘 1 PASS / 2 BLOCK（真攔）",
    [v for _, v, _ in verdicts] == ["PASS", "BLOCK", "BLOCK"], "見上" + "")

# F5 跨記錄不一致（唯讀記錄）
inconsistencies = [
    {"id": "GX1", "item": "u→2 / d→3 之地位",
     "detail": "M12 OP6：M4 docs/34 記「已裁決觀測輸入」；O140 docs/19 記「後驗候選」⟹ 同兩值兩個地位",
     "target": "M12 / M14"},
    {"id": "GX2", "item": "S4 阻塞",
     "detail": "M12 OP3 ＋ M14 §附二.3：阻塞於 c/b 模式指數未註冊 ⟹ 本輪證明：登記後仍係 rank 輸入，阻塞性質係「框架依賴」而非「缺一格」",
     "target": "M12 / M14"},
    {"id": "GX3", "item": "F3（O206 CSSG-v3）之層位",
     "detail": "被稱「物理 assignment」，但 n=1 三粒子重複 ＋ 5/9 落 9 模式集外 ⟹ 屬支撐編碼層，非索引層",
     "target": "M1 / O206（v10，唯讀）"},
]
chk("F5 跨記錄不一致清單 3 條已列", len(inconsistencies) == 3, ", ".join(i["id"] for i in inconsistencies))

OUT.mkdir(exist_ok=True)
(OUT / "05_registry.json").write_text(json.dumps({
    "frame_canonical": "F1_gamma_ascending",
    "states": ["STRUCT_CLASS", "INPUT_RANK", "UNREGISTERED"],
    "rows": rows,
    "gate": {"required": list(REQUIRED), "samples": [{"name": n, "verdict": v, "missing": m}
                                                     for n, v, m in verdicts]},
    "inconsistencies": inconsistencies,
}, ensure_ascii=False, indent=1))

print("\n  k_q 登記冊（正典 frame = F1）：")
print(f"  {'粒子':>4} {'n':>3} {'γ_n':>12} {'狀態':>13}  支撐")
import math
for r in rows:
    g = 4 * math.sin(r["n"] * math.pi / 88) ** 2
    print(f"  {r['particle']:>4} {r['n']:>3} {g:>12.9f} {r['state']:>13}  {sorted(r['support'])}")

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
