#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M15 scripts/09 — 第十候選之降級審計：後驗 ＋ 反演恆真 ＋ 有理身份閘
================================================================================
觸發：docs/09 判第十候選（`N_gen = B − 1`）為 🟡。本輪依兩條已存在之硬約束重審：
  (a) M13 `docs/13 §三` **自認**：迹 law `3/8 = (B−1)/(2B)` 係「見到 0.047σ 之後寫出」
      ⟹ 後驗（post-hoc），且「provenance 仍未閉 ⟹ 不得升 🟢」。
  (b) prior session 之**有理身份閘**（本 session 不可重推已降級之路線）：
      `dev/基準 = |迹−3/8| ÷ 1/(2q²)`，q = 16 ⟹ 1.339 ≥ 1 ⟹ 「相容，未識別」。

§1  來源層：兩處引文（檔案級 grep；v11 只讀）。
§2  ★ 反演恆真審計：r → B → N_gen = B−1 → r_pred，證明 r_pred ≡ r（零資訊）。
§3  有理身份閘：對 k = 2/3/4 之預言值，計 dev/基準（q = 分母）。
§4  判詞：第十候選 🔴（後驗 ＋ 恆真）；剩餘內容之精確表述。
§5  M15 verdict 回復 negative。

零新觀測；`N_trials = 0`。
"""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

M13 = Path("/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/"
           "main_problems/m13_flavour_orthogonality")
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
say("M15 scripts/09 — 第十候選降級審計（後驗／恆真／有理身份閘）")
say("=" * 78)

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§1  來源層：M13 docs/13 之自認（檔案級引文）")
say("-" * 78)
QUOTES = [
    ("迹 law 之結構形式係後驗讀出", "後驗：結構形式見數據後寫出"),
    ("不可宣稱為事前預言", "不可稱事前預言"),
    ("provenance 仍未閉", "provenance 未閉 ⟹ 不得升 🟢"),
]
txt = (M13 / "docs/13_symmetric_trace_closure.md").read_text()
hits = []
for q, desc in QUOTES:
    ok = q in txt
    hits.append(ok)
    say(f"  [{'✓' if ok else '✗'}] 「{q}」—— {desc}")
chk("§1 M13 自認後驗（三處實證）", all(hits),
    "⟹ `N_gen = B − 1` 及 `3/8 = (B−1)/(2B)` 屬**後驗**，不可作事前預言（docs/09 之 §二 立論失效）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§2  ★ 反演恆真審計：反演係恆等，零資訊")
say("-" * 78)
say("  程序（docs/09 §3b）：由實測 r 解 (B−1)/B² = r 得 B，再報 N_gen = B − 1。")
say("  檢：把該 N_gen 代回結構式 N_gen/B²，是否得回 r？—— 依構造必然。")
r_obs = Fr("0.1888076409")
sol = []
for cand in (Fr("3.95837"), Fr("3.958376"), Fr("3.9583765")):
    lhs = (cand - 1) / cand ** 2
    sol.append(abs(float(lhs) - float(r_obs)))
say(f"  三個數值 B 之回代殘差：{['%.2e' % x for x in sol]} ⟹ 收斂至 r（非獨立檢定）")


def solve_B(r):
    d = 1 - 4 * r
    return (1 + d ** 0.5) / (2 * r)


B_back = solve_B(float(r_obs))
Ng = B_back - 1
r_pred = (Ng) / B_back ** 2
say(f"  閉式反演：r = {float(r_obs):.10f} ⟹ B = {B_back:.8f} ⟹ N_gen = {Ng:.8f}")
say(f"  回代 N_gen/B² = {r_pred:.10f} ⟹ 差 = {abs(r_pred - float(r_obs)):.2e}")
chk("§2 反演恆真（回代差 < 1e−9）", abs(r_pred - float(r_obs)) < 1e-9,
    "⟹ 反演只係重寫；「N_gen = 2.958 ± 0.185」不含任何獨立信息")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§3  有理身份閘（prior session 判準；本 session 不得重推）")
say("-" * 78)
say("  判準：dev = |v_pred − r_obs|/v_pred；基準 = 1/(2q²)（q = v_pred 之最簡分母）；")
say("        dev/基準 ≥ 1 ⟹ 該有理值「無身份」（與隨機有理不可分）。")
r_obs_f = float(r_obs)
s_r = 0.005841
say("")
say(f"  {'k':>2} | {'預言值':>12} | {'q':>2} | {'dev%':>7} | {'基準%':>7} | {'dev/基準':>8} | 判")
gate = {}
for k in (2, 3, 4):
    v = Fr(k, (k + 1) ** 2)
    q = v.denominator
    dev = abs(float(v) - r_obs_f) / float(v)
    bench = 1.0 / (2 * q * q)
    ratio = dev / bench
    gate[k] = (float(v), q, dev, bench, ratio)
    mark = "✅ 過閘" if ratio < 1 else "🔴 無身份"
    say(f"  {k:>2} | {float(v):>12.6f} | {q:>2} | {dev * 100:>6.3f} | {bench * 100:>6.3f} | "
        f"{ratio:>8.3f} | {mark}")
say("")
say(f"  ★ k = 3（即 3/16，q = {gate[3][1]}）：dev/基準 = {gate[3][4]:.3f} ≥ 1 ⟹ **無有理身份**")
say("    ⚠ 與 prior session 記錄之 1.339 **數值不一致**（其自報 dev 0.6974% ÷ 基準 0.1953% = 3.570）；")
say("      兩者之**判詞相同**（≥ 1 ⟹ 無身份），但比值不同 ⟹ 其 1.339 疑屬算術滑落或基準定義不同。")
say("      本輪以可重算之定義為準（dev 對預言值、基準 1/(2q²)），標 CONV-DIFF，不改 prior 記錄。")
chk("§3 有理身份閘：3/16 未過閘", gate[3][4] >= 1.0,
    f"dev/基準 = {gate[3][4]:.3f} ⟹ 0.224σ 只可稱「相容」，不可稱「閉合／識別」"
    f"（prior 報 1.339，判詞同、比值異 ⟹ CONV-DIFF）")
say(f"  （附：k=2 之 dev/基準 = {gate[2][4]:.2f}、k=4 之 = {gate[4][4]:.2f} ⟹ 皆無身份；")
say("    故此閘**不**判別 k，只否定「3/16 已被識別」之宣稱。）")

# ─────────────────────────────────────────────────────────────────────────────
say("")
say("§4  判詞")
say("-" * 78)
say("  (1) docs/09 §二 之立論：「零參數預言 3/16、z = +0.224σ」⟹ 依 §1 為**後驗**；")
say("      依 §3 為「相容，未識別」。兩者皆禁止升格。")
say("  (2) docs/09 §3b 之「B = 3.958 ± 0.185 ⟹ 整數唯一 3」⟹ 依 §2 為**恆等反演**，零資訊。")
say("  (3) ⟹ 第十候選 **🔴 後驗 ＋ 恆真**（非擬合失敗，而係立論方向錯：由數據反寫恆等式再宣稱預言）。")
say("  (4) 剩餘之**合法**內容（可保留）：")
say("      「若採 M13 之恒等式 N_gen = B − 1，則 Δ_N/Δ_Ξ 之實測容許之整數世代數唯一為 3")
say("       （k=2 −5.2σ、k=4 +5.6σ）」—— 此為**條件性命題**，非 N_gen = 3 之推導。")
say("  (5) ⟹ M15 verdict 由 `partial` 回復 `negative`（十候選全滅；四載體 scoped no-go 不變）。")

with open(__file__.replace("scripts/", "data/").replace(".py", ".json"), "w") as f:
    json.dump({"results": RES, "r_obs": str(r_obs), "B_back": B_back, "Ngen_back": Ng,
               "r_pred": r_pred, "rational_gate": {str(k): gate[k] for k in gate},
               "M13_self_admission": bool(all(hits))}, f, ensure_ascii=False, indent=1)

nfail = sum(1 for r in RES if r["status"] == "FAIL")
say("")
say(f"§5  總結：檢查 {len(RES)} 項，FAIL = {nfail}")
sys.exit(1 if nfail else 0)