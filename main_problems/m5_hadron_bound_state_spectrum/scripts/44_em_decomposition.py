"""M5 docs/42 配套 — 「通用電荷分裂 A_EM」之真偽檢定（Dashen 分解）
================================================================================
觸發：`docs/39` 提出 `|ΔM|(π/K/D) = 4.4492 ± 0.4603 MeV` 為「通用 EM 項」並列為新首務。
本輪自審：符號係 **−, +, −, +**（交替）⟹ 可能係「EM + 夸克質量」兩項抵消之巧合，非通用常數。
判法＝**Dashen 分解**（`A` 由外部固定，不擬合）：`ΔM_f = s_f·A_EM + q_f`，`A_EM = 4.5936 MeV`（＝`|ΔM_π|`，
Dashen 定理之首階內容），`s_f = −1`（電荷態較重，π 定符號）。反解 `q_f`，檢其係否隨 `1/m_h` 走。

判準（事前寫死）：
  D1 `A_EM` 取 `4.5936 MeV`（Dashen，**不擬合**）；敏感度另掃 `A ∈ [4.3, 5.0]`。
  D2 `ΔM_f` 一律 `M(neutral) − M(charged)`，PDG 2024 直讀；`s_f = −1` 全族（EM 令電荷態較重）。
  D3 夸克項 `q_f = ΔM_f − s_f·A_EM`。
  D4 重夸克標度檢：`q_B/q_D` 應 ≈ `m_c/m_b = 1.27/4.18 = 0.304`（`1/m_h`）；`q_D/q_K` 應 ≈ `m_s/m_c = 0.495/1.27 = 0.390`。
  D5 判：若任一 `|ln(實測比/預期比)| > ln 3` ⟹ 該族與 `1/m_h` 標度不相容 ⟹ 標為 anomaly。
  D6 不得以擬合 `A_EM` 令 `q_f` 變靚；本輪只判生死。
輸出：`data/44_em_decomposition.json`
"""
import json
import math
from pathlib import Path

A_EM = 4.5936
M = {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    for k in ids:
        if k:
            M[int(k)] = mv * 1000.0

FAM = [("π", 111, 211, 0.005), ("K", 311, 321, 0.495), ("D", 421, 411, 1.27), ("B", 511, 521, 4.18)]
print("=" * 100)
print("M5 docs/42 — 通用 A_EM 之 Dashen 分解檢定")
print("=" * 100)
print(f"\n[D1] A_EM = {A_EM} MeV（Dashen，不擬合）")
print(f"\n[D2/D3] 分解表（ΔM = neutral − charged；s = −1）")
rows = []
for nm, i0, i1, mh in FAM:
    dM = M[i0] - M[i1]
    q = dM - (-1) * A_EM
    rows.append({"fam": nm, "dM": dM, "q": q, "m_h": mh})
    print(f"  {nm:2s}  ΔM = {dM:+8.4f}  ⟹ q_f = {q:+8.4f} MeV（m_h ≈ {mh} GeV）")

print("\n[D4] 夸克項之 1/m_h 標度檢")
qd = {r["fam"]: r["q"] for r in rows}
checks = [("q_D/q_K", qd["D"] / qd["K"], 0.495 / 1.27),
          ("q_B/q_D", qd["B"] / qd["D"], 1.27 / 4.18)]
verdicts = {}
for name, actual, expect in checks:
    lnr = math.log(abs(actual) / abs(expect)) if actual and expect else float("inf")
    ok = abs(lnr) <= math.log(3)
    verdicts[name] = {"actual": actual, "expect": expect, "ln_ratio": lnr, "consistent": bool(ok)}
    print(f"  {name}: 實測 {actual:+.4f} vs 1/m_h 預期 {expect:+.4f}"
          f"（ln 比 {lnr:+.3f}）⟹ {'✅ 相容' if ok else '🔴 不相容'}")

print("\n[D5] 敏感度（A_EM 掃 4.3–5.0，看 q_B 能否靠 A 救回）")
for A in (4.30, 4.5936, 4.80, 5.00):
    qb = (M[511] - M[521]) + A
    print(f"  A = {A:.4f} ⟹ q_B = {qb:+7.4f}；q_B/q_D = {qb / ((M[421] - M[411]) + A):+7.3f}"
          f"（1/m_h 預期 0.304）")

print("\n[D5b] |ΔM| 聚類之偶然性")
vals = [abs(r["dM"]) for r in rows[:3]]
mean = sum(vals) / 3
sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / 2)
print(f"  |ΔM|(π,K,D) = {vals}、均值 {mean:.4f}、樣本 sd {sd:.4f}（{100*sd/mean:.1f}%）")
print(f"  ⟹ 以 3 個數、sd/mean = {100*sd/mean:.1f}% 論「通用」，統計上不成立（3 點無法區分常數與抵消）")

anomaly = [k for k, v in verdicts.items() if not v["consistent"]]
out = {"A_EM": A_EM, "rows": rows, "scaling_checks": verdicts, "anomalies": anomaly,
       "clustering": {"vals": vals, "mean": mean, "sd": sd}}
with open("data/44_em_decomposition.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print(f"\n[總結] 與 1/m_h 不相容者：{anomaly if anomaly else '無'}")
print("-> data/44_em_decomposition.json")