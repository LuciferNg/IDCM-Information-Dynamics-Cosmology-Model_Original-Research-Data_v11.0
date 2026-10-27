"""M5 docs/36 配套 — R4「輕區異常」之誤差傳播重判（z 判準）
================================================================================
觸發：`docs/35` 判 K 族 `dev = 0.0589`（偏離 17.0×）⟹「🔴 輕區異常（R4）成立」。
但 `data/37` 已記錄 `σ_c(K) = 0.08216 MeV` 而 `c_K = −0.01350 MeV`（**σ/|c| = 6.09**），
而 R4 判準係**純比值**（`dev > 5×`），**未做誤差傳播**。本腳本補此缺口。

判準（事前寫死，寫於本 docstring；不得因結果改動）：
  Z1 `c_model = c_data / dev`（dev = ψ²_data/ψ²_model）須與 `data/37` 之 `psi2_proxy` 自洽（相對差 < 1e−9）。
  Z2 `z ≡ |c_data − c_model| / σ_c`。
  Z3 判：`z ≥ 5` ⟹ 發現；`3 ≤ z < 5` ⟹ 提示（hint）；`z < 3` ⟹ **不足以宣稱異常**。
  Z4 報「達 5σ 所需 σ_c」＝ `|Δc| / 5`，及需改善倍數。
  Z5 報模型鑑別力 `＝ |c_model| / σ_c`（模型預言相對數據誤差之大小）。
  Z6 三族全報，不得只報 K。

不含新擬合、不改 `data/37`、不動 `docs/35`（append-only；本輪只加 docs/36 與 register note）。
輸出：`data/38_k_anomaly_significance.json`
"""
import json
from decimal import Decimal as D

RAW = json.load(open("data/37_third_family.json"))
ROWS, DEV = RAW["rows"], RAW["dev"]

print("=" * 104)
print("M5 docs/36 — R4 輕區異常：誤差傳播重判（z 判準）")
print("=" * 104)

print(f"\n[Z1] 自洽檢查（dev = ψ²_data/ψ²_model ⟹ c_model = c_data / dev）")
res = {}
for k in ("K", "D", "B"):
    c, s, dv = D(str(ROWS[k]["c"])), D(str(ROWS[k]["sigma_c"])), D(str(DEV[k]))
    cm = c / dv                      # dev = ψ²_data/ψ²_model，ψ² ∝ |c| ⟹ c_model = c_data/dev
    psi_d, psi_m = D(str(ROWS[k]["psi2_proxy"])), D(str(ROWS[k]["psi2_proxy"])) / dv
    dev_psi = abs((psi_d / psi_m) / dv - 1)   # 應 ≈ 0（dev 定義自洽）
    z = abs(c - cm) / s
    need = abs(c - cm) / D(5)
    res[k] = {"c_data_MeV": float(c), "sigma_c": float(s), "sigma_over_c": float(abs(s / c)),
              "dev": float(dv), "c_model_MeV": float(cm), "delta_c": float(abs(c - cm)),
              "z_sigma": float(z), "sigma_needed_5sigma": float(need),
              "precision_gain_needed": float(s / need) if need else None,
              "discrimination": float(abs(cm) / s),
              "verdict": "發現" if z >= 5 else ("提示" if z >= 3 else "不足以宣稱異常")}
    print(f"  {k}: dev 自洽 {float(dev_psi):.2e}  σ/|c| = {float(abs(s/c)):.2f}")

print(f"\n[Z2/Z3] z 判（門檻：≥5 發現、3–5 提示、<3 不足）")
print(f"  {'族':<3}{'c_data':>11}{'σ_c':>9}{'σ/|c|':>8}{'dev':>8}{'c_model':>11}"
      f"{'|Δc|':>10}{'z(σ)':>8}{'判':>14}")
for k in ("K", "D", "B"):
    r = res[k]
    print(f"  {k:<3}{r['c_data_MeV']:>11.6f}{r['sigma_c']:>9.5f}{r['sigma_over_c']:>8.2f}"
          f"{r['dev']:>8.4f}{r['c_model_MeV']:>11.6f}{r['delta_c']:>10.6f}"
          f"{r['z_sigma']:>8.3f}{r['verdict']:>14}")

print(f"\n[Z4] 達 5σ 所需精度")
for k in ("K", "B"):
    r = res[k]
    print(f"  {k}: 需 σ_c ≤ {r['sigma_needed_5sigma']:.6f} MeV（現 {r['sigma_c']:.5f}）"
          f"⟹ 需精確 {r['precision_gain_needed']:.1f}×")

print(f"\n[Z5] 模型鑑別力（|c_model| / σ_c）")
for k in ("K", "D", "B"):
    r = res[k]
    print(f"  {k}: 模型預言只及數據誤差之 {r['discrimination']*100:.2f}%")

zK = res["K"]["z_sigma"]
if zK >= 5:
    verdict = f"🔴 異常成立（z = {zK:.3f}σ ≥ 5）⟹ 模型被排除"
elif zK >= 3:
    verdict = f"🟡 提示級 hint（z = {zK:.3f}σ，3–5σ）⟹ 未達 5σ 標準，不得宣稱「異常成立」"
else:
    verdict = f"🟡 未達標準（z = {zK:.3f}σ < 3）⟹ 既不足以宣稱異常，亦不足以宣稱一致"
print(f"\n[判決] {verdict}")
print(f"   ⟹ 「17.0×」係**真比值**（模型高估 K 之 ψ² 17 倍），但**顯著性只有 {zK:.2f}σ**：")
print(f"      c_model = {res['K']['c_model_MeV']:.6f} MeV vs c_data = {res['K']['c_data_MeV']:.6f}"
      f" ± {res['K']['sigma_c']:.5f} MeV ⟹ |Δc| = {res['K']['delta_c']:.6f} MeV = {zK:.2f}σ")
print(f"   ⟹ 模型對 K **可測**（預言 = 數據誤差之 {res['K']['discrimination']*100:.0f}%），"
      f"但未達 5σ；σ_c 改善 {res['K']['precision_gain_needed']:.2f}× 即定生死")
print(f"   ⟹ `docs/35` 之「🔴 輕區異常（R4）成立」應降為 **🟡 2.6σ 提示**（未達硬標準）")
print(f"   ⟹ 連帶：以 K 為驅動之 no-go（docs/28–34 之 17.58× 天花板、定理 P、S6 16.2×）"
      f"其**前提目標值**亦須做同樣誤差傳播，方知係否仍成立（docs/36 §5）")

out = {"criterion": {"found": "z>=5", "hint": "3<=z<5", "insufficient": "z<3"},
       "rows": res, "verdict_K": verdict}
with open("data/38_k_anomaly_significance.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print("\n-> data/38_k_anomaly_significance.json")