"""M5 docs/37 配套 — K 側結論之 σ 傳播總審（執行 docs/36 §5）
================================================================================
觸發：`docs/36`（外部 session）指出 K 族 `σ_c/|c| = 6.09` ⟹ `docs/35` 之「17× 異常」只得 2.63σ，
並列 §5「docs/28–34 之 no-go 全部以 K 為驅動，須重算顯著性」為未做項。本腳本執行該項。

判準（事前寫死）：
  S1 σ 一律由 PDG 原始誤差按象限合成；`A = (Δ0+3Δ1)/4 − X`，`X = 2.56 ± 0.094`（`docs/09`）。
  S2 相對誤差按一階傳播；比值 `r = A/c` 之 σ 用 `sqrt((σA/A)² + (σc/c)²)`。
  S3 分級（照 `docs/36` Z3）：z ≥ 5 發現｜3 ≤ z < 5 提示｜z < 3 不足以宣稱。
  S4 B 之 `Δ1` 來自 PDG 單一合併條目（`513 523`）⟹ 三個 σ 變體全報（見 §3），不得只揀一個。
  S5 不改任何既有數值，只加顯著性維度。
輸出：`data/39_sigma_propagation.json`
"""
import json
import math
from pathlib import Path

MASS, ERR = {}, {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    try:
        ev = abs(float(line[52:60]))
    except ValueError:
        ev = 0.0
    for k in ids:
        if k:
            MASS[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0

X, XS = 2.56, 0.094
FAM = {"K": ((311, 321), (313, 323)), "D": ((411, 421), (413, 423)), "B": ((511, 521), (513, 523))}

print("=" * 100)
print("M5 docs/37 — K 側結論 σ 傳播總審")
print("=" * 100)

rows = {}
print("\n[S1] 逐族輸入（PDG 2024 原始誤差）")
print(f"  {'族':3s}{'Δ0':>10s}{'σΔ0':>9s}{'Δ1':>10s}{'σΔ1':>9s}{'c':>10s}{'σ_c':>9s}{'σ/|c|':>8s}")
for k, (g0, g1) in FAM.items():
    d0, d1 = MASS[g0[0]] - MASS[g0[1]], MASS[g1[0]] - MASS[g1[1]]
    s0 = math.hypot(ERR[g0[0]], ERR[g0[1]])
    s1 = math.hypot(ERR[g1[0]], ERR[g1[1]])
    c = (d1 - d0) / 4
    sc = math.hypot(s0, s1) / 4
    A = (d0 + 3 * d1) / 4
    sA = math.hypot(math.hypot(s0, 3 * s1) / 4, XS)
    rows[k] = {"d0": d0, "s0": s0, "d1": d1, "s1": s1, "c": c, "sc": sc,
               "A_raw": A, "A": A - X, "sA": sA,
               "r": (A - X) / c if c else float("nan"),
               "sr_over_r": math.hypot(sA / (A - X), sc / c) if c else float("nan")}
    r = rows[k]
    print(f"  {k:3s}{d0:10.4f}{s0:9.4f}{d1:10.4f}{s1:9.4f}{c:10.5f}{sc:9.5f}{abs(sc/c):8.2f}")
for k, r in rows.items():
    r["sr"] = abs(r["r"]) * r["sr_over_r"]
    print(f"  {k}: A = {r['A']:.5f} ± {r['sA']:.5f}、r = A/c = {r['r']:.3f} ± {r['sr']:.1f}"
          f"（相對誤差 {r['sr_over_r']*100:.1f}%）")

print("\n[S2] K 側各項主張之顯著性（模型值 vs 數據 ± σ）")
CLAIMS = [
    ("c_K 模型 −0.229293 vs 數據", -0.229293, rows["K"]["c"], rows["K"]["sc"], "docs/35 R4／docs/36"),
    ("r_K 模型 1.4580 vs 數據", 1.4580, rows["K"]["r"], rows["K"]["sr"], "docs/23 D5／docs/25"),
    ("A_K 上界 1.24948 vs 數據", 1.24948, rows["K"]["A"], rows["K"]["sA"], "docs/34 定理 P"),
    ("4|c| 下限 0.8726 vs 數據 4|c|", 0.8726, 4 * abs(rows["K"]["c"]), 4 * rows["K"]["sc"], "docs/29 S6"),
]
print(f"  {'主張':38s}{'模型/門檻':>12s}{'數據':>12s}{'σ':>10s}{'z':>8s}{'分級':>12s}")
grades = {}
for name, model, data, sigma, src in CLAIMS:
    z = abs(model - data) / sigma if sigma else float("inf")
    g = "發現" if z >= 5 else ("提示" if z >= 3 else "不足")
    grades[name] = {"model": model, "data": data, "sigma": sigma, "z": z, "grade": g, "src": src}
    print(f"  {name:38s}{model:12.4f}{data:12.4f}{sigma:10.4f}{z:8.2f}{g:>12s}")

print("\n[S3] 天花板 R = r_K/r_D（比值之 σ 合成）")
srK, srD = rows["K"]["sr"], rows["D"]["sr"]
R = rows["K"]["r"] / rows["D"]["r"]
sR = abs(R) * math.hypot(srK / rows["K"]["r"], srD / rows["D"]["r"])
R_model = 1.6470
zR = abs(R - R_model) / sR
print(f"  數據 R = {R:.3f} ± {sR:.1f}（相對 {abs(sR/R)*100:.1f}%）；模型天花板 R_max = {R_model}")
print(f"  ⟹ z = {zR:.2f} ⟹ 分級 {'發現' if zR>=5 else ('提示' if zR>=3 else '不足')}"
      f"；「短缺 {abs(R)/R_model:.2f}×」係中央值比值，非統計結論")

print("\n[S4] B 之 Δ1 三個 σ 變體（PDG 513/523 單一合併條目）")
b0 = FAM["B"][0]
s0b = math.hypot(ERR[b0[0]], ERR[b0[1]])
variants = {"quadrature（docs/36 現用）": math.hypot(s0b, rows["B"]["s1"]) / 4,
            "同條目 ⟹ σ(Δ1)=0": s0b / 4,
            "未解析帶 ±0.2、均勻 σ=0.2/√3": math.hypot(s0b, 0.2 / math.sqrt(3)) / 4}
cB_model = -0.0796036682616
out_b = {}
for nm, sc in variants.items():
    z = abs(cB_model - rows["B"]["c"]) / sc
    out_b[nm] = {"sigma_c": sc, "z": z}
    print(f"  {nm:34s} σ_c = {sc:.5f} ⟹ z = {z:.3f} ⟹ {'一致' if z <= 2 else '張力'}")

out = {"rows": rows, "claims": grades, "R": {"value": R, "sigma": sR, "model": R_model, "z": zR},
       "B_sigma_variants": out_b, "X": X, "XS": XS}
with open("data/39_sigma_propagation.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/39_sigma_propagation.json")