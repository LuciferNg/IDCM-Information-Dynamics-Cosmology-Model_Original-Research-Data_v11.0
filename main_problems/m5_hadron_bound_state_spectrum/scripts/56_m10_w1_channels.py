"""M5 docs/54 配套 — M10 W1 接口：`em(n−p)` 之兩條通道能否補 1.38% 缺口？
================================================================================
M10（`docs/14` W1）要求：`em(n) − em(p) ≥ 1.2667 MeV`（FLAG 中性）或 `1.3963`（`γ` 比值），
而 M5 譜側兩體上限 `A_max(v2) = 1.2494 MeV` ⟹ 差 **1.38%**，M10 Tier II 判詞因此永久阻塞。

本輪問：**有無第二條通道可以補？** 唯一候選 = 自能（`∝ q²`，格點上以 `g(0)` 加權、有限）。
判準（事前寫死）：
  W1 兩體上限 `A_max = α(ħc/a₀)·|Δ(q_h q_l)|·⟨B⟩`，`⟨B⟩ ≤ 1` ⟹ 精確值 `α(ħc/a₀)/3`。
  W2 自能通道之號：`Δ(Σq²) ≡ Σq²(p) − Σq²(n)`（`p = uud`、`n = udd`），
     自能恆抬升 ⟹ `em(p) − em(n) ∝ Δ(Σq²)` ⟹ 若 `Δ(Σq²) > 0` 則 **p 較重** ⟹ `em(n) − em(p) < 0`。
  W3 若自能通道之號與所需相反 ⟹ 該通道**不可能**補缺口（不論係數大小、不論 `g(0)` 多大）。
  W4 若兩體上限 < 需求 ⟹ 兩體亦不足 ⟹ **兩通道皆無力 ⟹ W1 永久阻塞**（唔係「邊緣」）。
  W5 不引入新參數、不改 M10 任何結論；本檔只交「能否」之判定。
輸出：`data/56_m10_w1_channels.json`
"""
import json

ALPHA = 1 / 137.035999084
HBARC = 197.3269804
A0 = 0.384184
G0 = 2581 / 352            # g(0)，格點自能權重（有限）
NEED_FLAG, NEED_GAMMA = 1.2667, 1.3963

print("=" * 100)
print("M5 docs/54 — M10 W1：兩條通道能否補 em(n−p) 之缺口")
print("=" * 100)

q2 = {"u": 4 / 9, "d": 1 / 9}
sq_p = 2 * q2["u"] + q2["d"]      # uud
sq_n = q2["u"] + 2 * q2["d"]      # udd
print("\n[W2] Σq²（質子 uud vs 中子 udd）")
print(f"  Σq²(p) = 2×{q2['u']:.4f} + {q2['d']:.4f} = {sq_p:.6f}")
print(f"  Σq²(n) = {q2['u']:.4f} + 2×{q2['d']:.4f} = {sq_n:.6f}")
dsq = sq_p - sq_n
print(f"  Δ(Σq²) ≡ p − n = {dsq:+.6f} ⟹ 自能通道令 **p 較重** ⟹ `em(n)−em(p) < 0`")
print(f"  需求：`em(n) − em(p) = +{NEED_FLAG}` ⟹ [W3] 號相反 ⟹ **自能通道不可能補缺口**（與係數、g(0) 無關）")
amp = ALPHA * HBARC / A0
print(f"  （附量級：若硬用自能，(α/2)(ħc/a₀)|Δ|g(0) = {0.5*amp*abs(dsq)*G0:.4f} MeV —— 遠超需求，但號錯 ⟹ 不可用）")

print("\n[W1] 兩體上限（`docs/43`／`docs/37`）")
amax = amp * (1 / 3)
print(f"  α ħc/a₀ = {amp:.4f} MeV；A_max = α ħc/(3 a₀) = **{amax:.4f} MeV**")
for nm, need in (("FLAG 中性", NEED_FLAG), ("γ 比值", NEED_GAMMA)):
    print(f"  vs {nm} {need:.4f} ⟹ 差 {100*(need/amax-1):+.2f}%（缺口 {need-amax:.4f} MeV）")

print("\n[W4] 判")
print(f"  兩體：{amax:.4f} < {NEED_FLAG:.4f} ⟹ 不足 1.38%")
print("  自能：號相反 ⟹ 排除")
print("  ⟹ **兩通道皆無力 ⟹ W1 永久阻塞**（非『邊緣可補』）；M10 Tier II 之解鎖唔可能出自譜側 EM 通道。")
print("\n[W5] 附：此結論與 M5 docs/49 之『q² 自能對 π 恰零』同源 —— 自能通道在介子（號錯/零）與重子（號錯）皆被排除。")

out = {"sumq2_p": sq_p, "sumq2_n": sq_n, "delta_q2": dsq,
       "A_max": amax, "need_flag": NEED_FLAG, "need_gamma": NEED_GAMMA,
       "gap_pct": 100 * (NEED_FLAG / amax - 1),
       "self_channel_sign_ok": dsq < 0, "verdict": "W1 blocked by both channels"}
with open("data/56_m10_w1_channels.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/56_m10_w1_channels.json")