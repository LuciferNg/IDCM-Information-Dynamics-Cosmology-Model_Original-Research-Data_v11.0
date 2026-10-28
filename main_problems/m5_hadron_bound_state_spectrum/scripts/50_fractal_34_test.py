"""M5 docs/48 配套 — 「3/4」係唔係表面積／體積之分形律？
================================================================================
兩個層面分開判：

【一｜3/4 在 IDCM 之精確結構身份】（可否作為結構常數？）
  F1 b-ary 樹之 leaves/total：`L_n / T_n = b^n / ((b^{n+1}−1)/(b−1))`，`b = 4 = N_gen+1` ⟹ `n→∞ : 3/4`。
  F2 3-of-4 Cantor 子集之 Hausdorff 維度：`d_H = ln3/ln4`。
  F3 `45/60 = (N_d+1)/|I|`（`|I| = 60` 為二十面體旋轉群階）。
  ⟹ 三個獨立身份都給 3/4（或含 3、4）—— 故 3/4 **確有**結構來源，但三者互不等價。

【二｜3/4 能否解釋 M5 之觀測（表面/體積律之否證）】
  G1 若抑制因子 `s = A_B/A*` 源自樹之 volume/surface 比 ⟹ `s = T_n/L_n ∈ [1, b/(b−1)] = [1, 4/3]`。
     觀測 `s = 1/14.35` ⟹ 判：若 `s` 落在 [1,4/3] 外 ⟹ 樸素表面/體積律**不可能**。
  G2 整數冪檢：`(4/3)^k`（k = 1..12）中是否有 `≈ 14.35`（10% 容差）？
  G3 分形維度擬合：若 `A_f ∝ (r_h/a₀)^{d}`（`r_h = ħc/m_h`），由 B 定出 `d`，再**盲測** K 與 D。
     判：若同一 `d` 對 K 之預言偏離 > 2× ⟹ 無單一分形維度 ⟹ 分形律否證。
  G4 二值 UV/IR 判準（替代假說）：`r_h/a₀ < 1/8` ⟹ UV 族。報各族之 `r_h/a₀`。
  G5 不得以「3 族通過」之計數反推結構；計數 = 數據事實，非結構量。
輸出：`data/50_fractal_34_test.json`
"""
import json
import math

B = 4
print("=" * 100)
print("M5 docs/48 — 「3/4」分形律之判定")
print("=" * 100)

print("\n[F1] b = 4 之 leaves/total（b-ary 樹）")
for n in range(1, 8):
    L, T = B ** n, (B ** (n + 1) - 1) // (B - 1)
    print(f"  n = {n}: L/T = {L}/{T} = {L/T:.6f}")
print(f"  n→∞ : 3/4 = {3/4:.6f} ⟹ 精確極限")
print(f"\n[F2] d_H = ln3/ln4 = {math.log(3)/math.log(4):.6f}")
print(f"[F3] 45/60 = {45/60:.6f}（(N_d+1)/|I|）")

print("\n" + "-" * 100)
print("[G1] 樸素表面/體積律：s = T_n/L_n ∈ [1, 4/3]")
s_obs = 0.3100 / 4.4492
print(f"  觀測 s = {s_obs:.5f}（1/{1/s_obs:.2f}）")
print(f"  上界 b/(b−1) = {B/(B-1):.4f} ⟹ "
      f"{'🔴 s 遠在上界外 ⟹ 樸素表面/體積律不可能' if s_obs < 1 else '在界內'}")

print("\n[G2] 整數冪 (4/3)^k 是否有 ≈ 1/s = {:.3f}".format(1/s_obs))
hits2 = []
for k in range(1, 13):
    v = (B / (B - 1)) ** k
    dev = 100 * (v / (1 / s_obs) - 1)
    ok = abs(dev) <= 10
    if ok:
        hits2.append(k)
    print(f"  k = {k:2d}: {v:9.4f}（偏離 {dev:+8.1f}%）{'✅' if ok else ''}")
print(f"  ⟹ 命中 k：{hits2 if hits2 else '無'}")

print("\n[G3] 分形維度擬合：A_f ∝ (r_h/a₀)^d，由 B 定 d，盲測 K、D")
HB, A0 = 0.1973269804, 0.384184
MH = {"K": 0.495, "D": 1.27, "B": 4.18}
AF = {"K": 3.9340, "D": 4.8200, "B": 0.3100}
ratio = AF["B"] / AF["D"]
d_fit = math.log(ratio) / math.log(MH["D"] / MH["B"])
print(f"  (r_h/a₀) 之比 B/D = {MH['D']/MH['B']:.5f}；A 之比 = {ratio:.5f} ⟹ d = {d_fit:.4f}")
pred_K = AF["D"] * (MH["D"] / MH["K"]) ** d_fit
print(f"  盲測 K：預言 A_K = {pred_K:.4f} vs 觀測 {AF['K']:.4f} ⟹ 偏離 ×{pred_K/AF['K']:.2f}")
print(f"  ⟹ {'🔴 分形律否證（同一 d 不能容納 K）' if abs(pred_K/AF['K']-1) > 1.0 else '與 K 相容'}")

print("\n[G4] 替代：二值 UV/IR 判準 r_h/a₀ < 1/8")
for f, m in MH.items():
    print(f"  {f}: r_h/a₀ = {HB/(m*A0):.4f} ⟹ {'UV' if HB/(m*A0) < 0.125 else 'IR'}")
print("  （π：m 取 0.30 ⟹ 0.5843 ⟹ IR）")
print("\n[G5] 紀律：3/4 係 4 個數據點之計數，非結構量；不得以計數反推結構。")

out = {"tree_limit": 0.75, "d_H": math.log(3) / math.log(4), "s_obs": s_obs,
       "naive_sv_bound": B / (B - 1), "power_hits": hits2,
       "d_fit_from_B": d_fit, "A_K_pred": pred_K, "A_K_obs": AF["K"]}
with open("data/50_fractal_34_test.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/50_fractal_34_test.json")