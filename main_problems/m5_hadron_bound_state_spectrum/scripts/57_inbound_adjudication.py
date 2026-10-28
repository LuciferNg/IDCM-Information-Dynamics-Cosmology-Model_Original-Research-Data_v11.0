"""M5 docs/55 配套 — 三條入站裁定：M11 分支裁決（w=|v|²）、M6 docs/01、M10 v1/v2 並報
================================================================================
觸發：register 三條未處理入站
  M11（2026-10-28 22:52）：判 M5 `scripts/54` 之 `A0=0.384184` 係「v2 支污染」，要求統一至「正典 0.377422」；
      並附其 k₀^crit 交叉值。
  M6（2026-10-27 01:03）：`docs/01` prose 寫「g 遞增至 d=N/2」與閉式矛盾，要求 M5 自行更正。
  M10（2026-10-27 18:13 + 10-28 10:37）：要求 M5「確認是否已收 n−p 通道翻轉，及 **v1/v2 並報之規格**」。

判準（事前寫死）：
  N1 機率檢定（與勢能無關）：層基 `|d⟩_norm = N_d^(−1/2)Σ_{i∈layer d}|i⟩` 係正交歸一 ⟹
     任何物理態 `Ψ = Σ_d v_d|d⟩_norm` 之層機率 = `|v_d|²`（`Σ_d|v_d|² = 1`）。
     以「全 53 點均勻態」作零假設檢定：物理層機率 = `N_d/53`；若 `|v_d|² = N_d/53` 而 `N_d|v_d|² = N_d²/53` 則
     `Σ|v|² = 1`（唯一歸一者）⟹ **w = |v|²**；若用 `N_d·v²` 則總和 `Σ N_d v_d² = 27.5 ≠ 1` ⟹ 非機率。
  N2 M11 之「正典」須有結構依據，若只係追溯至 `docs/13`（v1 記錄）⟹ 屬「歷來如此」而非證明（`docs/27` §六.5）。
  N3 `A_max = k₀/3`，`k₀ = α ħc/a₀`；兩支並報。若 `A_max(v1) ≥ 1.2667 ≥ A_max(v2)` ⟹
     通道之生／死**完全由約定裁決決定**，須明文寫出。
  N4 `g(d)` 在 `0 ≤ d ≤ N/2` 之單調性以精確 `Fraction` 判（逐點差分），並定位整數零交叉。
  N5 不改 M11／M6／M10 任何檔案；只出判定。
輸出：`data/57_inbound_adjudication.json`
"""
import json
from fractions import Fraction as F

ALPHA = 1 / 137.035999084
HBARC = 197.3269804
A_V1 = 0.37742216816391858835   # M5 docs/13 記錄（w = N_d·v² 約定）
A_V2 = 0.384184                 # M5 docs/27 裁定（w = |v|²）
NEED = 1.2667
N = 88

print("=" * 100)
print("M5 docs/55 — 三條入站裁定")
print("=" * 100)

print("\n[N1] weight 約定：機率檢定（與勢能無關）")
layers = [1, 4, 12, 36]          # 4-ary 樹深度 3
tot = sum(layers)
print(f"  4-ary 樹層大小 {layers}，總節點 {tot}")
print(f"  均勻態物理層機率 N_d/{tot} = {[f'{l}/{tot}' for l in layers]}")
v2 = [l / tot for l in layers]           # |v_d|²
print(f"  若 v_d = √(N_d/{tot})：Σ|v_d|² = {sum(v2):.6f} ✓（＝1）")
print(f"  若用 N_d·|v_d|²：Σ N_d|v_d|² = {sum(l * p for l, p in zip(layers, v2)):.6f} ✗（≠1，且等於加權均值 {sum(l*l for l in layers)/tot:.6f}）")
print("  ⟹ 層基正交歸一 ⟹ 層機率 = |v_d|² ⟹ **w = |v|²（v2 支）**；N_d 加權係把層基當 per-site 態 ⟹ 雙重計數")
print("  ⟹ [N2] M11 之「正典 0.377422」追溯自 docs/13（v1 記錄）⟹ 屬『歷來如此』；")
print("     M5 已於 docs/27（scripts/29）以樹投影 8.9e−16 + 生產 H 之 OFF=[−√B,−√(B−1)…] 作**結構證明** ⟹ v2 為正典")

print("\n[N3] A_max 兩支並報（M10 明確要求之規格）")
rows = []
for tag, a in (("v1（w = N_d v²，已作廢）", A_V1), ("v2（w = |v|²，正典）", A_V2)):
    k0 = ALPHA * HBARC / a
    amax = k0 / 3
    gap = 100 * (NEED / amax - 1)
    rows.append({"branch": tag, "a0_fm": a, "hbarc_over_a0_MeV": HBARC / a, "k0_MeV": k0,
                 "A_max_MeV": amax, "gap_vs_need_pct": gap, "blocked": amax < NEED})
    print(f"  {tag}: a₀ = {a:.17f} fm ⟹ ħc/a₀ = {HBARC/a:9.4f} MeV ⟹ k₀ = {k0:.10f} MeV ⟹ "
          f"A_max = {amax:.4f} MeV ⟹ vs {NEED}: {gap:+.2f}% ⟹ {'超上限' if amax < NEED else '在上限內'}")
ratio = (ALPHA*HBARC/A_V1) / (ALPHA*HBARC/A_V2)
print(f"  兩支比 = {ratio:.6f}（＝a₀ 之比 {A_V2/A_V1:.6f}）⟹ 約定差 {100*(ratio-1):.2f}%")
print(f"  ★ A_max(v1) = {ALPHA*HBARC/A_V1/3:.4f} ≥ {NEED} ≥ {ALPHA*HBARC/A_V2/3:.4f} = A_max(v2)")
print("  ⟹ **n−p 通道之生／死完全由約定裁決決定**：v1 下在上限內（餘裕 0.40%）、v2 下超 1.39%")
print("  ⟹ M5 正典 = v2（N1 結構證明）⟹ W1 阻塞；但須明文標示此依賴（本檔即此規格）")

print("\n[N3b] docs/45 之『0.514 GeV』是否污染")
print(f"  ħc/a₀(v2) = {HBARC/A_V2/1000:.6f} GeV ⟹ 記作 0.514 GeV ✓ 與 v2 一致（非污染）；ħc/a₀(v1) = {HBARC/A_V1/1000:.6f} GeV ⟹ 0.523")
print("  ⟹ M11 之 ② 屬誤判（該數字係 v2 分支之正確值，而 v2 為 M5 正典）")

print("\n[N4] docs/01 g(d) 單調性（M6 之 docs/13 主張）")
def g(d, n=N):
    return F((n * n - 1) - 6 * d * (n - d), 12 * n)
vals = [g(d) for d in range(N // 2 + 1)]
strict_dec = all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
print(f"  g(0) = {vals[0]} = {float(vals[0]):+.6f}；g(44) = {vals[-1]} = {float(vals[-1]):+.6f}")
print(f"  0 ≤ d ≤ 44 嚴格遞減 = {strict_dec} ⟹ M6 判對、docs/01 原 prose『遞增』錯")
lastpos = max(d for d in range(N // 2 + 1) if vals[d] > 0)
print(f"  整數零交叉：g({lastpos}) = {float(vals[lastpos]):+.6f} > 0，g({lastpos+1}) = {float(vals[lastpos+1]):+.6f} < 0 ⟹ 最後正值 d = {lastpos}（與 M6 報 d=18 一致）")
print("  ⟹ M5 docs/01 已於 2026-10-27 加勘誤＋附註 A1（原文保留）⟹ 此項已處理，本檔只作核實")

out = {"branches": rows, "branch_ratio": ratio, "g_strictly_decreasing": strict_dec,
       "g_last_positive_d": lastpos, "w_convention": "w=|v|^2 (v2)",
       "m11_rejected": True, "m6_docs01_status": "fixed 2026-10-27"}
with open("data/57_inbound_adjudication.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/57_inbound_adjudication.json")