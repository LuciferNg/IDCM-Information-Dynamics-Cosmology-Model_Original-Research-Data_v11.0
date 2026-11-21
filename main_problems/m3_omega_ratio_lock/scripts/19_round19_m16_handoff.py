#!/usr/bin/env python3
"""M3 scripts/19 — 第十九輪：M16 handoff 處理（R5 理由更換 ／ MF2 閉合宣告）

兩條入向（M16）:
  ① 2026-11-03 20:22  R5 判詞：60=|A₅| 結論保留、理由改用 M16 定理 C(G) ⟺ G≅A₅
  ② 2026-11-07 14:53  補充：MF2 閉合（perfect 由 Ω_b 測定 k=1.0142±0.0176）；M3 無 M16 引用故無需重算

本輪獨立驗證（不採納對方數字）:
  A. M16 算術引理：distinct dims ≥2、Σ=12、∏=|G| ⟹ 唯一集 = {3,4,5}，故 |G|=60
  B. A₅ 之 dims = (1,3,3,4,5) ⟹ distinct 非平凡 = {3,4,5} ✓
  C. 「階 60 之 13 群中唯一」— 無 GAP／小群庫 ⟹ 本 M 未能獨立驗證（標外部未證）
  D. k = 1.0142 之來源：k ≡ (3/60)/Ω_b(obs) ⟹ 定義式、非獨立量測
  E. σ_k：M16 用 0.0176；由 Ω_b 傳播 = σ_b/(1/20) = 0.0060（2.93× 差）
  F. M3 內 'perfect' / |G:G'| 引用 = 0 ⟹ 該輸入不屬本 M
"""
from itertools import combinations
from math import prod
import mpmath as mp

mp.mp.dps = 30
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

# ── A. 算術引理 ─────────────────────────────────────────────────
print("§A M16 算術引理：相異 dims ≥2、Σ = 12")
sets = [c for r in range(1, 7) for c in combinations(range(2, 12), r) if sum(c) == 12]
print("  集 =", sets)
prods = {c: prod(c) for c in sets}
print("  積 =", prods)
uniq = [c for c, p in prods.items() if p == 60]
chk("Σ=12 且 ∏=60 之集唯一 = {3,4,5}", uniq == [(3, 4, 5)])
chk("⟹ 該類群階被逼為 |G| = 60", prod((3, 4, 5)) == 60)

# ── B. A₅ 之維數 ───────────────────────────────────────────────
print("\n§B A₅ 不可約維數")
dims = [1, 3, 3, 4, 5]
distinct = sorted({d for d in dims if d != 1})
chk("distinct 非平凡 dims = {3,4,5}", distinct == [3, 4, 5])
chk("∏ distinct = 60 = |A₅| ；Σ = 12", prod(distinct) == 60 and sum(distinct) == 12)

# ── C. 階 60 之 13 群 ─────────────────────────────────────────
print("\n§C 「階 60 之 13 群中唯一」")
print("  無 GAP／小群庫 ⟹ 本 M 未能獨立驗證；標「外部未證」，不作 🟢 依據")
chk("本 M 對該子命題之立場 = 未驗（誠實標記）", True, "external-unverified")

# ── D/E. k 之來源與 σ ─────────────────────────────────────────
print("\n§D/E MF2 之 k")
obs_b, sig_b = mp.mpf("0.0493"), mp.mpf("0.0003")
k = (mp.mpf(3)/60)/obs_b
print(f"  k ≡ (3/60)/Ω_b(obs) = {mp.nstr(k, 8)}  （M16 報 1.0142）")
chk("k 與 (3/60)/Ω_b 逐位相同 ⟹ 定義式非量測", abs(k - mp.mpf("1.0142")) < mp.mpf("1e-4"))
sig_prop = sig_b/(mp.mpf(1)/20)
z_prop = (k - 1)/sig_prop
z_m16 = (k - 1)/mp.mpf("0.0176")
print(f"  σ_k（由 Ω_b 傳播）= {mp.nstr(sig_prop, 5)} ；M16 用 0.0176 ⟹ 比 {mp.nstr(mp.mpf('0.0176')/sig_prop, 4)}×")
print(f"  z = {mp.nstr(z_prop, 4)}σ（傳播） vs {mp.nstr(z_m16, 4)}σ（M16）")
chk("M16 σ_k 非由 Ω_b 傳播（差 > 2×）", mp.mpf("0.0176")/sig_prop > 2)
chk("用 Ω_b 傳播 ⟹ z 即 R2 之 +2.33σ（同一張力重述）", abs(z_prop - mp.mpf("2.37")) < mp.mpf("0.1"))

# ── F. M3 有無 perfect 輸入 ───────────────────────────────────
print("\n§F M3 之依賴")
import subprocess
hits = [h for h in subprocess.run(["grep","-rl","--include=*.py","-e","perfect","-e","G:G","scripts"],
                                   capture_output=True, text=True).stdout.split()
        if "19_round19_m16_handoff.py" not in h]   # 排除本腳本自身（含該詞作說明）
chk("M3 內 'perfect'／|G:G'| 引用 = 0 ⟹ 不屬本 M 之輸入", not hits, f"{hits}")

# ── 判決 ──────────────────────────────────────────────────────
print("\n§判決")
print("  ① A：接受（算術引理，本 M 獨立重現）—— M16 之『60 = 3·4·5』在群選擇框架內係推論")
print("  ② B：接受（A₅ 滿足條件）")
print("  ③ C：標外部未證（階 60 唯一性未驗）")
print("  ④ D/E：🔴 MF2 之『由 Ω_b 測定』不成立 —— k 係 pred/obs 之定義式，非獨立量測")
print("  ⑤ F：M3 無 perfect 輸入 ⟹ 該輸入與本 M 閘門無關（M16 之歸因不適用）")
print("  ⑥ M3 之 §八 P3「|A₅| 與預算 60 不可互證」維持 —— 理由獲升級，身份仍未建立")
print(f"\n結算: PASS {PASS} / FAIL {FAIL}")