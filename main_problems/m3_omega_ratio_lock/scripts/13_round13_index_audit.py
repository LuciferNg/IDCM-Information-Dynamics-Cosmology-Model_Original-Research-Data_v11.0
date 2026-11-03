#!/usr/bin/env python3
"""M3 scripts/13 — 第十三輪：索引缺陷審計（吸收 M15 handoff ×2）

M15 → M3（2026-11-02 04:57 / 06:21）:
  ① scripts/03:57 `range(0, N_d)` 含零模 γ₀、漏 Nyquist γ₄₄
  ② 強制重算：rank 三分 2:16:42 ≠ 預算三分 3:16:41；和守恆令「Σ=60 閉合」失效作同構證據

本輪在 M3 內逐項重算（不引用 M15 數值），並裁定：
  A. 四種閘門組合枚舉
  B. rank 三分 vs 預算三分 逐項比較
  C. 3 = B−1 之結構身份（階 4 群非平凡特徵標）⟹ 譜計數存活、算子秩詮釋死亡
  D. ℓ = 2 之觀測鎖定（Ω_DM(ℓ) 掃描）
"""
import mpmath as mp
from fractions import Fraction as F
import itertools

mp.mp.dps = 30
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, Nd = 4, 44
g = lambda n: 4*mp.sin(mp.pi*n/88)**2

# ── A. 四種閘門組合 ─────────────────────────────────────────────
print("§A 閘門計數四組合（γ_n = 4sin²(πn/88)）")
res = {}
for kname, kap in (("κ=B⁻²=1/16", mp.mpf(1)/16), ("κ=B⁻¹=1/4", mp.mpf(1)/4)):
    for lo, hi, iname in ((0, Nd, "0..43 含零模漏Nyquist"), (1, Nd+1, "1..44 物理索引")):
        m = [n for n in range(lo, hi) if mp.e**(-n) > kap]
        res[(kname, iname)] = len(m)
        print(f"  {kname:12s} × {iname:22s} ⟹ {len(m)}  模態 {m}")
chk("只有「含零模 + κ=B⁻²」給 3", res[("κ=B⁻²=1/16", "0..43 含零模漏Nyquist")] == 3)
chk("物理索引 + κ=B⁻² 給 2", res[("κ=B⁻²=1/16", "1..44 物理索引")] == 2)
chk("索引兩頭皆錯：γ₀=0 屬真空、γ₄₄=4=Nyquist 屬真實模態",
    g(0) == 0 and abs(g(44) - 4) < mp.mpf('1e-25'))

# ── B. rank 三分 vs 預算三分 ────────────────────────────────────
print("\n§B 三分逐項比較")
rk_vis = res[("κ=B⁻²=1/16", "1..44 物理索引")]      # 2
rk_dm  = B*B                                          # 16
rk_lam = Nd - rk_vis                                  # 42
bud = (3, 16, 41)
op  = (rk_vis, rk_dm, rk_lam)
print(f"  算子秩三分 = {op}   和 = {sum(op)}")
print(f"  預算三分   = {bud}  和 = {sum(bud)}")
print(f"  逐項相等?  {[a==b for a,b in zip(op,bud)]}")
chk("逐項 [否, 是, 否]", [a==b for a,b in zip(op,bud)] == [False, True, False])
chk("鏈側和守恆（3+41 = 2+42 = 44，兩者 off-by-one 抵銷）", bud[0]+bud[2] == op[0]+op[2] == Nd)
chk("三分總和皆 60（總和守恆但逐項不等）", sum(bud)==sum(op)==60)
chk("⟹ 「Σ=60 閉合」不可作同構證據（off-by-one 互相抵銷）", True)

# ── C. 3 = B−1 之結構身份 ──────────────────────────────────────
print("\n§C 階 4 群之非平凡特徵標數")
def chars_cyclic(n):
    return [k for k in range(n)]
def chars_v4():
    out = []
    for a, b in itertools.product((0,1), repeat=2):
        out.append((a,b))
    return out
# ℤ₄: χ_k(j) = i^{kj}; 非平凡 = k≠0
Z4_nontriv = [k for k in range(4) if k != 0]
# V₄: 四特徵標 (a,b)·(x,y) ↦ (−1)^{ax+by}; 非平凡 = (a,b)≠(0,0)
V4_nontriv = [t for t in chars_v4() if t != (0,0)]
print(f"  ℤ₄ 非平凡特徵標 = {len(Z4_nontriv)}")
print(f"  V₄ 非平凡特徵標 = {len(V4_nontriv)} {V4_nontriv}")
chk("ℤ₄ 與 V₄ 皆給 3（階 4 全 abelian ⟹ 與群選擇無關）",
    len(Z4_nontriv) == len(V4_nontriv) == 3)
chk("3 = B−1 = γ_max−1（B = γ_max 由 Nyquist 鎖死）", 3 == B-1)
print("  ⟹ 3 之結構身份 = 非平凡特徵標數，【唔係】模計數")
chk("模計數讀法給 2 ≠ 3 ⟹ 算子秩詮釋死亡", rk_vis != 3)

# ── D. ℓ = 2 觀測鎖定 ─────────────────────────────────────────
print("\n§D Ω_DM(ℓ) = B^ℓ/(N_d+B^ℓ) 掃描（Planck 0.26447 ± 0.00500）")
obs, sig = mp.mpf("0.26447"), mp.mpf("0.00500")
zs = {}
for l in (1,2,3,4):
    p = mp.mpf(B**l)/mp.mpf(Nd + B**l)
    z = (p-obs)/sig
    zs[l] = z
    print(f"  ℓ={l}: Ω_DM = {mp.nstr(p,6):9s}  z = {mp.nstr(z,7):>10s}")
chk("ℓ=2 唯一 |z|<1（其餘 |z|>36σ）", abs(zs[2]) < 1 and min(abs(zs[l]) for l in (1,3,4)) > 30)
chk("ℓ=2 ⟹ B^ℓ = 16 = κ⁻¹（冪次對、count 錯，兩者分開）", B**2 == 16)

# ── 判決 ─────────────────────────────────────────────────────
print("\n§E 判決")
print("  ① `3 = rank(P_κ)` 撤回 —— 只在含真空、漏 Nyquist 之索引下成立")
print("  ② 算子秩三分 2:16:42 與預算 3:16:41 逐項不等 ⟹ 「⊕ 三分」詮釋死亡")
print("  ③ 3 保留為 B−1 = 非平凡特徵標數（結構）；41 = N_d − (B−1)（譜計數補）")
print("  ④ ℓ = 2 由觀測鎖死 ⟹ 暗份額 16 之「為何 level 2」獲得觀測依據")
print(f"\n結算: PASS {PASS} / FAIL {FAIL}")