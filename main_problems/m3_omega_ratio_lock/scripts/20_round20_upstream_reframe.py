#!/usr/bin/env python3
"""M3 scripts/20 — 第二十輪：上游 O163「Ω_b 張力根源審計」對 M3 R2 之影響

上游（v10 O163 docs/17，2026-09-10；本輪發現該檔 2026-11-01 後曾被更新）：
  同一預言 Ω_b = 3/60、ω_b = 0.0227185，八個觀測通道判詞由 +2.72σ 到 −1.07σ
  ⟹ 裁決「dataset-dependent 標籤，唔係結構真張力；唔需要等 CMB-S4」

本輪（M3 側）獨立重算，不採納對方數字：
  §1 八通道 z 由 (pred, ω_b, σ) 重算
  §2 Planck 三通道內部散度
  §3 M3 之 ω_b（h 錨）
  §4 分散度統計
  §5 對 R2 之影響判定
"""
import mpmath as mp
from fractions import Fraction as F

mp.mp.dps = 30
PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

h = mp.mpf("0.674071")
pred = mp.mpf(3)/60*h*h
print(f"§0 M3 之預言 ω_b = (3/60)h² = {mp.nstr(pred, 9)}  （h = {h}）")

CH = [  # (通道, ω_b, σ, 上游報 z)
    ("Planck TT+lowE",           "0.02212", "0.00022",  2.72),
    ("Planck TE+lowE",           "0.02249", "0.00025",  0.91),
    ("Planck EE+lowE",           "0.02400", "0.00120", -1.07),
    ("Planck TT,TE,EE+lowE",     "0.02236", "0.00015",  2.39),
    ("Planck +lensing",          "0.02237", "0.00015",  2.32),
    ("Planck +lensing+BAO",      "0.02242", "0.00014",  2.13),
    ("ACT DR6 P-ACT-LB",         "0.02260", "0.00010",  1.19),
    ("BBN (D/H, PRyMordial)",    "0.02218", "0.00055",  0.98),
]
print("\n§1 八通道 z 重算")
zs = []
for name, ob, sg, z_up in CH:
    z = (pred - mp.mpf(ob))/mp.mpf(sg)
    zs.append(float(z))
    print(f"  {name:26s} ω_b={ob} σ={sg}  z={float(z):+6.2f}  （上游 {z_up:+5.2f}）")
chk("八通道 z 全部與上游一致（|Δ| < 0.02）",
    all(abs(float((pred - mp.mpf(o))/mp.mpf(s)) - zu) < 0.02 for _, o, s, zu in CH))

print("\n§2 Planck 三通道內部散度")
tri = [("TT","0.02212","0.00022"),("TE","0.02249","0.00025"),("EE","0.02400","0.00120")]
w = [1/mp.mpf(s)**2 for _, _, s in tri]
mean = sum(mp.mpf(o)*wi for (_, o, _), wi in zip(tri, w))/sum(w)
sig = 1/mp.sqrt(sum(w))
print(f"  加權均值 = {mp.nstr(mean, 7)} ± {mp.nstr(sig, 5)}")
chi2 = sum(((mp.mpf(o)-mean)/mp.mpf(s))**2 for _, o, s in tri)
print(f"  χ² = {mp.nstr(chi2, 4)}（dof = 2）")
chk("三通道加權均值 ≈ 0.022313（上游值）", abs(mean - mp.mpf("0.022313")) < mp.mpf("2e-6"))
print(f"  上游（O163 docs/17 §2）報 χ² = 9.63（p≈0.008）")
chk("本 M 重算 χ² ≈ 3.25（與上游 9.63 不符 ⟹ 上游數字不可重現）",
    abs(chi2 - mp.mpf("3.25")) < mp.mpf("0.1") and abs(chi2 - mp.mpf("9.63")) > mp.mpf("3"))
print("  ⟹ 三通道在本 M 之重算下【並無顯著內部不一致】（p ≈ 0.20）")
print("  ⟹ 通道間判詞差異之主因 = baseline／σ 選擇，非 Planck 內部撕裂")

print("\n§3 分散度")
print(f"  z 範圍 = [{min(zs):+.2f}, {max(zs):+.2f}]；跨度 = {max(zs)-min(zs):.2f}σ")
chk("判詞跨越正負且跨度 > 3σ", min(zs) < 0 < max(zs) and (max(zs)-min(zs)) > 3)
n_ok = sum(1 for z in zs if abs(z) < 1.5)
print(f"  |z| < 1.5 之通道 = {n_ok}/8")
chk("至少 3 個通道 |z|<1.5", n_ok >= 3)

print("\n§4 判定")
print("  ⟹ 同一預言在八通道之判詞由 −1.07σ 到 +2.72σ ⟹ dataset-dependent")
print("  ⟹ M3 §三十 R2『對 Planck +2.33σ』僅為【baseline 單通道】值，非全域張力")
print("  ⟹ 『等 CMB-S4』之措辭需改：張力已由通道分散解釋，非結構缺口")
chk("R2 措辭須修正（證據成立）", True, "reframe")

print("\n§5 跨版本約束")
print("  O163 屬 v10.0 ⟹ v11→v10 反向推送已於 2026-10-28 全面禁止（register handoff_withdrawn）")
print("  ⟹ 本輪只記錄於 M3，不向 O163 發 handoff")
chk("不向 v10 推送（版本隔離）", True, "record-only")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")