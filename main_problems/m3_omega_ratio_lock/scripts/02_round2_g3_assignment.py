#!/usr/bin/env python3
"""M3 scripts/02 — 第二輪：G3 賦值結構（閉式、剛性、暗份額合法性）

核心閉式（本輪推導）:
    Ω = ( (B−1) : B² : N_d−B+1 ) / (N_d + B²)
    B   = γ_max = 4   (Nyquist 定理, 與鏈長無關)
    N_d = 44          (88-site 鏈半譜模數)
驗證: B=4, N_d=44 ⟹ (3:16:41)/60。

另核:
  (i)  60 份預算 ≡ 90 譜預算 × 2/3 嘅精確性（Fraction 級）
  (ii) B 剛性：B=3/5/6 對 Planck 嘅 z
  (iii) N_d 剛性：N_d=44 vs 45（60 vs 61 份）——數據釘唔釘得住?
  (iv) 暗份額候選枚舉（gate 派生量, 排除雙重計數）
"""
from fractions import Fraction as F

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_d = 4, 44
obs = (0.0493, 0.2645, 0.6847)          # Ω_b, Ω_dm, Ω_Λ (Planck 2018, O215 docs/04)
sig = (0.0003, 0.0050, 0.0074)

def omega(b, nd):
    tot = nd + b*b
    return ((b-1)/tot, b*b/tot, (nd-b+1)/tot), tot

def chi2(t):
    return sum(((v-o)/s)**2 for v, o, s in zip(t, obs, sig))

# ── §1 閉式與精確性 ─────────────────────────────────────────────────
print("§1 閉式 Ω = ((B−1), B², N_d−B+1)/(N_d+B²)")
om, tot = omega(B, N_d)
print(f"  B={B}, N_d={N_d} ⟹ ({B-1}:{B*B}:{N_d-B+1})/{tot}")
chk("B² + B − 1 = 19（19 之單多項式身份）", B*B + B - 1 == 19, f"{B*B+B-1}")
chk("N_d − B + 1 = 41", N_d - B + 1 == 41, f"{N_d-B+1}")
chk("(B−1)+B²+(N_d−B+1) = N_d+B²", (B-1)+B*B+(N_d-B+1) == N_d+B*B, f"{tot}")
chk("41/60 ≡ 61.5/90 (EXACT, Fraction)",
    F(41,60) == F(F(123,2), 90), f"{F(41,60)} = {F(123,2)//F(1)}/…")
chk("19/60 ≡ 28.5/90 (EXACT, Fraction)", F(19,60) == F(F(57,2), 90))
print(f"  χ²(60 份) = {chi2(om):.3f}   z = " +
      "(" + ",".join(f"{((v-o)/s):+.2f}" for v,o,s in zip(om,obs,sig)) + ")")

# ── §2 B 剛性 ───────────────────────────────────────────────────────
print("\n§2 B 剛性（N_d=44 固定）")
for b in (3, 4, 5, 6):
    o, t = omega(b, N_d)
    z = [ (v-oo)/s for v,oo,s in zip(o,obs,sig) ]
    print(f"  B={b}: ({b-1}:{b*b}:{N_d-b+1})/{t} = " +
          "(" + ",".join(f"{v:.4f}" for v in o) + ")   z=(" +
          ",".join(f"{x:+.1f}" for x in z) + f")  χ²={chi2(o):.1f}")
print("  ⟹ B≠4 全部 >14σ 排除：B 剛性成立（B 由 γ_max=4 釘死）")

# ── §3 N_d 剛性（60 vs 61 份） ──────────────────────────────────────
print("\n§3 N_d 剛性 / 分母張力（B=4 固定）")
for nd in (42, 43, 44, 45, 46):
    o, t = omega(B, nd)
    print(f"  N_d={nd}: 分母 {t}   χ²={chi2(o):7.3f}   Ω_b={o[0]:.5f} (obs .0493)")
print("  ★ 數據偏好 N_d=45（分母 61）而非 44（60）：χ² 0.6 vs 5.7")
print("  ⟹ 已知 Ω_b +2.33σ 張力 = 【分母張力】, 唔係分子張力")

# ── §4 暗份額候選枚舉（gate 派生量, 排除雙重計數） ──────────────────
print("\n§4 暗份額 H 候選（gate 派生量; 需 ≠ 可見 3、≠ 未投影 41）")
cands = {
    "κ⁻¹ = B²（門檻倒數）":       B**2,
    "B²−1":                       B**2 - 1,
    "B(B−1) = N_gen·B":           B*(B-1),
    "N_d−B²":                     N_d - B**2,
    "γ_max = B":                  B,
    "Σγ_even−Σγ_odd (=2)":        2,
    "K_tree = M":                 33,
    "隱藏模數 = N_d−N_gen (=41)":  N_d - 3,
}
ok = []
for k, v in cands.items():
    clash = (v == 3) or (v == 41)
    total = 3 + v + 41
    o = (3/total, v/total, 41/total)
    recog = total in (60, 44, 90, 123, 178)
    print(f"  {k:26} H={v:3d}  total={total:3d}  {'✗ 撞 3/41' if clash else '○'}  "
          f"{'★ total∈結構數' if recog else ''}  χ²={chi2(o):8.1f}")
    if not clash and recog and chi2(o) < 6:
        ok.append((k, v))
print(f"  ⟹ 同時滿足（唔撞 / total 係結構數 / χ²<6）嘅候選: {ok}")
print("  判定: 只有 κ⁻¹=B² 過濾；屬【排除法】而非獨立定理（G3 仍未由第一性原理導出）")

# ── §5 張力歸因：比值 vs 標度 ────────────────────────────────────────
print("\n§5 張力歸因（比值通道 vs 標度通道）")
r_struct = F(B*B, B-1)                       # 16/3
r_obs    = obs[1]/obs[0]
s_r = r_obs * ((sig[1]/obs[1])**2 + (sig[0]/obs[0])**2) ** 0.5
print(f"  結構比 Ω_dm/Ω_b = B²/(B−1) = {r_struct} = {float(r_struct):.4f}")
print(f"  觀測比          = {r_obs:.4f} ± {s_r:.4f}   z = {(float(r_struct)-r_obs)/s_r:+.2f}σ")
D_need = (B-1)/obs[0]
print(f"  Ω_b=3/D 反解所需分母 D = {D_need:.2f}（結構給 60 ⟹ 差 {D_need-60:+.2f}）")
print("  ⟹ 比值通道 −0.3σ 命中；張力全部落在標度（分母 60 vs 60.85）")
print("     = 等價於 N_d=45（§3）⟹ 唔可以靠 Ω 數據釘死 N_d")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")