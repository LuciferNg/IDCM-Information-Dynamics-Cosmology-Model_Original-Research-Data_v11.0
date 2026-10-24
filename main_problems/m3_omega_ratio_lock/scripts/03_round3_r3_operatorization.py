#!/usr/bin/env python3
"""M3 scripts/03 — 第三輪：R3 算子化（三份額升格為顯式算子的 rank/dim）

核心認同（本輪）:
    K_chain  ≅ Laplacian of cycle C_88  (L = D − A)
    spec(L) = { 4 sin²(πk/88) : k = 0..87 }  ⟹ 半譜 k=1..44 = 44 modes = N_d
    tr L = 2N = 176 ,  γ_max = γ_44 = 4 = B
    K_tree  = 4-ary tree (branching B) ; level-k 節點數 = B^k
    gate    : κ = B⁻² ; 可見模態數 = #{n≥0 : e^{−n} > κ} = 3

三份額:
    N_gen       = rank(P_κ)          (chain 上嘅閘門投影, 秩 3)
    B²          = rank(P_level2)     (tree level-2 節點數 = 16)
    N_d − N_gen = rank(1 − P_κ)      (chain 未投影, 秩 41)
    total       = N_d + B² = 60
"""
import numpy as np
from mpmath import mp, mpf, sin, pi, e, log
mp.dps = 50

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

B, N_d, N = 4, 44, 88

# ── §1 K_chain = Laplacian(C_88)：譜 = γ_n ───────────────────────────
print("§1 K_chain ≅ L(C_88)")
A = np.zeros((N, N), dtype=int)
for i in range(N):
    A[i, (i+1) % N] = 1; A[(i+1) % N, i] = 1
D = np.diag(A.sum(axis=1))
L = D - A
ev = np.linalg.eigvalsh(L.astype(float))
chk("tr L = 2N = 176", int(np.trace(L)) == 2*N, f"{int(np.trace(L))}")
g = np.sort([4*np.sin(np.pi*k/N)**2 for k in range(N)])
chk("spec(L) = {4sin²(πk/88)}", float(np.max(np.abs(np.sort(ev) - g))) < 1e-10,
    f"max|Δ|={float(np.max(np.abs(np.sort(ev)-g))):.2e}")
chk("γ_max = γ₄₄ = 4 = B", abs(g[-1] - 4) < 1e-12, f"{g[-1]:.12f}")
chk("zero mode 唯一 (C_88 連通)", int(np.sum(ev < 1e-9)) == 1)

# ── §2 半譜求和 = (tr L + γ_max)/2（算子形式） ───────────────────────
print("\n§2 半譜求和（算子形式）")
gh = [4*sin(mpf(k)*pi/N)**2 for k in range(1, N//2 + 1)]
S_half = sum(gh)
chk("Σγ(half) = (tr L + γ_max)/2 = 90",
    abs(S_half - (mpf(2*N) + 4)/2) < mpf('1e-40'), f"{mp.nstr(S_half,12)}")
chk("Σγ(half) = 2N_d + 2（T-M1, B-普遍）", abs(S_half - (2*N_d + 2)) < mpf('1e-40'))
print(f"  ⚠️ O211 之「Σγ = B(N_d+1)/2 = {B*(N_d+1)//2}」與 2N_d+2 只喺 B=4 相等")
print(f"     （(B−4)(N_d+1)=0 ⟹ 屬 B=4 改寫, 非獨立律）")

# ── §3 閘門：κ = B⁻²，rank(P_κ) = 3 ─────────────────────────────────
print("\n§3 閘門投影 rank")
kappa = mpf(1)/B**2
vis = [n for n in range(0, N_d) if e**(-n) > kappa]
chk("κ = B⁻² = 1/16", kappa == mpf(1)/16, f"{mp.nstr(kappa,8)}")
chk("rank(P_κ) = #{e^{−n} > κ} = 3 = N_gen", len(vis) == 3, f"n={vis}")
chk("rank(1−P_κ) = N_d − 3 = 41", N_d - len(vis) == 41, f"{N_d-len(vis)}")
chk("⌈ln B²⌉ = 3（解析同值）", int(np.ceil(np.log(B*B))) == 3,
    f"ln16={np.log(16):.4f}")

# ── §4 K_tree level 計數：level-2 = B² = 16 ─────────────────────────
print("\n§4 K_tree level 計數")
level = [1]
for k in range(1, 3):
    level.append(level[-1]*B)
chk("level-0/1/2 = 1/4/16", level == [1, 4, 16], f"{level}")
chk("rank(P_level2) = B² = 16", level[2] == B**2, f"{level[2]}")
print(f"  ⟹ 可見 3 = 閘門（tree 閾值 ln B²）; 暗 16 = tree level-2; 真空 41 = chain 未投影")

# ── §5 組裝與可測比（chain-獨立性 = 真預言） ─────────────────────────
print("\n§5 組裝：60 = 3 + 16 + 41，且比值 chain-獨立")
obs = (0.0493, 0.2645, 0.6847); sig = (0.0003, 0.0050, 0.0074)
def omega(nd):
    tot = nd + B**2
    return ((B-1)/tot, B**2/tot, (nd-(B-1))/tot), tot
om, tot = omega(N_d)
chk("3 + 16 + 41 = 60", (len(vis)) + B**2 + (N_d-len(vis)) == 60, f"{tot}")
chk("份額 = rank 三元組", om == (0.05, 16/60, 41/60), f"{om}")
r = B**2/(B-1)
print(f"  可測比 Ω_dm/Ω_b = B²/(B−1) = {r:.4f}（**唔含 N_d**）")
for nd in (42, 43, 44, 45, 46):
    o, t = omega(nd)
    print(f"    N_d={nd}: 分母 {t}  ratio={o[1]/o[0]:.4f}  Ω_b={o[0]:.5f}")
chk("ratio 對 N_d 完全不變", all(abs(omega(nd)[0][1]/omega(nd)[0][0] - r) < 1e-15
                            for nd in range(40, 50)))
print("  ⟹ Ω_dm/Ω_b 係純 tree 量（chain-independent）: 預言級")
print("  ⟹ 絕對標度 (Ω_b=3/60) 需要 N_d ⟹ 呢個就係 Ω_b 張力所在")

# ── §6 B 掃描：比值單獨已可釘死 B ───────────────────────────────────
print("\n§6 B 掃描（只用比值通道，唔用 N_d）")
r_obs = obs[1]/obs[0]
s_r = r_obs*((sig[1]/obs[1])**2 + (sig[0]/obs[0])**2)**0.5
for b in (3, 4, 5, 6):
    ng = int(np.ceil(np.log(b*b)))
    rb = b*b/ng
    print(f"  B={b}: N_gen=⌈ln B²⌉={ng},  ratio={b*b}/{ng}={rb:.4f},  "
          f"z={(rb-r_obs)/s_r:+.1f}σ")
print(f"  觀測 ratio = {r_obs:.4f} ± {s_r:.4f}")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")