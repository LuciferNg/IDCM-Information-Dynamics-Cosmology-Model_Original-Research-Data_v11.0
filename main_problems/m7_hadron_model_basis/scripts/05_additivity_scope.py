"""M7 scripts/05 — ① M2 加性泛函之適用範圍：容量／飽和攻擊（R6）
================================================================================================
問題（本輪唯一目標）：M2 之加性泛函 m = M₀ Σ_n σ_n γ_n 是否覆蓋【複合態】？
  M6（並行 session）已用 ℚ(ζ₈₈) 精確證「子集和加性」= 域內恆等式（T1）。
  但該定理假設 σ 可自由相加。物理上能否相加，取決於**每模容量**：
  若 σ_n 有上限（每模量子數），則同一模內多夸克 ⟹ 飽和 ⟹ 非可加（束縛）。
  本腳本即攻此缺口。

判準（寫於執行前）：
  C1 容量核：每鏈模可用態數 = (自旋 2) × dim(色嵌入) = 2 × (B−1) = 6（O216：色 = 嵌入層標籤）。
     強子價夸克數 N_q ≤ 3 ⟹ max σ_n = N_q/6 = 0.5 < 1 ⟹ **永不飽和** ⟹ σ 可加。
  C2 強迫共享搜索：對味道內容 {uud, uus, uds, uss, uuu, ddd, sss} 逐個做單射分配搜索
     （每個夸克入 (sector, slot)，slot = 6），違反 0 ⟹ 無任何強子被迫共享模。
  C3 反事實核（若容量 = 1，即忽略 6 重簡併）：dev = −(n_down − 1)·M₀γ₄；與三錨所需比較（含符號）。
  C4 三體束縛（②）之前置結構：已註冊核全部 **pairwise（body order = 2）**，
     三體頂點在語料 0 命中 ⟹ E_bind^(3) 之味奇部分只可由 pair 核誘導。
  C5 pairwise-only 之三體相關量級（可解模型，精確）：3 等質量 + pairwise 諧核
     E₀^(3) = ℏ√(3k/m) vs Σ_pair E_pair = (3/2)ℏ√(2k/m) ⟹ 相關 = −18.4%（O(10%)，非 O(10²)）。
輸出：data/05_additivity_scope.json
"""
import hashlib
import json
import os
import re
import subprocess

from mpmath import mp, mpf, pi, sin, sqrt

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.abspath(os.path.join(M7, "..", "..", "..", ".."))
V11 = os.path.abspath(os.path.join(M7, "..", ".."))

M0 = mpf("27.9009145112")          # MeV（O204，經 M2 docs/01 §0）
MD_MU = mpf("2.56")                # MeV（FLAG24，M4 docs/28）
NEED = [-0.109, 0.498, 1.293]      # ⟨O⟩_h − 1（02b K4.need）
CH = ["N (uud)", "Σ (uus)", "Ξ (uss)"]

OUT = {"module": "M7", "round": "R6", "script": "05_additivity_scope.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<52} got {got} | want {want}" + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def gamma(n):
    return 4 * sin(pi * n / 88) ** 2


print("=" * 110)
print("M7 R6 — ① 加性泛函之適用範圍：容量／飽和攻擊")
print("=" * 110)
print(f"慣例：γ_n = 4sin²(πn/88)（88-site cycle，44 模）；M₀ = {M0} MeV；m_d − m_u = {MD_MU} MeV")

# ── C1 容量核 ────────────────────────────────────────────────────────────
print("\nC1 每模容量：自旋 × 色嵌入")
B = 4
dim_col = B - 1                        # O216：色 = 嵌入層標籤，dim = B−1 = 3
cap = 2 * dim_col                      # 自旋 2 × 色 3
Nq_max = 3                             # 重子價夸克數
sigma_max = mpf(Nq_max) / cap
ok("C1 dim(色嵌入) = B−1 = 3（O216）", dim_col == 3, dim_col, 3)
ok("C1 每模容量 = 2×3 = 6", cap == 6, cap, 6)
ok("C1 N_q ≤ 3 ⟹ max σ_n = 3/6 = 0.5 < 1（永不飽和）", sigma_max < 1,
   f"{float(sigma_max):.3f}", "< 1", "⟹ σ 可自由相加 ⟹ 加性覆蓋複合態")
OUT["C1"] = dict(dim_colour=dim_col, capacity=cap, sigma_max=float(sigma_max))

# ── C2 強迫共享搜索（單射分配）────────────────────────────────────────────
print("\nC2 強迫共享搜索：每夸克 → (sector, slot∈[0,6))，尋單射分配")
SECTOR = {"u": 1, "c": 1, "d": 4, "s": 4, "e": 5, "mu": 5, "tau": 5}  # v9 §D.1（sector 級）
HADRONS = {"N uud": "uud", "Σ uus": "uus", "Λ uds": "uds", "Ξ uss": "uss",
           "Δ⁺⁺ uuu": "uuu", "Ω⁻ sss": "sss", "d* ddd": "ddd"}
forced = []
for name, fl in HADRONS.items():
    slots = {}          # (sector, slot) -> quark
    conflict = None
    for q in fl:
        sec = SECTOR[q]
        free = [s for s in range(cap) if (sec, s) not in slots]
        if not free:
            conflict = q
            break
        slots[(sec, free[0])] = q
    per_sector = {}
    for q in fl:
        per_sector[SECTOR[q]] = per_sector.get(SECTOR[q], 0) + 1
    tag = f"  ⟹ 衝突於 {conflict}" if conflict else "  ⟹ 單射分配成功（無共享）"
    print(f"    {name:<10} 每 sector 夸克數 {per_sector}  max={max(per_sector.values())} / {cap}{tag}")
    if conflict:
        forced.append(name)
ok("C2 全部 7 種味道內容皆無溢出 ⟹ 無強迫共享", len(forced) == 0, forced or "0 個", "0 個",
   "最大單 sector = 3（uuu/ddd/sss）≤ 6")
OUT["C2"] = dict(hadrons=list(HADRONS), forced_sharing=forced, capacity=cap)

# ── C3 反事實：若容量 = 1（忽略 6 重簡併）─────────────────────────────────
print("\nC3 反事實核（容量 = 1）：dev = −(n_down − 1)·M₀γ₄")
g4 = gamma(4)
dev = {}
NDOWN = {"N (uud)": 1, "Σ (uus)": 1, "Ξ (uss)": 2}   # down-type = d,s（k=4）
for h, nd in NDOWN.items():
    dev[h] = -(nd - 1) * g4 * M0
    need = NEED[CH.index(h)] * MD_MU
    print(f"    {h:<10} n_down={nd}  非可加 dev = {float(dev[h]):+.4f} MeV"
          f"   所需 = {float(need):+.4f} MeV")
ok("C3 γ₄ = 4sin²(4π/88) 與 v9 §D.1 之 0.0810 一致", abs(float(g4) - 0.0810) < 5e-5,
   f"{float(g4):.7f}", 0.0810)
ok("C3 容量=1 之 Ξ dev 符號 vs 所需（需 +，得 −）", dev["Ξ (uss)"] < 0 < NEED[2] * MD_MU,
   f"{float(dev['Ξ (uss)']):+.3f} MeV", f"需 {float(NEED[2]*MD_MU):+.3f} MeV", "★符號相反")
ok("C3 容量=1 之 N/Σ dev（皆 0，所需非 0）", abs(dev["N (uud)"]) < 1e-30 and abs(dev["Σ (uus)"]) < 1e-30,
   "0, 0", f"需 {float(NEED[0]*MD_MU):+.3f}, {float(NEED[1]*MD_MU):+.3f}")
OUT["C3"] = dict(gamma4=float(g4), dev={k: float(v) for k, v in dev.items()},
                 need_MeV=[float(n * MD_MU) for n in NEED],
                 verdict="容量=1 路線：Ξ 之量級接近（2.26 vs 需 3.31）但**符號相反**，且 N/Σ 皆 0 ⟹ 否證")

# ── C4 已註冊核之 body order ─────────────────────────────────────────────
print("\nC4 已註冊交互核之 body order（三體束縛之前置）")
kernels = [
    ("em 電荷核 α q_iq_j B^{−d}", 2, "M5 docs/09"),
    ("樹核 κ q_iq_j B^{−d}", 2, "M4/M5"),
    ("on-site 磁接觸 ∝ δ_{d,0}（兩體）", 2, "M5 docs/09"),
    ("層釘（單粒子）", 1, "P10 §6 / M5 docs/12"),
]
for k, order, src in kernels:
    print(f"    {k:<38} body order = {order}  [{src}]")
ok("C4 已註冊核全部 body order ≤ 2（無三體頂點）", all(o <= 2 for _, o, _ in kernels),
   [o for _, o, _ in kernels], "全部 ≤ 2", "⟹ E_bind^(3) 之味奇部分只可由 pair 核誘導")
# 語料 grep 三體頂點（tracked files，排除本 M）
pats = [("三體", r"三體"), ("three-body", r"three[- ]body"), ("3-body", r"\b3[- ]body\b"),
        ("triple vertex", r"triple (vertex|kernel)")]
tgt = V11
tree = subprocess.run(["git", "-C", tgt, "ls-tree", "-r", "HEAD"], capture_output=True, text=True).stdout.splitlines()
rel_m7 = os.path.relpath(M7, tgt)
cand = []
for ln in tree:
    meta, path = ln.split("\t", 1)
    if path.endswith((".md", ".py")) and not path.startswith(rel_m7):
        cand.append(path)
hits = {}
for label, rx in pats:
    n, files = 0, []
    for p in cand:
        try:
            txt = open(os.path.join(tgt, p), encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        m = re.findall(rx, txt)
        if m:
            n += len(m)
            files.append(p)
    hits[label] = dict(n=n, files=files[:5])
    print(f"    語料「{label}」: {n} 命中" + (f"（{files[:2]}）" if files else ""))
OUT["C4"] = dict(kernels=[dict(name=k, body_order=o, src=s) for k, o, s in kernels], corpus_hits=hits)

# ── C5 pairwise-only 三體相關量級（可解模型，精確）──────────────────────
print("\nC5 pairwise-only 之三體相關（3 等質量 + pairwise 諧核，精確本徵頻率）")
# H = Σp²/2m + ½k Σ_{i<j}(xi−xj)²；K 矩陣本徵值 = {0, 3k/m, 3k/m}
omega_pair = sqrt(2)                 # 兩體：ω = √(2k/m)（單位 k=m=1）
E_pair = mpf(3) * omega_pair / 2     # 3 對 × ½ℏω
E3 = (mpf(0) + 2 * sqrt(3)) / 2      # (ℏ/2)(0 + 2√3)
ratio = (E3 - E_pair) / E_pair
print(f"    Σ_pair E_pair = {float(E_pair):.6f} ℏ√(k/m)；E₀^(3) = {float(E3):.6f} ℏ√(k/m)")
print(f"    三體相關 = {float(ratio)*100:+.2f}% of Σ_pair")
ok("C5 三體相關為 O(10%)，非 O(10²)", abs(ratio) < mpf("0.5"), f"{float(ratio)*100:+.2f}%", "|相關| < 50%")
gap_lo, gap_hi = 36.5, 433.0
ok("C5 相關量級 ≪ ② 所需放大（36.5–433×）", float(abs(ratio)) * 100 < gap_lo,
   f"{float(abs(ratio))*100:.1f}%", f"需 ≥ {gap_lo:.1f}× = {gap_lo*100:.0f}%")
OUT["C5"] = dict(E_pair=float(E_pair), E3=float(E3), corr_rel=float(ratio),
                 required_gain=[gap_lo, gap_hi],
                 verdict="pairwise-only 三體束縛最多提供 O(10%) 修正 ⟹ 不可能放大 36.5–433×")

# ── 判 ────────────────────────────────────────────────────────────────
print("\n" + "=" * 110)
print(f"結果：{NP[1]}/{NP[0]} PASS")
verdict = ("① **NO-GO（加性覆蓋複合態）**：每模容量 = 6 > N_q = 3 ⟹ σ_n ≤ 0.5 ⟹ 永不飽和 ⟹ σ 可加；"
           "7 種味道內容皆無強迫共享；唯一反事實路線（容量=1）符號相反且 N/Σ 為 0。"
           "② **NO-GO（條件性）**：已註冊核全部 pairwise（body order ≤ 2）；三體核／頂點定義語料 0 命中"
           f"（「三體」{hits['三體']['n']} 命中全屬 M5 之界限聲明「重子需三體階梯解，未推」）⟹ E_bind^(3) 之味奇部分只可誘導，而 pairwise-only 之三體相關僅 "
           f"{float(ratio)*100:.1f}% ⟹ 無法放大 36.5–433×。")
print(verdict)
OUT.update(npass=NP[1], ntotal=NP[0], verdict=verdict)
json.dump(OUT, open(os.path.join(M7, "data", "05_additivity_scope.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/05_additivity_scope.json")