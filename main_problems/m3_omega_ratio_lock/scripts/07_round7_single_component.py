#!/usr/bin/env python3
"""M3 scripts/07 — 第七輪：單分量（為何 W ∈ ℂ 而非 ℂ^k）

命題: 序參量分量數 k = dim ker L_phase（相位扇區之零模重數）
      相位只來自 K_chain（O204: K_chain = 形狀/相位；K_tree = 標度，無相位）
      ⟹ k = dim ker L(C₈₈) = 1 ⟹ W ∈ ℂ ⟹ 單一相位 ⟹ 1D 鏈（第六輪）

檢定（全精確有理算術，Fraction Gaussian elimination）:
  §1 C₈₈ Laplacian: rank = 87, nullity = 1
  §2 連通路徑 C₈₈−e（仍連通）: nullity = 1 ⟹ 驅動因素係【連通性】, 唔係 1D 性
  §3 斷開 2×C₈₈: nullity = 2 ⟹ k = 連通分量數（定理）
  §4 凝聚子空間 = ker(T−I) ∩ ker(R−I) 之維 = 1（唯一 n=0）
  §5 誠實負面: 預算比例係 projective ⟹ (3:16:41) 不能判 k
"""
from fractions import Fraction as F

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

def rank_exact(M):
    """Fraction Gaussian elimination; returns rank."""
    M = [row[:] for row in M]
    rows, cols = len(M), len(M[0])
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i][c] != 0: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        r += 1
        if r == rows: break
    return r

def cycle_lap(N, edges=None):
    if edges is None:
        edges = [(i, (i+1) % N) for i in range(N)]
    M = [[F(0)]*N for _ in range(N)]
    for a, b in edges:
        M[a][a] += 1; M[b][b] += 1
        M[a][b] -= 1; M[b][a] -= 1
    return M

N = 88

# ── §1 精確零模重數 ────────────────────────────────────────────────
print("§1 C₈₈ Laplacian 精確零模重數（ℚ 上）")
L = cycle_lap(N)
rk = rank_exact(L)
print(f"  rank = {rk}, nullity = {N - rk}")
chk("rank(L(C₈₈)) = 87 ⟹ nullity = 1", rk == 87 and N-rk == 1, f"rank={rk}")

# ── §2 連通路徑（非 1D 環，仍連通） ───────────────────────────────
print("\n§2 控制組：C₈₈ 去一條邊（變路徑 P₈₈，仍連通）")
path_edges = [(i, i+1) for i in range(N-1)]
Lp = cycle_lap(N, path_edges)
rkp = rank_exact(Lp)
print(f"  rank = {rkp}, nullity = {N - rkp}")
chk("路徑 P₈₈ nullity 仍 = 1 ⟹ 驅動因素 = 連通性（非 1D 性）", N-rkp == 1, f"{N-rkp}")

# ── §3 斷開：兩個分量 ─────────────────────────────────────────────
print("\n§3 斷開 2×C₈₈（兩個不相交世界）")
M2 = [[F(0)]*(2*N) for _ in range(2*N)]
for blk in (0, 1):
    c = cycle_lap(N)
    for i in range(N):
        for j in range(N):
            M2[blk*N+i][blk*N+j] = c[i][j]
rk2 = rank_exact(M2)
print(f"  rank = {rk2}, nullity = {2*N - rk2}")
chk("2 分量 ⟹ nullity = 2 ⟹ k = 連通分量數", 2*N-rk2 == 2, f"{2*N-rk2}")

# ── §4 凝聚子空間：平移不變 ∩ 反射不變（精確矩陣） ────────────────
print("\n§4 凝聚子空間維數（精確：ker(T−I) ∩ ker(R−I)）")
# T = shift, R = reflection n ↦ −n mod N
Tm = [[F(0)]*N for _ in range(N)]
for i in range(N):
    Tm[(i+1) % N][i] = F(1)              # (T v)[n] = v[n−1]
Rm = [[F(0)]*N for _ in range(N)]
for i in range(N):
    Rm[(-i) % N][i] = F(1)               # (R v)[n] = v[−n]

def nullity(M):
    return len(M[0]) - rank_exact(M)

nT = nullity([[Tm[i][j] - (1 if i == j else 0) for j in range(N)] for i in range(N)])
nR = nullity([[Rm[i][j] - (1 if i == j else 0) for j in range(N)] for i in range(N)])
stacked = [[Tm[i][j] - (1 if i == j else 0) for j in range(N)] for i in range(N)] + \
          [[Rm[i][j] - (1 if i == j else 0) for j in range(N)] for i in range(N)]
nTR = nullity(stacked)
print(f"  nullity(T−I) = {nT} ; nullity(R−I) = {nR} ; nullity(交集) = {nTR}")
chk("ker(T−I) 維 = 1（平移不變空間僅均勻模）", nT == 1, f"{nT}")
chk("ker(R−I) 維 = 45（n=0,22,44 + 21 對）", nR == 45, f"{nR}")
chk("交集維 = 1 ⟹ 凝聚方向唯一", nTR == 1, f"{nTR}")
chk("nR = N/2 + 1", nR == N//2 + 1, f"{nR}")

# ── §5 真正判據：b₁（獨立相位數 = 環秩） ──────────────────────────
print("\n§5 ★ 獨立相位數 = b₁（相位係 holonomy，需環）")
def edges_cycle(N):  return [(i, (i+1) % N) for i in range(N)]
def edges_path(N):   return [(i, i+1) for i in range(N-1)]
def edges_2cyc(N):   return [(b*N+i, b*N+(i+1) % N) for b in (0, 1) for i in range(N)]
def edges_torus(m, n):
    E = []
    for a in range(m):
        for b in range(n):
            E.append((a*n+b, ((a+1) % m)*n+b))
            E.append((a*n+b, a*n+(b+1) % n))
    return E

def components(V, E):
    par = list(range(V))
    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in E:
        ra, rb = find(a), find(b)
        if ra != rb: par[ra] = rb
    return len({find(i) for i in range(V)})

def b1(V, E):
    return len(E) - V + components(V, E)

graphs = [
    ("C₈₈（單環，世界鏈）", 88, edges_cycle(88)),
    ("P₈₈（路徑，去一邊）", 88, edges_path(88)),
    ("2×C₈₈（兩斷開世界）", 176, edges_2cyc(88)),
    ("C₈×C₁₁（2D 環面）", 88, edges_torus(8, 11)),
    ("C₂×C₄₄（2D 環面）", 88, edges_torus(2, 44)),
]
res = {}
for name, V, E in graphs:
    b0 = components(V, E); b = b1(V, E)
    res[name.split("（")[0]] = (V, len(E), b0, b)
    print(f"  {name:24s} V={V:3d} E={len(E):3d} b₀={b0} ⟹ b₁ = {b}")
chk("C₈₈: b₁ = 1 ⟹ 恰好一個獨立相位 ⟹ k=1", res["C₈₈"][3] == 1)
chk("P₈₈: b₁ = 0 ⟹ 冇環 ⟹ 冇相位（鏈必須係環）", res["P₈₈"][3] == 0)
chk("C₈×C₁₁: b₁ = 89 ⟹ 相位數遠多於 1 ⟹ 排除", res["C₈×C₁₁"][3] == 89)
chk("C₂×C₄₄: b₁ = 89 ⟹ 同上", res["C₂×C₄₄"][3] == 89)
chk("2×C₈₈: b₀=2 且 b₁=2", res["2×C₈₈"][2] == 2 and res["2×C₈₈"][3] == 2)
print("  ⟹ k = b₁(相位扇區)；C₈₈ 之 b₁=1 係【單環】拓撲之後果")

# ── §6 誠實負面：預算 projective ──────────────────────────────────
print("\n§6 誠實負面：預算能否判 k？")
for k in (1, 2, 3):
    tot = 60*k
    ratio = (3*k, 16*k, 41*k)
    norm = tuple(F(x, tot) for x in ratio)
    print(f"  k={k}: 分子×{k} = {ratio}, 除以 {tot} ⟹ {norm[0]}, {norm[1]}, {norm[2]}")
chk("比例對 k 完全不變（projective）⟹ 預算【不能】判 k",
    len({(F(3*k,60*k), F(16*k,60*k), F(41*k,60*k)) for k in (1,2,3)}) == 1,
    "此路線無效，記錄為 FAIL-ROUTE")

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")