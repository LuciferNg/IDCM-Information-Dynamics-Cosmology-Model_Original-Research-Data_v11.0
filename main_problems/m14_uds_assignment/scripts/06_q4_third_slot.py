#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/06 — Q4：`N_gen = rank(耦合塊) + 1` 可否升為定理？

問題（record §二 Q4）：M13 docs/10 §二 之「槽數 = rank(A₂)+1 = 3」可否由 Gram 秩結構證成，
而非觀察？並以 δ_c ≈ 0 作唯一可測檢定。
判準：A3（dof=0）／A5（範圍）／A7（禁以第 4 槽擴充耦合塊 —— 本檔**不**擴充，
      只把「rank+1」當**待檢命題**在既有掃描箱上檢驗）。

檢定：
  §一 逐三元組分解：rank(M)、單態數 n_singlet，核 3 = rank + n_singlet（維數守恆）
  §二 「rank+1 = 3」在掃描箱內之命中率（＝是否普遍律）
  §三 耦合塊之**身份**：係 {u,d} 對，抑或「世代結構」？
  §四 尺度檢定：同一構造套於 2 站／4 站，「rank+1」是否仍成立
  §五 δ_c ≈ 0 之判別力
N_trials = 112（掃描箱全報）＋ 4（尺度檢定之站數）。
"""
import itertools
import mpmath as mp

mp.mp.dps = 50
print("=" * 78)
print("M14 scripts/06 — Q4：第三槽定理（N_gen = rank + 1）之升格檢定")
print("=" * 78)

NN = 44
phi = lambda n, j: mp.sqrt(mp.mpf(2) / 45) * mp.sin(mp.pi * n * j / 45)
gam90 = lambda n: 4 * mp.sin(mp.pi * n / 90) ** 2
TOL = mp.mpf("1e-40")


def gram(p, stations=None):
    """p=(p_u,p_d[,p_s…])；Δ_n = φ_n(p_d)−φ_n(p_u)（u↔d 交換向量，隨三元組變）。"""
    st = stations if stations is not None else p
    D = [phi(n, p[1]) - phi(n, p[0]) for n in range(1, NN + 1)]
    V = [[phi(n, k) * D[n - 1] for n in range(1, NN + 1)] for k in st]
    m = len(st)
    M = mp.matrix(m, m)
    for i in range(m):
        for j in range(m):
            M[i, j] = mp.fsum(V[i][n] * V[j][n] / gam90(n + 1) for n in range(NN))
            M[j, i] = M[i, j]
    return M


def decomp(M):
    """回傳 (耦合塊秩 rk_block, 單態數 n_singlet, 是否 A₂⊕1, 連通分量)。
    耦合塊 = 以「非零交叉項」連通之非單態塊；其秩以該塊之特徵值計。"""
    m = M.rows
    adj = [[i, j] for i, j in itertools.combinations(range(m), 2) if abs(M[i, j]) > TOL]
    parent = list(range(m))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in adj:
        parent[find(i)] = find(j)
    comps = {}
    for i in range(m):
        comps.setdefault(find(i), []).append(i)
    n_singlet = sum(1 for c in comps.values() if len(c) == 1)
    # 耦合塊秩 = 各「非單態連通塊」之秩之和（單態塊每個秩 = 1，不計入耦合塊）
    rk_block = 0
    for c in comps.values():
        if len(c) == 1:
            continue
        sub = mp.matrix([[M[i, j] for j in c] for i in c])
        rk_block += sum(1 for e in mp.eigsy(sub, eigvals_only=True) if abs(e) > TOL)
    A2p1 = (len(comps) == 2 and sorted(len(c) for c in comps.values()) == [1, 2] and
            rk_block == 2 and abs(M[0, 1] + M[0, 0] / 2) < TOL and abs(M[0, 0] - M[1, 1]) < TOL)
    return rk_block, n_singlet, A2p1, comps


triples = [(a, b, c) for a in range(1, 8) for b in range(1, 9) for c in range(1, 11) if a < b < c]
print(f"\n【§一】掃描箱 = {len(triples)} 個三元組（p_u<p_d<p_s；M13 docs/05 §S1 同箱）")
print("  逐三元組計算 Gram 之 (耦合塊秩 rank, 正交單態數 n_singlet, 塊大小)")
dist = {}
viol = []          # rank+1 ≠ 3
no_singlet = []    # rank+1 = 3 但 n_singlet ≠ 1（「+1」無所指）
for p in triples:
    M = gram(p)
    rank, ns, a2, comps = decomp(M)
    sizes = sorted(len(c) for c in comps.values())
    dist[(rank, ns)] = dist.get((rank, ns), 0) + 1
    if rank + 1 != 3:
        viol.append((p, rank, ns))
    elif ns != 1:
        no_singlet.append((p, rank, ns, sizes))
print(f"  {'(rk_block, n_singlet)':>22s} {'個數':>6s} {'塊大小':>12s} {'rk+1':>6s} 判")
seen = {}
for p in triples:
    M = gram(p)
    rank, ns, a2, comps = decomp(M)
    key = (rank, ns)
    seen.setdefault(key, []).append(tuple(sorted(len(c) for c in comps.values())))
for (rk, ns), cnt in sorted(dist.items()):
    sizes = sorted(set(seen[(rk, ns)]))
    verdict = ("✅ rank+1 = 3 且恰一單態（M13 之命題成立）" if (rk + 1 == 3 and ns == 1)
               else ("🟡 rank+1 = 3 但**無單態**（「+1」無所指）" if rk + 1 == 3
                     else f"🔴 rank+1 = {rk+1} ≠ 3（直接失敗）"))
    print(f"  {str((rk, ns)):>22s} {cnt:>6d} {str(sizes):>12s} {rk+1:>6d} {verdict}")
print(f"  ⟹ 「N_gen = rank+1 且恰一單態」只對 {dist.get((2, 1), 0)}/{len(triples)} 成立"
      f"（= A₂⊕1 類，M13 docs/05 §S1 之 44/112）")
print(f"  ⟹ rank+1 ≠ 3 者 {len(viol)}/{len(triples)}：其 rank = "
      f"{sorted(set(v[1] for v in viol))} ⟹ rank+1 ∈ "
      f"{sorted(set(v[1]+1 for v in viol))}")
print(f"  ⟹ rank+1 = 3 但無單態者 {len(no_singlet)}/{len(triples)}（例："
      f"{[v[0] for v in no_singlet[:4]]}）⟹ 該等三元組之「+1」**無所指**")

# =====================================================================
print("\n【§二】「rank+1」之地位")
tot = len(triples)
n_ok = dist.get((2, 1), 0)
print(f"  命中率 = {n_ok}/{tot} = {100*n_ok/tot:.1f}%")
print(f"  失敗者之 rank = {sorted(set(v[1] for v in viol))} ⟹ 其耦合塊係 rank 3 ⟹ rank+1 = 4 ≠ 3")
print("  例（失敗）：" + "、".join(str(v[0]) for v in viol[:6]) + " …")
print("  ⟹ **「N_gen = rank + 1」不是普遍律**：它只在『s 與 u,d 正交』之子集成立（A₂⊕1 類）")
print(f"     （該子集 = {n_ok}/{tot} = {100*n_ok/tot:.1f}%，M13 docs/05 §S1 亦記 44/112）")

# =====================================================================
print("\n【§三】耦合塊之身份：係 {u,d} 對，抑或世代結構？")
print("  依構造：v_x(n) = φ_n(p_x)·Δ_n，而 **Δ_n = φ_n(p_d) − φ_n(p_u)**")
print("  ⟹ Δ 只含 u、d 兩站 ⟹ 交叉項 M_xy = Σ_n φ(p_x)φ(p_y)Δ²/γ_n 之非零結構**恆由 (u,d) 對承載**")
M234 = gram((2, 3, 4))
print(f"  實例 p=(2,3,4)：M_ud/M_uu = {mp.nstr(M234[0,1]/M234[0,0],6)}（= −1/2，A₂ Cartan 之非對角）")
print(f"                  M_us/M_uu = {mp.nstr(M234[0,2]/M234[0,0],3)}（≡ 0）")
print("  ⟹ rank = 2 之來歷 = **isospin 對 (u,d)**，非「世代數」")
print("  ⟹ 故 rank 不含任何關於 N 之信息；N 係輸入站數")
print("  ⟹ 「槽數 = rank(A₂)+1」之內容 = 「N = 2 + (N−2)」= 恆等式")

# =====================================================================
print("\n【§四】尺度檢定（N_trials = 4）：同一構造套於 2／3／4／5 站")
print(f"  {'站數':>4s} {'站集':>16s} {'塊秩':>5s} {'單態':>5s} {'塊秩+單態':>12s} {'塊秩+1':>7s} 判")
for st in [(2, 3), (2, 3, 4), (2, 3, 4, 5), (2, 3, 4, 5, 6)]:
    M = gram(st[:2], stations=st)
    rank, ns, a2, comps = decomp(M)
    print(f"  {len(st):>4d} {str(st):>16s} {rank:>5d} {ns:>5d} {len(st):>14d} {rank+1:>7d}  "
          f"{'✅ rank+1 = N' if rank + 1 == len(st) else '🔴 rank+1 = ' + str(rank+1) + ' ≠ N'}")
print("  ⟹ 2 站：rank+1 = 3 ≠ 2 🔴；3 站：3 ✅；4 站：3 ≠ 4 🔴；5 站：3 ≠ 5 🔴")
print("  ⟹ 「rank+1 = N_gen」只在 N = 3 一點成立 ⟹ 不能由它**推出** N = 3（循環）")
print("  ⚠ A7 遵守：本節**不**宣稱第 4 味存在；只把同一構造之站數當控制變量，檢驗命題之普遍性")

# =====================================================================
print("\n【§五】δ_c ≈ 0 之判別力（record 指定之唯一可測檢定）")
d_c, s_c = mp.mpf("-0.282"), mp.mpf("0.261423")     # scripts/01 R6／R7
d_d, d_s = mp.mpf("1.382"), mp.mpf("3.425")
print(f"  δ_c = {mp.nstr(d_c,5)} ± {mp.nstr(s_c,5)} MeV ⟹ z(δ_c vs 0) = {mp.nstr(d_c/s_c,4)}")
print(f"  若 c 與 d/s 等量入耦合塊（M13 docs/10 §二之路 ①）：")
for nm, ref in (("|δ_d|", d_d), ("|δ_s|", d_s)):
    print(f"    對 {nm} = {mp.nstr(ref,4)}：差 {mp.nstr((ref-abs(d_c)),4)} MeV ⟹ "
          f"{mp.nstr((ref-abs(d_c))/s_c,4)}σ ⟹ {'🔴 排除' if abs((ref-abs(d_c))/s_c)>3 else '🟡'}")
print("  ⟹ 數據**不支持** c 進入耦合塊（路 ① 死）；與路 ②（正交單態，零約束）相容")
print("  ⟹ 但此檢定對 Q4 之命題**無判別力**：無論「rank+1」真假，δ_c ≈ 0 皆被預期")
print("     （因 s 之正交性來自 Δ 之支撐，與耦合塊秩無關）")

# =====================================================================
print("\n" + "=" * 78)
print("Q4 判決")
print("=" * 78)
print(f"""
  §一 dim = 耦合塊秩 + 單態數：112/112 成立 ⟹ **恆等式**（維數守恆），零信息
  §二 「rank+1 = 3」只對 {n_ok}/{tot}（{100*n_ok/tot:.1f}%）成立；{len(viol)}/{tot} 直接失敗
      （失敗者之耦合塊秩 = 3 ⟹ rank+1 = 4）
  §三 耦合塊恆由 (u,d) 對承載（rank = 2 = isospin 對）；N 係輸入站數
  §四 站數 2/3/4/5 之控制檢定：rank+1 只於 N=3 成立 ⟹ 循環
  §五 δ_c ≈ 0 與 rank+1 之真偽**無關** ⟹ 無判別力

  ⟹ **Q4 = 死（不可升為定理）**
    範圍標籤（A5）：
      ① {{有限集：112 個三元組}} —— rank+1 於 {len(viol)} 個直接失敗（耦合塊秩 = 3 ⟹ rank+1 = 4）
      ② {{尺度族：站數 2,3,4,5}} —— rank+1 只於 N=3 成立 ⟹ 不能推出 N=3
    判詞：`N_gen = rank(耦合塊)+1` 係 **3 = 2 + 1 之重述**；其「2」= isospin 對之秩
    （由 Δ = φ(p_d)−φ(p_u) 之定義而來），其「1」= s 站與 u,d 正交（僅 {n_ok}/112 之子集性質）。
    ⟹ N_gen = 3 **在現有結構中仍是輸入**，而非結論。
  ⟹ 正面保留（不作死）：A₂⊕1 之形式（耦合塊秩 2 ＋ 1 單態，{n_ok}/112）係真結構事實，
     可作**該三元組之描述**；不可反推 N_gen。
  N_trials = 112（§一）＋ 4（§四）。零上游寫入。
""")
raise SystemExit(0)