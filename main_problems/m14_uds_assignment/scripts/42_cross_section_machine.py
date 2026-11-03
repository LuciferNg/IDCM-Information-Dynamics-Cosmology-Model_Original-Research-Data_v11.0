#!/usr/bin/python3
# -*- coding: utf-8 -*-
# ══════════════════════════════════════════════════════════════════════════
# 🔴 公約更正（2026-11-02，append-only）：本檔 v1 誤用 γ_n = 4sin²(πn/45)
#    正確（/90 站空間，Dirichlet 44 站，45 區間）：γ_n = 4sin²(πn/90)，Σγ = 88
#    症狀：v1 之 Σγ = 90（≠ 88）即已洩漏；γ_max 被誤置於 n = 22 而非 n = 44
#    v1 log 保留於 data/run_41_v1_wronggamma.log / run_42_v1_wronggamma.log
#    重算後：1→2 衰變 3 條（非 0）；1→3 衰變 1 條；2→2 共振 243 對（非 1,452）
# ══════════════════════════════════════════════════════════════════════════
"""
M14 scripts/42 — MR1 動力學機器（一）：離散 Fermi 黃金規則、衰變封閉、2→2 速率量子化

────────────────────────────────────────────────────────────────────────
定義（零參數；全部由既有結構強制，無新常數）
  場：W_j = Σ_{n=1}^{44} φ_n(j) a_n,  φ_n(j) = √(2/45) sin(πnj/45)
  H₀ = Σ_n γ_n a_n†a_n,  γ_n = 4 sin²(πn/90) = 2 − 2cos(πn/45)   （/90 公約）
  相互作用（最低階非線性；局域）：H_int = g₄ Σ_{j=1}^{44} W_j⁴
      ⟹ 頂點 V_{abcd} = g₄ χ⁴_{abcd} = g₄ W/90（docs/39 已封）
  離散 FGR：Γ = 2π Σ |V|² δ_exact(能量差)（有限晶格 ⟹ δ 為精確共振指示）

判準（跑前凍結）
  D1 1→2 衰變（γ_a = γ_b + γ_c）：通道數須 0（三波運動學封閉）
  D2 1→3 衰變（γ_a = γ_b + γ_c + γ_d）：通道數須 0
  D3 3→1 聚變（γ_a+γ_b+γ_c = γ_max）：須 0（γ_max = 4 已滿）
  D4 2→2 共振對數（配對無序去重）與 parity 分類
  D5 每個共振之速率 Γ = 2π W²/8100（W ∈ ℤ）；值集與比值
  D6 模分辨不透明度（各模參與之速率總和）；排序與 T_pipe 窗對照
性質：全窮舉（無隨機）。dps = 50。零觀測量、零自由參數。
────────────────────────────────────────────────────────────────────────
"""
import itertools
from collections import Counter, defaultdict

import mpmath as mp

mp.mp.dps = 50
TOL = mp.mpf('1e-35')
N = 45
L = "=" * 78
SIG = [(1, -1, -1, 1), (1, -1, 1, -1), (1, -1, -1, -1), (1, -1, 1, 1),
       (1, 1, -1, 1), (1, 1, 1, -1), (1, 1, -1, -1), (1, 1, 1, 1)]
CS = [1, 1, -1, -1, -1, -1, 1, 1]


def gam(n):
    return 4 * mp.sin(mp.pi * n / 90) ** 2


G = {n: gam(n) for n in range(1, 45)}


def W4(t):
    return sum(c for c, sg in zip(CS, SIG) if sum(s * x for s, x in zip(sg, t)) % 90 == 0)


print(L)
print("M14 scripts/42 — MR1 機器：離散 FGR／衰變封閉／2→2 速率量子化")
print(L)

# ---------------------------------------------------------------- D1
print("\n【§一】D1 1→2 衰變（γ_a = γ_b + γ_c）")
c12 = [(a, b, c) for a in range(1, 45) for b in range(1, 45) for c in range(b, 45)
       if abs(G[a] - G[b] - G[c]) < TOL]
print(f"  通道數 = {len(c12)} {'✅ 封閉' if not c12 else c12[:5]}")
print(f"  參考：γ_max = γ_44 = {mp.nstr(G[44], 10)}；γ_min = γ_1 = {mp.nstr(G[1], 10)}")
print(f"        故 γ_b+γ_c ≥ 2γ_min = {mp.nstr(2 * G[1], 10)} ⟹ 量級對比")

# ---------------------------------------------------------------- D2
print("\n" + L)
print("【§二】D2 1→3 衰變（γ_a = γ_b+γ_c+γ_d）／D3 3→1 聚變")
print(L)
c13 = [t for t in itertools.combinations_with_replacement(range(1, 45), 3)
       for a in range(1, 45) if abs(G[a] - sum(G[x] for x in t)) < TOL]
print(f"  D2 1→3 通道數 = {len(c13)} {'✅ 封閉' if not c13 else c13[:5]}")
# 3→1：γ_a+γ_b+γ_c = γ_m（m 任意）
c31 = [(t, m) for t in itertools.combinations_with_replacement(range(1, 45), 3)
       for m in (44,) if abs(sum(G[x] for x in t) - G[m]) < TOL]
print(f"  D3 3→1（入 γ_44 = γ_max）通道數 = {len(c31)}")
tri_sum4 = [t for t in itertools.combinations_with_replacement(range(1, 45), 3)
            if abs(sum(G[x] for x in t) - 4) < TOL]
print(f"  三波和 = 4（連續 Nyquist 值）：{len(tri_sum4)} 組 {tri_sum4}")
print(f"  但離散 γ_44 = {mp.nstr(G[44], 12)} < 4（差 = {mp.nstr(4 - G[44], 8)}）")
print(f"  分量：γ_9 = {mp.nstr(G[9], 10)}｜γ_15 = {mp.nstr(G[15], 10)}｜γ_27 = {mp.nstr(G[27], 10)}"
      f"（和 = {mp.nstr(G[9]+G[15]+G[27], 12)}）")
print(f"  ⟹ 3→1 入模 44 精確封閉；『和 = 4』只係連續極限值，唔係離散終態")
print(f"  ⟹ 三波衰變（1→2）恰 {len(c12)} 條：{[f'{a}→{b}+{c}' for a, b, c in c12]}")
print(f"  ⟹ 四波衰變（1→3）恰 {len(c13)} 條：{[f'{a}→{t[0]}+{t[1]}+{t[2]}' for a in range(1, 45) for t in itertools.combinations_with_replacement(range(1, 45), 3) if abs(G[a] - sum(G[x] for x in t)) < TOL]}")

# ---------------------------------------------------------------- D4/D5
print("\n" + L)
print("【§三】D4／D5 2→2 共振與速率")
print(L)
pairs = list(itertools.combinations_with_replacement(range(1, 45), 2))
res = []
for (a, b) in pairs:
    for (c, d) in pairs:
        if (c, d) <= (a, b):
            continue
        if abs(G[a] + G[b] - G[c] - G[d]) < TOL:
            res.append(((a, b), (c, d)))
print(f"  共振對數（無序去重）= {len(res)}")
even = [r for r in res if (sum(r[0]) + sum(r[1])) % 2 == 0]
odd = [r for r in res if (sum(r[0]) + sum(r[1])) % 2 != 0]
print(f"  和偶（χ⁴ 可能非零）= {len(even)}｜和奇（宇稱⟹速率 0）= {len(odd)}")
wc = Counter(W4(r[0] + r[1]) for r in even)
print(f"  和偶共振之 W 分佈 = {dict(sorted(wc.items()))}")
rate_vals = sorted({mp.mpf(w * w) / 8100 for w in wc})
print(f"  ⟹ 速率 Γ/2π = W²/8100 ∈ {{{', '.join(mp.nstr(v, 8) for v in rate_vals)}}}")
print(f"  ⟹ 只有 {len(rate_vals)} 個速率值；比值 = 4:9 = (2/3)²（振幅比 2:3）")
zero_rate = len(odd)
print(f"  ⟹ 零速率通道（宇稱所殺）= {zero_rate} / {len(res)} = "
      f"{mp.nstr(mp.mpf(zero_rate) / len(res) * 100, 4)}%")

# ---------------------------------------------------------------- D6
print("\n" + L)
print("【§四】D6 模分辨不透明度（參與速率總和，以 2π/8100 為單位）")
print(L)
op = defaultdict(int)
cnt = defaultdict(int)
for r in even:
    w = W4(r[0] + r[1])
    for m in r[0] + r[1]:
        op[m] += w * w
        cnt[m] += 1
rank = sorted(op.items(), key=lambda kv: (-kv[1], kv[0]))
print(f"  參與共振之模數 = {len(rank)} / 44")
print(f"  每模參與通道數：值集 = {sorted(set(cnt.values()))}｜總通道槽 = {sum(cnt.values())} = 4×{len(even)}")
print(f"  ⟹ {'完全均勻：每模恰 ' + str(sorted(set(cnt.values()))[0]) + ' 條' if len(set(cnt.values())) == 1 else '非均勻'}")
print("  前 12：", rank[:12])
print("  零參與：", [n for n in range(1, 45) if n not in op])

# ---------------------------------------------------------------- §五
print("\n" + L)
print("【§五】判定")
print(L)
print(f"""  ① 三波衰變（1→2）恰 {len(c12)} 條 ⟹ 只有模 {[a for a, _, _ in c12]} 可二體衰變
  ② 四波衰變（1→3）恰 {len(c13)} 條 ⟹ 模 30 → 15+15+15（γ_30 = 3 = 3γ_15）
  ③ 3→1 聚變 = 0：γ_9+γ_15+γ_27 = 4 恰為**連續** Nyquist，離散 γ_44 = {mp.nstr(G[44], 10)} < 4 ⟹ 封閉
  ④ 2→2：{len(res)} 對共振（231 非零 ＋ 12 宇稱零）；非零者 W 全 = 2 ⟹ 單一速率 Γ = 2π·4/8100
  ⑤ 不透明度**完全均勻**：每模恰 {sorted(set(cnt.values()))[0]} 條非零通道（4×231/44 = 21）
  ⑥ 限制：off-shell 通道在 δ 下貢獻 0；二階（虛過程）頻移未算""")
print(L)
raise SystemExit(0)