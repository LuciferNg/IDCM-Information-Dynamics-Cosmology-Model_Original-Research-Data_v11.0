#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/15 — 通用共振定律（N-無關）：γ_n + γ_(N−n) = 4 逼出全部在殼價對

────────────────────────────────────────────────────────────────────────
主張（待檢）：
  鏈長 N，γ_n = 4sin²(nπ/(2N))，n = 1..N（44-鏈 N=44；88-鏈 N=88）。
  價對 a<b：
    sum 通道在殼  γ_a + γ_b = γ_m  ⟺  a + b = N （落 m = N，γ = 4 = Nyquist）
    beat 通道在殼 γ_b − γ_a = γ_m  ⟺  b = N     （落 m = N − a）
  證明：γ_n + γ_(N−n) = 4。
  (i)  a+b = N ⟹ γ_a+γ_b = γ_a+γ_(N−a) = 4 = γ_N
  (ii) b = N   ⟹ γ_N−γ_a = 4−γ_a = γ_(N−a)
  本腳本：(a) 逐值核恆等式；(b) 兩條鏈**掃盡** (a,b,m) 之精確共振；
          (c) 檢查有無「不符定律」之額外共振。

性質：確定性掃盡，N_trials = 0。零觀測量。dps = 40，閾 1e−30。
⚠ 陷阱（已實測）：float64 預篩用 1e−9 會漏 (21,67,N=88) —— 必須高精度全掃。
────────────────────────────────────────────────────────────────────────
"""
import mpmath as mp

mp.mp.dps = 40
L = "=" * 78
TOL = mp.mpf('1e-30')


def run(N):
    U = mp.pi / (2 * N)
    g = [mp.mpf(0)] * (N + 1)
    for n in range(1, N + 1):
        g[n] = 4 * mp.sin(n * U) ** 2
    print(f"\n{'-'*78}\n【N = {N}】鏈長 {N}，γ_n = 4sin²(nπ/{2*N})，Nyquist γ_N = {mp.nstr(g[N], 8)}")
    # (a) 恆等式
    mx = max(abs(g[n] + g[N - n] - 4) for n in range(1, N))
    print(f"  γ_n + γ_(N−n) = 4 ：max|dev| = {mp.nstr(mx, 3)}（n=1..{N-1}）")
    # (b) 掃盡
    exS, exD = [], []
    pairs = [(a, b) for a in range(1, N + 1) for b in range(a + 1, N + 1)]
    for (a, b) in pairs:
        s = g[a] + g[b]
        d = g[b] - g[a]
        for m in range(1, N + 1):
            if abs(s - g[m]) < TOL:
                exS.append((a, b, m))
            if abs(d - g[m]) < TOL:
                exD.append((a, b, m))
    setS = {(a, b) for (a, b, _) in exS}
    setD = {(a, b) for (a, b, _) in exD}
    print(f"  價對 = {len(pairs)}；三格 = {len(pairs)*N}")
    print(f"  sum  精確解 = {len(exS)}（價對 {len(setS)}）  命中 m 集 = {sorted({m for *_, m in exS})}")
    print(f"  beat 精確解 = {len(exD)}（價對 {len(setD)}）  命中 m 數 = {len({m for *_, m in exD})}")
    # (c) 對照定律預測
    predS = {(a, b) for (a, b) in pairs if a + b == N}
    predD = {(a, b) for (a, b) in pairs if b == N}
    print(f"  定律預測：sum ⟺ a+b=N ⟹ {len(predS)} 對；beat ⟺ b=N ⟹ {len(predD)} 對")
    print(f"  相符？sum {setS == predS} ；beat {setD == predD}")
    print(f"  互斥？交集 = {len(setS & setD)}；活價對 = {len(setS | setD)} / {len(pairs)}"
          f"（兩通道皆 off-shell = {len(pairs) - len(setS | setD)}）")
    return len(predS), len(predD), len(pairs)


print(L)
print("M14 scripts/15 — 通用共振定律（N-無關）")
print(L)
r44 = run(44)
r88 = run(88)

print("\n" + L)
print("【§三】定律與推論")
print(L)
print(f"""  定律（精確，來自 γ_n+γ_(N−n)=4）：
    sum 在殼 ⟺ a + b = N ；beat 在殼 ⟺ b = N ；兩者互斥（a ≥ 1 ⟹ b=N 與 a+b=N 不可同真）。
    N = 44 ： 21（sum）＋ 43（beat）＝ 64 活價對 / 946
    N = 88 ： {r88[0]}（sum）＋ {r88[1]}（beat）＝ {r88[0]+r88[1]} 活價對 / {r88[2]}
  ⟹ 定律與 N 無關（只依賴同一條恆等式）⟹ 唔係 44 鏈之巧合。
  推論：
    (1) 44-鏈之活價對 = 88-鏈活價對之「偶×偶」子集（even a,b ⟹ 44-鏈模）。
    (2) 每條鏈之 Nyquist 模 n=N 係 beat 共振之**必要**成分（b=N 唯一）。
    (3) sum 共振落於 γ=4（Nyquist 能量），beat 共振落於 γ_(N−a)（鏡像能量）。
    (4) 活價對佔比 = {r44[0]+r44[1]}/{r44[2]} = {(r44[0]+r44[1])/r44[2]:.4f}（N=44）；
        N=88 為 {(r88[0]+r88[1])/r88[2]:.4f} ⟹ 隨 N 趨近 2/N（稀疏）。
""")
print(L)
raise SystemExit(0)