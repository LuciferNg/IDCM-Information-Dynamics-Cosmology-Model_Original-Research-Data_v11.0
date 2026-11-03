#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/26 — /88 奇模共振（12-F3）＋ 5 個 3 類分割在共振集之可分辨性（15-T2）

────────────────────────────────────────────────────────────────────────
殘留：
  12-F3 「88-鏈奇模共振」—— docs/12 只做了頻率 no-go 與能量 64 活價對；
        **奇模扇區（a 奇）從未分離計算**。
  15-T2 docs/24 §七 T2 —— 5 個互不等價之 3 類分割，在共振三元組上有無分別。

判準（跑前凍結）
  P1 /88 共振三元組：γ^(88)_a ± γ^(88)_b = γ^(88)_m，a,b,m ∈ [1,44] 互異，容差 1e−40。
  P2 奇模 = 三元組內**三元素皆奇**；半奇 = 恰兩奇。
  P3 χ³ 亮 ⟺ a+b+m 奇（docs/16）；χ⁴ 亮 ⟺ a+b+c+d 偶（docs/18）。
  P4 可分辨 = 存在一個不變量，使 3 個 /90 三元組之標籤集合兩兩不相交。
性質：確定性窮舉。N_trials = 0。dps = 60。零觀測量。
────────────────────────────────────────────────────────────────────────
"""
import itertools
import mpmath as mp

mp.mp.dps = 60
L = "=" * 78
N = 44
TOL = mp.mpf('1e-40')


def gam(x, n):
    return 4 * mp.sin(n * mp.pi / x) ** 2


g88 = [mp.mpf(0)] + [gam(88, n) for n in range(1, N + 1)]
g90 = [mp.mpf(0)] + [gam(90, n) for n in range(1, N + 1)]

print(L)
print("M14 scripts/26 — /88 奇模共振（12-F3）＋ 3 類分割可分辨性（15-T2）")
print(L)

# ============================================================================
print("\n【§一】/88 三波共振全枚舉（sum 與 beat，a,b,m 互異）")
res = set()
for a in range(1, N + 1):
    for b in range(1, N + 1):
        if a == b:
            continue
        for m in range(1, N + 1):
            if m in (a, b):
                continue
            if abs(g88[a] + g88[b] - g88[m]) < TOL or abs(g88[a] - g88[b] - g88[m]) < TOL:
                res.add(tuple(sorted((a, b, m))))
res = sorted(res)
print(f"  共振三元組（無序，去重）= **{len(res)}** 個")
print(f"  全部：{res}")
has44 = sum(1 for t in res if 44 in t)
print(f"  含 Nyquist 44 者 = {has44}/{len(res)} ⟹ **{'全部' if has44 == len(res) else '非全部'}**")
print(f"  形態檢查：每組皆為 {{a, 44−a, 44}} ？ "
      f"{all(sorted(t) == sorted((a, 44 - a, 44)) for t in res for a in [min(t)])}")

# ============================================================================
print("\n【§二】parity 分類（P2 修正：44 恆偶 ⟹ 按**價對**判，非按三元素）")
odd_pair, even_pair, all_odd = [], [], []
for t in res:
    a = min(t)
    pair = (a, 44 - a)
    if all(x % 2 == 1 for x in pair):
        odd_pair.append(t)
    else:
        even_pair.append(t)
    if all(x % 2 == 1 for x in t):
        all_odd.append(t)
print(f"  價對皆奇（奇模扇區）= **{len(odd_pair)}** 個：a = {[min(t) for t in odd_pair]}")
print(f"  價對皆偶 = **{len(even_pair)}** 個：a = {[min(t) for t in even_pair]}")
print(f"  三元素皆奇 = {len(all_odd)}（必然 0：44 恆偶 ⟹ 任何三元組至少一偶；此即先前分類之錯處）")
k2 = sum(1 for t in res if sum(1 for x in t if x % 2 == 1) == 2)
k0 = sum(1 for t in res if sum(1 for x in t if x % 2 == 1) == 0)
print(f"  按元素計：恰兩奇 = {k2} 個、全偶 = {k0} 個、全奇 = 0 ⟹ 恰兩奇 ⟺ 價對皆奇")
print(f"  ⟹ 奇模扇區 {len(odd_pair)} 組 = {[min(t) for t in odd_pair]}；11 = L₅（僅記錄，不作因果宣稱）")

# ============================================================================
print("\n【§三】χ³／χ⁴ 判定（P3）")
for t in res:
    s3 = sum(t)
    assert s3 % 2 == 0, t
print(f"  全部 {len(res)} 組之 a+b+m = 88（偶）⟹ χ³ **全暗**（χ³ 亮 ⟺ 和奇）")
print("  ⟹ /88 三波共振全部宇稱禁戒（與 /90 之 3 組同命運）")
print("  χ⁴（四波）則亮：和偶 ⟹ 可自耦合（docs/18 已立）")

# ============================================================================
print("\n【§四】/90 對照")
res90 = set()
for a in range(1, N + 1):
    for b in range(1, N + 1):
        if a == b:
            continue
        for m in range(1, N + 1):
            if m in (a, b):
                continue
            if abs(g90[a] + g90[b] - g90[m]) < TOL:
                res90.add(tuple(sorted((a, b, m))))
res90 = sorted(res90)
print(f"  /90 三波共振 = {len(res90)} 個：{res90}（全偶）")
print(f"  對照：/88 = {len(res)} 個（全含 44）；兩者 χ³ 皆全暗、χ⁴ 皆亮")
print(f"  ⟹ 兩鏈之三波扇區**質上同型**（暗），差異只在數目（{len(res)} vs {len(res90)}）與成分")

# ============================================================================
print("\n【§五】15-T2：5 個 3 類分割在共振三元組上能否分辨")
INV5 = {
    "i4 Δ%3": lambda a, b: (b - a) % 3,
    "i5 Σ%3": lambda a, b: (a + b) % 3,
    "i7 |∩{9,15,27}|": lambda a, b: len({a, b} & {9, 15, 27}),
    "i8 |∩{18,30,36}|": lambda a, b: len({a, b} & {18, 30, 36}),
    "i9 T_pipe": lambda a, b: (a <= 27, b <= 27),
}
print("  對象：/90 之 3 個三元組；每個三元組取 3 個價對之標籤集合")
sets = {}
for t in res90:
    pairs = list(itertools.combinations(t, 2))
    sets[t] = {k: {f(*p) for p in pairs} for k, f in INV5.items()}
for t in res90:
    print(f"    {t}：")
    for k in INV5:
        print(f"        {k:<18} 標籤集 = {sets[t][k]}")
sep = {}
for k in INV5:
    A, B, C = (sets[t][k] for t in res90)
    if A.isdisjoint(B) and B.isdisjoint(C) and A.isdisjoint(C):
        sep[k] = "三分離"
    elif A.isdisjoint(B) or A.isdisjoint(C) or B.isdisjoint(C):
        sep[k] = "部分分離"
    else:
        sep[k] = "不可分"
print("\n  可分辨性判定：")
for k, v in sep.items():
    print(f"    {k:<18} → {v}")
good = [k for k, v in sep.items() if v == "三分離"]
print(f"\n  能三分離之不變量 = {good if good else '**無**'}")
print("  ⟹ 15-T2 得答：" + ("5 個分割**皆不能**在共振三元組上三分離 "
                            "⟹ 5 者對共振結構完全等價（無優選）"
                            if not good else f"僅 {good} 可三分離，但其屬指標算術（docs/24 §三 canonical=無）"))

# ============================================================================
print("\n【§六】奇模扇區之進一步結構（12-F3 附帶）")
print(f"  奇模 11 組：a = {[a for a in range(1, 44) if a % 2 == 1 and a < 44 - a]}")
print(f"  每組之 γ 和 = 4 = γ_44（Nyquist）⟹ 奇模扇區即『和 = Nyquist』之子族")
print(f"  γ^(88)_44 精確 = {mp.nstr(g88[44], 12)}（= 4，因 sin(π/2) = 1）")
dev = max(abs(g88[a] + g88[44 - a] - 4) for a in range(1, 44))
print(f"  max|γ_a + γ_(44−a) − 4| = {mp.nstr(dev, 3)}")

print("\n" + L)
print("【§七】判定")
print(L)
print(f"""  ① **12-F3 得答**：/88 三波共振恰 **{len(res)}** 組，形態一律 {{a, 44−a, 44}}（全部含 Nyquist 44）；
     按**價對**分扇區：奇模（a,44−a 皆奇）= **{len(odd_pair)}** = L₅（a = 1,3,…,21）；
     偶模 = **{len(even_pair)}**（a = 2,4,…,20）；「三元素皆奇」= 0（44 恆偶，分類須按價對）。
  ② χ³ **全暗**（a+b+m = 88 偶）；χ⁴ 全亮 ⟹ /88 之三波扇區與 /90（3 組，全偶，全暗）**質上同型**。
  ③ **15-T2 得答（否定式）**：5 個互不等價之 3 類分割，在 /90 之 3 個共振三元組上
     {'皆不可三分離' if not good else '部分可分離'} ⟹ 冇任何一個可作「三世代」判別器
     ⟹ 與 docs/24 §三（canonical = 無）一致：即使只睇共振集，5 者亦無優選。
  ④ 奇模 11 = L₅ 之三重身份（/88 折疊 odd core 11 ＝ 本扇區 11 ＝ Fibonacci 11）——
     **僅記錄，不作因果宣稱**（docs/21 §五 紀律）。
  ⑤ 自我更正：初版把「奇模扇區」定義為「三元素皆奇」⟹ 恆 0（44 偶）；
     已改為按價對判。此類錯會令結果顯得「奇模扇區不存在」，屬**假陰性**方向。
""")
print(L)
raise SystemExit(0)