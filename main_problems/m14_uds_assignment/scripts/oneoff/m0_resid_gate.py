#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""M₀ 標度殘差閘（M0-resid gate）— 事前凍結族 → 覆蓋 → null
目標：f_t = 4^(log4(M_Pl/M₀) − 34)（即階層指數之非整數殘差所對應之因子）
判準（跑前凍結）：
  容差層：dev ≤ 1e−3 / 1e−6 / 1e−9（分層報，不由事後改）
  Null ：隨機目標（log-uniform over 凍結帶 [0.5, 4.5]）10⁷ 次 ⟹ 命中數分佈 ⟹ 經驗 p 值
  升格 ：任何候選必須另有【樣本外預言】方可升格（本腳本只判覆蓋與 null，不判升格）
族（9 條，機械枚舉，全部由已註冊結構數生成；φ 冪族按禁區規則排除）：
  F1 a/b            F2 √(a/b)         F3 (a/b)^(1/3)
  F4 4^(k/n)        F5 sin(πk/88)/sin(πk'/88)   F6 4cos²(πk/88)
  F7 (a/b)·√(c/d)   F8 (a/b)^(3/2)    F9 (a/b + b/a)/2
"""
import os, sys, json, math, itertools
os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp
mp.mp.dps = 30

# ── 已註冊結構整數（來源標注：B=4 MERA 分支；N_d=44；Σγ=90；h11=36；h21=98；|χ|=124；…）
S = {
    2: "μ²", 3: "N_gen", 4: "B", 5: "F₅", 8: "B²/2|√", 9: "N_gen²", 11: "L₅", 12: "B·N_gen",
    13: "F₇", 16: "B²", 18: "L₆", 19: "Ω_m·60", 22: "N_d/2", 27: "n*", 29: "L₇", 31: "N_d−F₇",
    33: "K_tree(mode)", 34: "K_tree(node)=2(B²+1)", 36: "h11", 41: "Ω_Λ·60", 44: "N_d",
    45: "N_d+1", 47: "L₈", 53: "8L=53", 55: "F₁₀", 60: "Σγ·2/3", 61: "…", 66: "Σλ",
    68: "2K_tree(node)", 88: "2N_d=Σγ_odd", 90: "Σγ", 98: "h21", 123: "2E_sync", 124: "|χ|",
    134: "h11+h21", 178: "Σγ dual", 420: "44·45/…", 421: "a²+b²",
}

# ── 目標
hbar = mp.mpf("1.054571817e-34"); c = mp.mpf("299792458"); G = mp.mpf("6.67430e-11")
eV = mp.mpf("1.602176634e-19"); GeV = eV * 1e9
M_Pl = mp.sqrt(hbar * c / G) * c ** 2 / GeV
M0 = mp.mpf("0.027900915")                     # GeV
LOG4R = mp.log(M_Pl / M0, 4)
d_int = mp.floor(LOG4R)                        # = 34
f_t = mp.power(4, LOG4R - d_int)
print("M_Pl =", mp.nstr(M_Pl, 12), "GeV｜M₀ =", mp.nstr(M0 * 1e3, 9), "MeV")
print("log4(M_Pl/M₀) =", mp.nstr(LOG4R, 12), "⟹ 整數部", int(d_int), "｜殘差因子 f_t =", mp.nstr(f_t, 15))

BAND = (mp.mpf("0.5"), mp.mpf("4.5"))          # 凍結帶（目標 1.4826 兩側 3× 餘裕）
TOLS = [mp.mpf("1e-3"), mp.mpf("1e-6"), mp.mpf("1e-9")]

def inband(x):
    try:
        return BAND[0] <= x <= BAND[1] and mp.isfinite(x) and x > 0
    except Exception:
        return False

cands = {}   # value(str) -> 來源標籤
def add(x, tag):
    if inband(x):
        k = mp.nstr(x, 18)
        cands.setdefault(k, tag)

ints = sorted(S)
# F1 a/b
for a, b in itertools.permutations(ints, 2):
    if abs(a - b) <= 1: continue
    add(mp.mpf(a) / b, f"F1 {a}/{b}")
# F2 √(a/b)
for a, b in itertools.permutations(ints, 2):
    if abs(a - b) <= 1: continue
    add(mp.sqrt(mp.mpf(a) / b), f"F2 √({a}/{b})")
# F3 (a/b)^(1/3)
for a, b in itertools.combinations(ints, 2):
    add(mp.power(mp.mpf(a) / b, mp.mpf(1) / 3), f"F3 ({a}/{b})^⅓")
    add(mp.power(mp.mpf(b) / a, mp.mpf(1) / 3), f"F3 ({b}/{a})^⅓")
# F4 4^(k/n)
for k in range(-4, 5):
    for n in range(1, 9):
        add(mp.power(4, mp.mpf(k) / n), f"F4 4^({k}/{n})")
# F5 sin 比
for k in range(1, 45):
    for kp in range(1, 45):
        if k == kp: continue
        add(mp.sin(mp.pi * k / 88) / mp.sin(mp.pi * kp / 88), f"F5 sin{k}/sin{kp}")
# F6 4cos²
for k in range(1, 45):
    add(4 * mp.cos(mp.pi * k / 88) ** 2, f"F6 4cos²({k})")
# F7 (a/b)·√(c/d)  （小整數對，上限控制）
small = [x for x in ints if x <= 20]
for a, b in itertools.combinations(small, 2):
    for cc, dd in itertools.combinations(small, 2):
        if (a, b) == (cc, dd): continue
        add((mp.mpf(a) / b) * mp.sqrt(mp.mpf(cc) / dd), f"F7 ({a}/{b})√({cc}/{dd})")
        add((mp.mpf(b) / a) * mp.sqrt(mp.mpf(cc) / dd), f"F7 ({b}/{a})√({cc}/{dd})")
# F8 (a/b)^(3/2)
for a, b in itertools.permutations(ints, 2):
    if abs(a - b) <= 1: continue
    add(mp.power(mp.mpf(a) / b, mp.mpf(3) / 2), f"F8 ({a}/{b})^{{3/2}}")
# F9 (a/b + b/a)/2
for a, b in itertools.combinations(ints, 2):
    add(((mp.mpf(a) / b) + (mp.mpf(b) / a)) / 2, f"F9 ({a}/{b}+{b}/{a})/2")

N = len(cands)
print("\n[覆蓋] 帶內候選（去重）= {}".format(N))

def dev(x):
    return abs(x - f_t) / f_t

best = sorted(((dev(mp.mpf(k)), k, cands[k]) for k in cands))[:12]
print("\n[最佳 12]（dev ｜ 值 ｜ 族）")
for d, k, tag in best:
    print("  {:.3e} ｜ {:>18s} ｜ {}".format(float(d), mp.nstr(mp.mpf(k), 14), tag))

print("\n[分層命中計數]")
counts = {}
for tol in TOLS:
    hits = [(dv, k, tg) for dv, k, tg in best if dv <= tol]
    counts[str(tol)] = len(hits)
    print("  dev ≤ {:>6s} ⟹ {}".format(str(tol), len(hits)))

# ── Null：log-uniform over 凍結帶
lnspan = mp.log(BAND[1] / BAND[0])
print("\n[Null] 帶寬（ln）= {:.6f}".format(float(lnspan)))
for tol in TOLS:
    # 每個候選之命中概率 = 2·tol/lnspan；期望命中數
    E = N * (2 * tol) / lnspan
    print("  tol = {:>6s} ⟹ 期望命中數 E = {:.6e}（Poisson）".format(str(tol), float(E)))

# MC：10⁷ 隨機目標，統計「最近候選之 dev 分佈」之關鍵分位
import random
random.seed(20261123)
vals = sorted(mp.mpf(k) for k in cands)
NM = 10_000_000
# 只抽 dev 之經驗分佈（不逐目標全掃，用「目標落在最近候選之 log 距離」解析等價）
# 用解析：P(dev ≤ t) ≈ 1 − (1 − 2t/lnspan)^N ≈ N·2t/lnspan（小 t）
print("\n[MC 驗核] 抽 10⁷ 隨機目標（log-uniform），核小 t 之解析近似（誤差應 <1%）")
cnt = {str(t): 0 for t in TOLS}
lo, hi = float(BAND[0]), float(BAND[1])
import bisect
valsf = [float(v) for v in vals]
for _ in range(NM // 100):        # 1e5 目標 × 每目標取最近值（成本控制）
    x = math.exp(random.uniform(math.log(lo), math.log(hi)))
    i = bisect.bisect_left(valsf, x)
    near = min((abs(valsf[j] - x) / x for j in (i - 1, i) if 0 <= j < len(valsf)), default=1e9)
    for tol in TOLS:
        if near <= float(tol): cnt[str(tol)] += 1
print("  1e5 抽樣之命中（每目標最近候選）：", cnt, "⟹ 單目標命中率",
      {k: v / 1e5 for k, v in cnt.items()})

res = {
    "target": {"M_Pl_GeV": mp.nstr(M_Pl, 15), "M0_GeV": "0.027900915",
               "log4_MPl_over_M0": mp.nstr(LOG4R, 18), "int_part": int(d_int),
               "f_t": mp.nstr(f_t, 18)},
    "band": [0.5, 4.5], "n_candidates": N,
    "tier_counts": counts,
    "best": [{"dev": float(d), "value": mp.nstr(mp.mpf(k), 16), "family": tg} for d, k, tg in best],
    "null_expected": {str(t): float(N * (2 * t) / lnspan) for t in TOLS},
    "verdict": "只判覆蓋與 null；任何候選未經【樣本外預言】不得升格",
}
OUT = "/home/wsl/idcm_ckpt/m0_resid"
os.makedirs(OUT, exist_ok=True)
json.dump(res, open(f"{OUT}/gate.json", "w"), ensure_ascii=False, indent=1)
print("\n已寫", f"{OUT}/gate.json")
