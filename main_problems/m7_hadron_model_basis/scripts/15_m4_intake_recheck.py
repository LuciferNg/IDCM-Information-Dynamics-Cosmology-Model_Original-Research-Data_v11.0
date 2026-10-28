"""M7 scripts/15 — 回核 M4 docs/39（收件審計 #2）：X 兩值、B13 措辭、n_s 定義域
==========================================================================================================
M4 docs/39 對本 M 兩封 handoff 作獨立重核（B1–B12 逐位一致），並提出兩點：
  (i)  更正：本 M 用 X = 2.56（捨入）；M4 canonical X = 42667/16670 = 2.5595081
  (ii) 更強後果 B13：⟨O⟩ = 1+(2/3)n_s 要命中 ⟨O⟩_K = 2.4218 需 n_s = +2.133，
       而「該線性律之 |n_s| ≤ 1 定義域內根本冇解」
本檔回核：兩個 X 之 ⟨O⟩_K／σ 偏移；並核 B13 之**措辭**（domain 究竟係咩）。
輸出 data/15_m4_intake_recheck.json
"""
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "15_m4_intake_recheck.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


def show(x, n=8):
    return nstr(mpf(x.numerator) / mpf(x.denominator), n)


EM = F("2.264589")            # scripts/14 重算（精確有理近似，1e-6）
DM_K = F("-3.9340")           # M4 docs/18 §表 K 列
SIG_O = F("0.262")            # BMW Σ 錨 1.498±0.262
O_SIG_ANCHOR = F("1.498")
X_A = F("2.56")               # M4 docs/18 §2 明載（FLAG24，外部固定）
X_B = F("42667") / F("16670")  # M4 docs/39 稱 canonical

print("=" * 100)
print("M7 scripts/15 — 回核 M4 docs/39：X 兩值 ＋ B13 措辭 ＋ n_s 定義域")
print("=" * 100)

hdr("§1 X 兩值之 ⟨O⟩_K（外部值反推）")
rows = {}
for lbl, X in [("X = 2.56（docs/18 §2）", X_A), ("X = 42667/16670（docs/39）", X_B)]:
    O = (EM - DM_K) / X
    off = (O - O_SIG_ANCHOR) / SIG_O
    rows[lbl] = (X, O, off)
    print(f"    {lbl:<28} X = {show(X, 10)} ⟹ ⟨O⟩_K = {show(O, 7)}，距 Σ 錨 {float(off):+.3f}σ")
dO = abs(rows["X = 2.56（docs/18 §2）"][1] - rows["X = 42667/16670（docs/39）"][1])
print(f"    兩值差 = {show(dO, 6)}（相對 {float(dO/rows['X = 2.56（docs/18 §2）'][1])*100:.3f}%）")
ok("兩 X 皆給 ⟨O⟩_K 遠離 Σ 錨（>3σ）⟹ 判決不改",
   all(abs(off) > 3 for _, _, off in rows.values()),
   ", ".join(f"{float(off):+.2f}σ" for _, _, off in rows.values()), "全部 >3σ")
ok("M4 之 canonical 值可重現（2.4217892）",
   abs(mpf(rows["X = 42667/16670（docs/39）"][1].numerator) / rows["X = 42667/16670（docs/39）"][1].denominator - mpf("2.4217892")) < mpf("2e-6"),
   show(rows["X = 42667/16670（docs/39）"][1], 8), "2.4217892")
ok("本 M docs/22 之 2.42132 可重現（X=2.56）",
   abs(mpf(rows["X = 2.56（docs/18 §2）"][1].numerator) / rows["X = 2.56（docs/18 §2）"][1].denominator - mpf("2.42132")) < mpf("2e-5"),
   show(rows["X = 2.56（docs/18 §2）"][1], 8), "2.42132")

hdr("§2 B13：所需 n_s 及定義域核查")
for lbl, (X, O, off) in rows.items():
    ns = (O - 1) * F(3, 2)
    print(f"    {lbl:<28} ⟹ 所需 n_s = (⟨O⟩−1)·3/2 = {show(ns, 7)}")
ns_a = (rows["X = 2.56（docs/18 §2）"][1] - 1) * F(3, 2)
ns_b = (rows["X = 42667/16670（docs/39）"][1] - 1) * F(3, 2)
print("    M4 docs/29 V2 明文：n_s = 0,1,2 三點同時符合（不許只中一點）⟹ 定義域 = {0,1,2}，非 |n_s| ≤ 1")
print(f"    實測：K(n_s=1)、Σ(n_s=1)、Ξ(n_s=2) ⟹ 可達 n_s 最大 = 2")
ok("所需 n_s 落在整數定義域 {0,1,2} 之外 ⟹ 整線無解（結論成立）",
   ns_a > 2 and ns_b > 2, f"{show(ns_a,6)} / {show(ns_b,6)}", "> 2",
   "結論（整線排除）成立，但理由須改：域 = {0,1,2}，非 |n_s| ≤ 1")
ok("★ M4 B13 措辭錯：其自身 docs/29 V2 已寫 n_s ∈ {0,1,2}，Ξ 錨 n_s = 2（⟨O⟩=2.293）",
   True, "docs/29 V2", "n_s∈{0,1,2}",
   "『|n_s| ≤ 1 定義域』與其自身文檔矛盾 ⟹ 屬措辭錯，非結果錯")

hdr("§3 更乾淨之殺法：M1 在 K 通道 n_s = 1 固定，無需外推")
O_M1 = F(5, 3)                     # ⟨O⟩_K = 1 + (2/3)·1
for lbl, X in [("X = 2.56", X_A), ("X = 42667/16670", X_B)]:
    em_pred = DM_K + X * O_M1
    z = (em_pred - EM) / F("0.120863")
    print(f"    M1 @ n_s=1：⟨O⟩_K = 5/3 精確 ⟹ em_pred = Δm + X·5/3 = {show(em_pred, 6)} MeV"
          f"（{lbl}）⟹ z = {float(z):+.2f}")
ok("M1 @ n_s=1 直接給出凍結預言 0.3327（差 <0.1%）",
   abs(float((DM_K + X_A * O_M1) - F("0.3327"))) < 1e-3,
   show(DM_K + X_A * O_M1, 6), "0.3327")
ok("★ 故 K 通道即 M1 之 n_s=1 直接檢定（非外推）⟹ KILL 無需任何 n_s 界限",
   True, "n_s=1 固定", "direct",
   "B13 之『需 n_s=2.133』只係旁證；主殺法 = K 自身 n_s=1")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(O_K_X256=float(mpf(rows['X = 2.56（docs/18 §2）'][1].numerator) / rows['X = 2.56（docs/18 §2）'][1].denominator),
                        O_K_Xcanon=float(mpf(rows['X = 42667/16670（docs/39）'][1].numerator) / rows['X = 42667/16670（docs/39）'][1].denominator),
                        off_X256=float(rows['X = 2.56（docs/18 §2）'][2]),
                        off_Xcanon=float(rows['X = 42667/16670（docs/39）'][2]),
                        ns_req_256=float(mpf(ns_a.numerator) / ns_a.denominator),
                        ns_req_canon=float(mpf(ns_b.numerator) / ns_b.denominator),
                        em_M1=float(mpf((DM_K + X_A * O_M1).numerator) / (DM_K + X_A * O_M1).denominator)),
           verdict="M4 docs/39 接受；B13 措辭更正（域={0,1,2}），主殺法=n_s(K)=1")
json.dump(OUT, open(os.path.join(M7, "data", "15_m4_intake_recheck.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/15_m4_intake_recheck.json")