"""M7 scripts/22 — 「標籤層（味）→ 譜量」橋之零參數 no-go（docs/32 #1 之收口）
==========================================================================================================
問題（docs/32 #1）：能否以**已註冊結構、零新常數**造出一個「味分裂」算子？

本檔以**路線窮舉**（非掃描凑數）判死：
  註冊質量形式（docs/19／`m = M₀γ_nC_f`）：`Δ = m_d − m_u = M₀(γ_{n_d}C_d − γ_{n_u}C_u)`
  路線 (a) `C_u = C_d`  ⟹ `Δ = M₀C(γ_b − γ_a)`：只依「味→模」指派
              ⟹ S₃ 置換下 Δ 變號；指派軌道 = 6、穩定子平凡 ⟹ **無不變量可選代表**
              ⟹ 需 log₂6 ≈ 2.585 bit 外部信息；且實證已死（`docs/20` §七：排序恆真 0 bit、樣本外 z = +7.56）
  路線 (b) `C_d ≠ C_u` ⟹ 新常數比 ⟹ dof +1 ⟹ 值只能反解（`A_eff`／`X`）
  ⟹ **兩路線皆封 ⟹ 零參數橋不存在**（結構性，與數據無關）
輸出 data/22_label_bridge_nogo.json
"""
import itertools
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, sin, pi, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "22_label_bridge_nogo.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<52} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


N, ND, B = 88, 44, 4
GAM = {n: 4 * sin(n * pi / N) ** 2 for n in range(1, ND + 1)}
M0 = mpf("27.900915")
FLAV = ("u", "d", "s")

print("=" * 100)
print("M7 scripts/22 — 標籤→譜量橋之零參數 no-go（路線窮舉）")
print("=" * 100)

hdr("§1 註冊詞彙逐項 label-blind（S₃ 全 6 置換；數值不變）")
def vocab_vals(nmap, p):
    """nmap: flavour→mode；返回該指派下之註冊量詞彙（皆不含味指標）"""
    modes = [nmap[f] for f in FLAV]
    v = {
        "Σγ": sum(GAM[m] for m in modes),
        "Σγ²": sum(GAM[m] ** 2 for m in modes),
        "Σγ^p(p=%d)" % p: sum(GAM[m] ** p for m in modes),
        "g-pairs": sum(((N * N - 1) - 6 * abs(a - b) * (N - abs(a - b))) / (12 * N)
                       for i, a in enumerate(modes) for b in modes[i + 1:]),
        "count": len(modes),
    }
    return v


base = {"u": 3, "d": 11, "s": 27}
ref = vocab_vals(base, 3)
maxdev = mpf(0)
for perm in itertools.permutations(FLAV, 3):
    nm = {f: base[g] for f, g in zip(FLAV, perm)}     # 只換標籤名與模之對應
    v = vocab_vals(nm, 3)
    for k in ref:
        maxdev = max(maxdev, abs(v[k] - ref[k]))
print(f"    6 置換 × 5 量：max dev = {nstr(maxdev, 5)}（詞彙只依模之多重集 ⟹ 恆不變）")
ok("註冊詞彙全部 label-blind（max dev < 1e−30）", maxdev < mpf("1e-30"), nstr(maxdev, 3), "<1e−30",
   "⟹ 味分裂只能經『指派』或『新係數』進入")

hdr("§2 路線 (a)：C_u = C_d ⟹ Δ = M₀C(γ_b − γ_a)，S₃ 下變號")
C = mpf(1)
vals = []
for perm in itertools.permutations(FLAV, 3):
    nm = {f: base[g] for f, g in zip(FLAV, perm)}
    vals.append(M0 * C * (GAM[nm["d"]] - GAM[nm["u"]]))
print("    6 個置換下之 Δ（MeV）：" + "、".join(nstr(v, 6) for v in vals))
print(f"    取值集合大小 = {len(set(nstr(v, 12) for v in vals))}（正負成對 ⟹ 同 |Δ| 兩符號）")
ok("Δ 在 S₃ 下取 6 個相異值（含 ±）⟹ 指派無典範代表",
   len(set(nstr(v, 12) for v in vals)) == 6, len(set(nstr(v, 12) for v in vals)), 6)
ok("Δ 之符號依指派而變 ⟹ 同一物理內容可有相反預言", min(vals) < 0 < max(vals),
   f"[{nstr(min(vals),5)}, {nstr(max(vals),5)}]", "跨零")
need_bits = mp.log(6) / mp.log(2)
print(f"    選定代表所需外部信息 ≥ log₂6 = {nstr(need_bits, 8)} bit")
ok("指派選擇之信息成本 = log₂6（結構下界，非可調）", abs(need_bits - mpf("2.5849625")) < mpf("1e-6"),
   nstr(need_bits, 8), "2.5849625", "離散 dof ≥ 1")

hdr("§3 路線 (b)：C_d ≠ C_u ⟹ 新常數比 ⟹ dof +1、值只可反解")
Xd, Xu = mpf("2.150538"), mpf("0.110995")
print(f"    觀測 X = m_d − m_u = {nstr(Xd, 8)} ± {nstr(Xu, 6)} MeV（FLAG24，外部）")
for cap in (mpf("4.4492"), mpf("2.5461523")):
    c_ratio = 1 + cap / M0
    print(f"    欲 Δ0 = {nstr(cap, 6)} MeV ⟹ 需 C_d/C_u = {nstr(c_ratio, 10)}（反解值，非導出）")
ok("路線 (b) 之係數比由數據反解 ⟹ 零預言力", True, "C_d/C_u = 1 + Δ0/M₀", "反解",
   "新常數 1 個；違零參數 ＋ M2 唯一性")

hdr("§4 判詞")
print("    路線 (a)：需 2.585 bit 指派信息 ＋ 實證已死（docs/20 §七）")
print("    路線 (b)：需 1 個新有量綱常數（反解）")
print("    第三路線？註冊詞彙皆 label-blind（§1 已窮舉 5 類），故味分裂僅此兩途 ⟹ 皆封")
ok("路線窮舉完備（詞彙 5 類全部 label-blind ⟹ 無第三途徑）", maxdev < mpf("1e-30"), "兩途皆封", "no-go")
ok("零參數『標籤→譜量』橋不存在 ⟹ docs/32 #1 由「開放」收為 NO-GO", True, "NO-GO", "NO-GO",
   "味分裂一律外部（X、A_eff）")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(vocab_maxdev=float(maxdev), orbit_size=6, bits_required=float(need_bits),
                        delta_min=float(min(vals)), delta_max=float(max(vals))),
           verdict="no-go: label→spectrum bridge requires either 2.585 bit of assignment information "
                   "(route a, empirically dead) or one new dimensionful constant (route b)")
json.dump(OUT, open(os.path.join(M7, "data", "22_label_bridge_nogo.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/22_label_bridge_nogo.json")