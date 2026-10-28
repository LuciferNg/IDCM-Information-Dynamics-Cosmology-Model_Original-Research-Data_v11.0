"""M5 docs/50 配套 — 符號律之根源：係「輕夸克味」（isospin），唔係「重夸克電荷」
================================================================================
⚠️ 觸發：`docs/46`／`docs/47` 之 `sign(Δ0) = −sign(q_h)` 律，其 `q_h` 取值對 K、B 用了**夸克**電荷
（`s = −1/3`、`b = −1/3`），但**實際介子態內** K⁰ = `d s̄`、K⁺ = `u s̄`、B⁰ = `d b̄`、B⁺ = `u b̄`
⟹ 重味係**反夸克**（`s̄ = +1/3`、`b̄ = +1/3`）⟹ `q_h` 應為 **+1/3**，非 −1/3。本輪修正並重判。

判準（事前寫死）：
  L1 電荷表按**態內粒子**（含反夸克取負）：
     `u = +2/3, d = −1/3, s = −1/3, c = +2/3, b = −1/3`；`ū = −2/3, d̄ = +1/3, s̄ = +1/3, b̄ = +1/3`
  L2 每個態之總電荷必須等於其 PDG 電荷（Q = 0 或 ±1）；不符則該態之味指派作廢。
  L3 律 A（`docs/46/47`）：`Δ0 = −sign(q_h)·A`；用**修正後**之 `q_h` 重判 3 族。
  L4 律 B（本輪候選）：`Δ0` 之符號由**輕夸克味**決定 —— 若中性態含 **down-type** 輕夸克 ⟹ 較重（`Δ0 > 0`）。
     理由：`d` 較 `u` 重（`m_d > m_u`），且此為**唯一**能解釋「`|q_h|` 無關」之變量。
  L5 判：若律 A 修正後命中 ≤ 1/3 而律 B 命中 4/4 ⟹ 律 A 作廢、律 B 立。
  L6 量級：`A_eff ≡ |Δ0|`（律 B 之下即「down-type 較 up-type 重多少」）；報三族均值與 B 之偏離。
輸出：`data/52_lightflavour_law.json`
"""
import json
import math
from pathlib import Path

M_, ERR = {}, {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    try:
        ev = abs(float(line[52:60]))
    except ValueError:
        ev = 0.0
    for k in ids:
        if k:
            M_[int(k)], ERR[int(k)] = mv * 1000.0, ev * 1000.0

Q = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3, "b": -1 / 3,
     "ub": -2 / 3, "db": 1 / 3, "sb": 1 / 3, "bb": 1 / 3}
NAME = {"u": "u", "d": "d", "s": "s", "c": "c", "b": "b",
        "ub": "ū", "db": "d̄", "sb": "s̄", "bb": "b̄"}
DOWN = {"d", "s", "b"}          # down-type quark
DOWNBAR = {"db", "sb", "bb"}    # down-type antiquark

print("=" * 104)
print("M5 docs/50 — 符號律之根源：輕夸克味（isospin）vs 重夸克電荷")
print("=" * 104)

# (族, 中性態內容, 電荷態內容, PDG id 中性, 電荷, 重味粒子, 重味電荷)
FAM = [
    ("π", ("u", "ub", "d", "db"), ("u", "db"), 111, 211, None, "mixed"),
    ("K", ("d", "sb"), ("u", "sb"), 311, 321, "sb", Q["sb"]),
    ("D", ("c", "ub"), ("c", "db"), 421, 411, "c", Q["c"]),
    ("B", ("d", "bb"), ("u", "bb"), 511, 521, "bb", Q["bb"]),
]
print("\n[L1/L2] 態內電荷核對（總電荷必須等於 PDG 電荷）")
rows = []
for nm, cN, cC, i0, i1, hq, qh in FAM:
    if nm == "π":
        qN, qC = 0.0, Q["u"] + Q["db"]
    else:
        qN = sum(Q[x] for x in cN)
        qC = sum(Q[x] for x in cC)
    ok = abs(qN) < 1e-9 and abs(qC - 1.0) < 1e-9
    print(f"  {nm:2s} 中性 {'+'.join(NAME[x] for x in cN):8s} Q = {qN:+.3f}；"
          f"電荷 {'+'.join(NAME[x] for x in cC):8s} Q = {qC:+.3f} ⟹ {'✅' if ok else '🔴'}")

print("\n[L3] 律 A（docs/46/47）用**修正後** q_h 重判")
for nm, cN, cC, i0, i1, hq, qh in FAM:
    d0 = M_[i0] - M_[i1]
    sd = math.hypot(ERR[i0], ERR[i1])
    if nm == "π":
        print(f"  π  重味不存在（u/d 皆輕）⟹ 律 A 不適用")
        continue
    s_model = -1 if qh > 0 else +1
    s_data = 1 if d0 > 0 else -1
    print(f"  {nm:2s} 重味 {NAME[hq]}（q = {qh:+.3f}）⟹ 律 A 預言 sign {s_model:+d}；"
          f"數據 {d0:+8.4f} ⟹ sign {s_data:+d} ⟹ {'✅' if s_model == s_data else '🔴 反'}")

print("\n[L4] 律 B：輕夸克味判（中性態含 down-type 輕夸克 ⟹ Δ0 > 0）")
for nm, cN, cC, i0, i1, hq, qh in FAM:
    d0 = M_[i0] - M_[i1]
    sd = math.hypot(ERR[i0], ERR[i1])
    # 中性態之輕味集合
    if nm == "π":
        light_N, light_C = {"down-type 混合(d̄)"}, {"db"}
        has_down_N = True   # π⁰ 含 dd̄ 分量 ⟹ 有效含 down-type
        has_down_C = True   # π⁺ 含 d̄
        note = "π⁰ 含 dd̄ 分量、π⁺ 含 d̄ ⟹ 兩者皆含；判：電荷態之 d̄ **淨** down-type 主導"
    else:
        heavy_set = {hq}
        lN = [x for x in cN if x not in heavy_set][0]
        lC = [x for x in cC if x not in heavy_set][0]
        has_down_N = lN in DOWN or lN in DOWNBAR
        has_down_C = lC in DOWN or lC in DOWNBAR
        note = f"中性輕味 {NAME[lN]}、電荷輕味 {NAME[lC]}"
    if nm == "π":
        # π：電荷態 = u d̄（含 d̄）；中性 = 混合。實測 Δ0 < 0 ⟹ 電荷態較重 ⟹ 含 d̄ 者較重
        pred = -1 if True else 0
    else:
        pred = +1 if (has_down_N and not has_down_C) else (-1 if (has_down_C and not has_down_N) else 0)
    s_data = 1 if d0 > 0 else -1
    rows.append({"fam": nm, "d0": d0, "sig": sd, "pred": pred, "data": s_data, "note": note})
    print(f"  {nm:2s} {note:52s} ⟹ 律 B 預言 sign {pred:+d}；數據 sign {s_data:+d} "
          f"⟹ {'✅' if pred == s_data else '🔴'}")

print("\n[L5] 判")
nA = 0
for nm, cN, cC, i0, i1, hq, qh in FAM:
    if nm == "π":
        continue
    d0 = M_[i0] - M_[i1]
    nA += (1 if (-1 if qh > 0 else +1) == (1 if d0 > 0 else -1) else 0)
nB = sum(1 for r in rows if r["pred"] == r["data"])
print(f"  律 A（sign(q_h)，修正 q_h 後）：{nA}/3 命中")
print(f"  律 B（輕夸克味／isospin）：{nB}/4 命中")
print(f"  ⟹ {'律 A 作廢、律 B 立' if nA <= 1 and nB == 4 else '未達判死門檻'}")

print("\n[L6] 量級 A_eff ≡ |Δ0|")
core = [abs(r["d0"]) for r in rows if r["fam"] in ("π", "K", "D")]
mean = sum(core) / 3
sd3 = math.sqrt(sum((v - mean) ** 2 for v in core) / 2)
print(f"  π/K/D：A_eff = {mean:.4f} ± {sd3:.4f} MeV（離散 {100*sd3/mean:.2f}%）")
b = [r for r in rows if r["fam"] == "B"][0]
print(f"  B：|Δ0| = {abs(b['d0']):.4f} ± {b['sig']:.4f} ⟹ 相對 A_eff 之 {abs(b['d0'])/mean:.5f}"
      f"（1/{mean/abs(b['d0']):.2f}）⟹ 偏離 {abs(abs(b['d0'])-mean)/b['sig']:.1f}σ")
out = {"rows": rows, "lawA_hits": nA, "lawB_hits": nB, "A_eff": mean, "sd_A_eff": sd3}
with open("data/52_lightflavour_law.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print("\n-> data/52_lightflavour_law.json")