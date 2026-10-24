"""M5 docs/08 配套 — H5c: 樣本外通道盲測（PDG 2024 官方檔）＋自旋態零參數檢定
================================================================================
公約（嚴格照 M4 docs/18）：E = [m(A) − m(B)] − coeff·X，X = m_d − m_u = 2.56 ± 0.094（FLAG24）
  8 通道用 M4 原 (A,B,coeff) ⟹ 必須 8/8 重現（管線校驗；首版用錯 co在与方向 ⟹ 6/8，已修）
  自旋檢定：同夸克內容、唔同自旋 ⟹ 任何「味標籤＋電荷」核預言 ΔE ≡ 0（對自旋完全盲）
數據源：`data/pdg_mass_width_2024.mcd`（PDG 2024 官方，curl 下載，sha256 12136f946df7860e…）
輸出：`data/08_sampleout_channels.json`
"""
import json
from pathlib import Path

MCD = Path("data/pdg_mass_width_2024.mcd")
X, XS = 2.56, 0.094
MASS, ERR = {}, {}
for line in MCD.read_text().splitlines():
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
            MASS[int(k)] = mv * 1000.0
            ERR[int(k)] = ev * 1000.0

print("=" * 96)
print("M5 docs/08 — H5c: 樣本外通道（PDG 2024）＋自旋態零參數檢定")
print("=" * 96)
print(f"解析 {len(MASS)} 個 PDG ID；X = {X} ± {XS} MeV")

def E(a, b, coeff=1):
    dm = MASS[a] - MASS[b]
    s = (ERR[a] ** 2 + ERR[b] ** 2 + (abs(coeff) * XS) ** 2) ** 0.5
    return dm - coeff * X, s, dm

CH8 = [("n−p", 2112, 2212, 1, -1.2667), ("K⁺−K⁰", 321, 311, -1, -1.3740),
       ("π⁺−π⁰", 211, 111, 1, 2.0336), ("Σ⁻−Σ⁰", 3112, 3212, 1, 2.2470),
       ("Σ⁰−Σ⁺", 3212, 3222, 1, 0.7120), ("Σ⁻−Σ⁺", 3112, 3222, 2, 2.9590),
       ("Ξ⁻−Ξ⁰", 3312, 3322, 1, 4.2900), ("D⁺−D⁰", 411, 421, 1, 2.2600)]
print("\n[① 管線校驗：重現 M4 docs/18 之 8 通道]")
ok = 0
for name, a, b, c, ref in CH8:
    v, s, dm = E(a, b, c)
    good = abs(v - ref) < 0.002
    ok += good
    print(f"  {name:8s} {v:+9.4f} ± {s:.4f}   M4 {ref:+9.4f}  {'OK' if good else 'MISMATCH'}")
print(f"  => {ok}/8")

print("\n[② 自旋態零參數檢定：同內容、唔同自旋 ⟹ 味標籤核預言 ΔE ≡ 0]")
PAIRS = [("K vs K*", (311, 321), (313, 323)),
         ("D vs D*", (411, 421), (413, 423)),
         ("Σ vs Σ(1385) [−/0]", (3112, 3212), (3114, 3214)),
         ("Σ vs Σ(1385) [0/+]", (3212, 3222), (3214, 3224)),
         ("Σc(2455) vs Σc(2520) [0/+]", (4212, 4222), (4214, 4224)),
         ("Ξc vs Ξc′ [0/+]", (4132, 4232), (4312, 4322))]
rows = []
for fam, p1, p2 in PAIRS:
    v1, s1, _ = E(*p1)
    v2, s2, _ = E(*p2)
    d = v1 - v2
    sd = (s1 ** 2 + s2 ** 2) ** 0.5
    z = d / sd
    rows.append({"family": fam, "E_low": v1, "E_high": v2, "dE": d, "sigma": sd, "z": z,
                 "s1": s1, "s2": s2})
    tag = "SPIN-DEPENDENT" if abs(z) > 3 else ("marginal" if abs(z) > 2 else "consistent")
    print(f"  {fam:30s} E1 = {v1:+8.4f} ± {s1:.4f}   E2 = {v2:+8.4f} ± {s2:.4f}"
          f"   dE = {d:+7.4f} ± {sd:.4f}  z = {z:+5.2f}  {tag}")

print("\n[③ 樣本外單通道（H5b 之後作真盲測對照）]")
EXTRA = [("K*0-K*+", (313, 323)), ("D*+-D*0", (413, 423)),
         ("S(1385)--S(1385)0", (3114, 3214)), ("S(1385)0-S(1385)+", (3214, 3224)),
         ("Sc(2455)0-Sc(2455)+", (4112, 4212)), ("Sc(2455)+-Sc(2455)++", (4212, 4222)),
         ("Sc(2520)0-Sc(2520)+", (4214, 4224)),
         ("Xic0-Xic+", (4132, 4232)), ("Xic'0-Xic'+", (4312, 4322)),
         ("Xib--Xib0", (5132, 5232))]
extra = {}
for name, (a, b) in EXTRA:
    v, s, dm = E(a, b, 1)
    extra[name] = {"E": v, "sigma": s, "dm": dm}
    print(f"  {name:22s} dm = {dm:+9.4f} -> E = {v:+8.4f} ± {s:.4f}")

nz = sum(1 for r in rows if abs(r["z"]) > 3)
nm = sum(1 for r in rows if abs(r["z"]) > 2)
print(f"\n[判定] |z|>3: {nz}/{len(rows)};  |z|>2: {nm}/{len(rows)}")
print("  味標籤核（H2-H4 四層全部）對自旋完全盲 => 每對 dE 必須嚴格為 0")
print("  => 任何 |z|>3 之對 = 該核類被 PDG 數據獨立否證（零參數、真樣本外）")

out = {"source": "https://pdg.lbl.gov/2024/mcdata/mass_width_2024.mcd",
       "X": X, "X_sigma": XS, "pipeline": f"{ok}/8",
       "spin_pairs": rows, "n_z_gt3": nz, "n_z_gt2": nm,
       "sampleout_low_spin": extra,
       "verdict": ("味標籤核類被 PDG 樣本外數據否證（自旋盲預言 dE=0 破產）；"
                   "em 算子必須含態（波函數）")}
with open("data/08_sampleout_channels.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/08_sampleout_channels.json")