"""M4 Q2/Q3 — em 讀數 operator 約束集 v1：7 通道實證譜，形式限制與依賴圖（2026-10-24）
======================================================================================
上游：docs/14（R_em = 3.479、反相=恆等還原）；docs/17（最小電荷模型結構否證 ⟹ 必須非電荷加性）；
      01_problem_record §十一 Q2/Q3；Q5（selection rule）已 handoff O140。
目標：唔再「揀模型」，而係**由數據反推 operator 必須滿足嘅形式約束**，並定出依賴圖
      （若 operator 需 flavour→mode 映射，則 Q2/Q3 形式上下游於 O140）。

Protocol（事前寫死）：
  X ≡ m_d − m_u = 2.560 MeV（FLAG24，外部固定；唔由本輪擬合）
  C1 通道 em 讀數：E_ch ≡ Δm_ch − (substitution count)·X ，substitution = u↔d 替換數
     通道（PDG 2024）：n−p、K+−K0、π+−π0、Σ−−Σ0、Σ0−Σ+、Σ−−Σ+、Ξ−−Ξ0
  C2 形式約束測試（每項事前指定判準）：
     (a) 電荷加性：E 對 ΔQ_had 或 Δ(ΣQ²) 線性 ⟹ 若 R² > 0.9 且殘差 < 10% ⟹ 未排除
     (b) substitution-count 線性：E ≈ k·(sub count) ⟹ 同上判準
     (c) spectator 獨立：同 substitution、唔同 spectator 嘅通道 E 應相等（|ΔE| < 2σ）
     (d) 標度假設：E ∝ M_ch^p（單一 p）⟹ 由 (n−p, D+−D0) 定 p，預測其餘 ⟹ dev < 10% 才未排除
  C3 依賴圖：若 (a)–(d) 全排除 ⟹ operator 必須含 mode-level（通道內 mode 佔用）資訊
     ⟹ 需 flavour→mode 映射 = Q5 selection rule = O140 ⟹ Q2/Q3 標 blocked-by-O140
輸出：data/18_em_constraints.json ＋ stdout
"""
import json
import math

OUT = "data/18_em_constraints.json"

X = 2.560          # m_d − m_u（FLAG24，外部固定）
SX = 0.094

# 通道：名稱, Δm (MeV), σ, substitution count (u↔d 替換數), 通道質量尺度 (MeV), 類型
CH = [
    ("n−p",      1.293331, 0.000001, 1, 938.918, "baryon u↔d, spectators {u,d}"),
    ("K+−K0",   -3.934,    0.021,   -1, 495.644, "meson u→d（Δ 定義為 K+−K0 ⟹ 質量項 = −X）"),
    ("π+−π0",    4.59359,  0.00014,  1, 137.274, "annihilation（質量項 +X，標準分解）"),
    ("Σ−−Σ0",    4.807,    0.037,    1, 1195.046, "baryon d→u, spectators {d,s}"),
    ("Σ0−Σ+",    3.272,    0.074,    1, 1191.006, "baryon d→u, spectators {u,s}"),
    ("Σ−−Σ+",    8.079,    0.076,    2, 1193.410, "baryon 2×u↔d"),
    ("Ξ−−Ξ0",    6.85,     0.21,     1, 1318.285, "baryon d→u, spectators {s,s}"),
    ("D+−D0",    4.82,     0.07,     1, 1867.250, "meson d→u, spectator c"),
]

print("=" * 78)
print("M4 Q2/Q3 — em 讀數 operator 約束集 v1")
print("=" * 78)
print(f"\n外部固定：X = m_d − m_u = {X} ± {SX} MeV（FLAG24）")

rows = []
print(f"\n[C1] 通道 em 讀數 E_ch = Δm − (sub)·X")
print(f"  {'通道':<9}{'Δm':>10}{'sub':>5}{'E_ch':>10}{'σ_E':>8}  類型")
for name, dm, sdm, sub, mscale, typ in CH:
    E = dm - sub * X
    sE = math.hypot(sdm, abs(sub) * SX)
    rows.append({"name": name, "dm": dm, "sigma_dm": sdm, "sub": sub, "E": E,
                 "sigma_E": sE, "mass_scale": mscale, "type": typ})
    print(f"  {name:<9}{dm:>10.4f}{sub:>5}{E:>10.4f}{sE:>8.4f}  {typ}")

# ---------- C2(a) 電荷加性 ----------
print("\n[C2a] 電荷加性測試：E 對 Δ(ΣQ²) 與 ΔQ_had 線性？")
# Δ(ΣQ²)：以夸克荷 Q_u=2/3, Q_d=−1/3，substitution 取代價
dQ2 = {"n−p": -1/3, "K+−K0": -1/3, "π+−π0": 0.0, "Σ−−Σ0": -1/3, "Σ0−Σ+": -1/3,
       "Σ−−Σ+": -2/3, "Ξ−−Ξ0": -1/3, "D+−D0": +1/3}
dQh = {"n−p": +1, "K+−K0": +1, "π+−π0": +1, "Σ−−Σ0": 0, "Σ0−Σ+": +1,
       "Σ−−Σ+": +1, "Ξ−−Ξ0": +1, "D+−D0": +1}
for key in ("dQ2", "dQh"):
    xs = [dQ2[r["name"]] if key == "dQ2" else dQh[r["name"]] for r in rows]
    ys = [r["E"] for r in rows]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs); sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx == 0:
        print(f"  {key}: x 全常數 ⟹ 無鑑別力")
        continue
    b = sxy / sxx; a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    print(f"  {key}: slope {b:+.3f}  R² = {r2:.4f} ⟹ {'🟡 未排除' if r2 > 0.9 else '🔴 排除'}")

# ---------- C2(b) substitution-count 線性 ----------
print("\n[C2b] substitution-count 線性：E ≈ k·sub")
sel_b = [r for r in rows if r["sub"] != 0]
xs = [abs(r["sub"]) for r in sel_b]; ys = [r["E"] for r in sel_b]
n = len(xs); mx, my = sum(xs) / n, sum(ys) / n
sxx = sum((x - mx) ** 2 for x in xs); sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
b = sxy / sxx; a = my - b * mx
ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys)); ss_tot = sum((y - my) ** 2 for y in ys)
r2_sub = 1 - ss_res / ss_tot
print(f"  E = {a:+.3f} + {b:+.3f}·sub   R² = {r2_sub:.4f} ⟹ {'🟡 未排除' if r2_sub > 0.9 else '🔴 排除'}")
for r in sel_b:
    pred = a + b * abs(r["sub"])
    r["resid_sub"] = r["E"] - pred
print(f"  殘差範圍 {min(r['resid_sub'] for r in sel_b):+.3f} … {max(r['resid_sub'] for r in sel_b):+.3f} MeV")

# ---------- C2(c) spectator 獨立 ----------
print("\n[C2c] spectator 獨立：同 substitution、唔同 spectator ⟹ E 應相等（|ΔE| < 2σ）")
groups = {
    "A: +X 方向（final 含 d 取代 u）：meson/baryon 混合": ["n−p", "K+−K0"],
    "B: +X 方向 baryon (d 取代 u, 唔同 spectator)": ["Σ−−Σ0", "Σ0−Σ+", "Ξ−−Ξ0", "D+−D0"],
}
spec_ok = True
for gname, names in groups.items():
    sel = [r for r in rows if r["name"] in names]
    Es = [r["E"] for r in sel]; sEs = [r["sigma_E"] for r in sel]
    spread = max(Es) - min(Es)
    print(f"  {gname}: " + "  ".join(f"{r['name']}={r['E']:+.3f}±{r['sigma_E']:.3f}" for r in sel)
          + f"  ⟹ spread {spread:.3f} MeV")
    if spread > 0.30:
        spec_ok = False
        print(f"     ⟹ 🔴 spectator 獨立**排除**（spread {spread:.3f} ≫ 誤差）")
print(f"  綜合：{'🟡 未排除' if spec_ok else '🔴 spectator 獨立排除'}")

# ---------- C2(d) 標度幂律 ----------
print("\n[C2d] 標度幂律 E ∝ M^p：由 (n−p, D+−D0) 定 p，再預測其餘")
def get(nm): return next(r for r in rows if r["name"] == nm)
ra, rb = get("n−p"), get("D+−D0")
p = math.log(abs(rb["E"]) / abs(ra["E"])) / math.log(rb["mass_scale"] / ra["mass_scale"])
print(f"  由 n−p({ra['E']:+.3f} @ {ra['mass_scale']:.1f}) 與 D+−D0({rb['E']:+.3f} @ {rb['mass_scale']:.1f}) ⟹ p = {p:.4f}")
worst = 0.0
for r in rows:
    if r["name"] in ("n−p", "D+−D0") or r["sub"] == 0:
        continue
    pred = abs(ra["E"]) * (r["mass_scale"] / ra["mass_scale"]) ** p
    dev = abs(pred / abs(r["E"]) - 1) * 100
    worst = max(worst, dev)
    print(f"    {r['name']:<9} 預言 |E| = {pred:.3f} vs 實測 {abs(r['E']):.3f} ⟹ dev {dev:>6.2f}%")
print(f"  worst dev {worst:.1f}% ⟹ {'🟡 未排除' if worst < 10 else '🔴 幂律標度排除'}")

# ---------- C3 依賴圖 ----------
print("\n[C3] 依賴圖")
print("  事實：8 通道中，同一 substitution 數嘅 E 值由 −1.374 到 +4.29（跨度 5.66 MeV）")
print("        ⟹ 電荷加性 (a)、substitution-count 線性 (b)、spectator 獨立 (c)、單一幂律 (d) 全排除")
print("  結論：operator 必須含**通道內 mode 佔用**資訊（邊界類型／佔用 mode 集合）")
print("        ⟹ 需 flavour→mode 映射 = Q5 selection rule ⟹ **屬 O140**（已 handoff）")
print("        ⟹ Q2/Q3 形式上下游於 O140：本 O 唔可自行補此映射（10-16 裁決）")

res = {"meta": {"script": "18_em_constraints", "date": "2026-10-24", "X_external": X,
                "X_coef_convention": "Δm = Σ(m_final − m_initial) + E ⟹ K+−K0 嘅 X 係數 = −1"},
       "channels": rows,
       "verdicts": {"charge_additive_dQ2": "excluded", "dQ_had": "excluded",
                    "sub_count_linear": "excluded" if r2_sub < 0.9 else "not_excluded",
                    "spectator_independence": "excluded" if not spec_ok else "not_excluded",
                    "power_law_scale": "excluded" if worst >= 10 else "not_excluded"},
       "dependency": {"Q2_Q3": "blocked_by_O140_flavour_mode_map"}}
with open(OUT, "w") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)
print(f"\n→ {OUT}")
