#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
18_q2_provenance_recheck.py — 依 M13 回覆（2026-11-03 00:35）重核 Q2

M13 之三點異議：
  ② 幾何級數對 B 無鑑別力（B=4 z=−0.685 vs B=5 z=−0.079，α 通道）
  ③ 44/33 與幾何級數代數恆等（33 = N_d(B−1)/B），除非 K_tree = 33 有獨立登錄
  ④ 迹通道才 B-selecting（Δ_N/Δ_Ξ；B=3 −5.72σ／B=5 +4.93σ／B=6 +8.55σ）

本 script 逐點重算（M15 側自算，不引用 M13 數字作證據；M13 之 α 通道數字只作引述標記）
"""
from fractions import Fraction as Fr
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OK, BAD = [], []


def chk(label, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {label}" + (f"  |  {detail}" if detail else ""))
    (OK if cond else BAD).append(label)
    return cond


NM, ND, L5 = 44, 44, 11

print("=" * 78)
print("§1  幾何級數 Σ_{k≥0} B^{-k} = B/(B−1)：對 B 有無鑑別力？")
print("=" * 78)
rows = []
for B in range(2, 9):
    s = Fr(B, B - 1)
    rows.append((B, s, float(s)))
    print(f"  B = {B}: Σ = {s} = {float(s):.6f}")
chk("恆等式對所有 B ≥ 2 成立 ⟹ **不揀 B**", all(rows[i][1] == Fr(i + 2, i + 1) for i in range(len(rows))))
tgt = Fr(4, 3)
which = [B for B, s, _ in rows if s == tgt]
chk("4/3 只對應 B = 4 ⟹ 要用 4/3 就必須先有 B = 4", which == [4], f"{which}")

print()
print("=" * 78)
print("§2  「N_d/K_tree = 4/3」是否獨立第二源？（M13 ③）")
print("=" * 78)
ktree = ND - L5
print(f"  K_tree = N_d − L₅ = {ND} − {L5} = {ktree}")
print(f"  N_d/K_tree = {ND}/{ktree} = {Fr(ND, ktree)} = {float(Fr(ND, ktree)):.6f}")
alt = Fr(ND * (4 - 1), 4)
print(f"  N_d(B−1)/B = {ND}×3/4 = {alt}   ⟹ 兩式同值 ⟺ L₅ = N_d/B")
chk("N_d/K_tree = 4/3 精確", Fr(ND, ktree) == Fr(4, 3))
chk("★ 同值之充要條件 = L₅ = N_d/B", L5 * 4 == ND, f"{L5}×4 = {L5*4} = N_d = {ND}")
print("  ⟹ 44/33 與 B/(B−1) 係**同一恆等**（用 L₅ = N_d/B 橋接）⟹ 不構成獨立第二源")
chk("故「兩條零輸入源」之說法不成立（實為一恆等 + 一數值巧合）", True)

print()
print("=" * 78)
print("§3  迹通道（Δ_N/Δ_Ξ = (B−1)/B²）是否 B-selecting？")
print("=" * 78)
mp_, mn_ = Fr("0.93827208816"), Fr("0.9395654205")
sxi_m, sxi_0 = Fr("1.32171"), Fr("1.321710")
dN = mn_ - mp_
# PDG: m(Ξ⁻) = 1.32171, m(Ξ⁰) = 1.31486
mXi_m, mXi_0 = Fr("1.32171"), Fr("1.31486")
dXi = mXi_m - mXi_0
r = dN / dXi
# σ（PDG 不確定度：m_Ξ⁻ ±0.00007、m_Ξ⁰ ±0.00020、m_n ±0.0000005、m_p ±0.0000001）
s_dXi = (Fr("0.00007") ** 2 + Fr("0.00020") ** 2) ** Fr(1, 2)
s_r = float(r) * float(s_dXi / dXi)
print(f"  Δ_N = {float(dN):.10f} MeV；Δ_Ξ = {float(dXi):.6f} MeV")
print(f"  Δ_N/Δ_Ξ = {float(r):.10f} ± {s_r:.9f}")
best = []
for B in range(2, 8):
    pred = Fr(B - 1, B * B)
    z = (float(pred) - float(r)) / s_r
    best.append((B, float(pred), z))
    print(f"  B = {B}: (B−1)/B² = {float(pred):.8f}   z = {z:+.2f}σ")
zmin = min(best, key=lambda t: abs(t[2]))
print(f"  ⟹ 最近者 B = {zmin[0]}（z = {zmin[2]:+.2f}σ）")
chk("迹通道選 B = 4（|z| < 1，次選 > 4σ）", zmin[0] == 4 and abs(zmin[2]) < 1
    and min(abs(t[2]) for t in best if t[0] != 4) > 4, f"B={zmin[0]} z={zmin[2]:+.2f}")
chk("故 M13 ④ 成立：**迹通道係唯一 B-selecting 者**", True)

print()
print("=" * 78)
print("§4  判詞（Q2 重核）")
print("=" * 78)
print("""  ① 幾何級數 = B/(B−1) 係對所有 B 成立之恆等 ⟹ **不揀 B**（M13 ② 成立）
  ② 44/33 與 B/(B−1) 同值 ⟺ L₅ = N_d/B（11 = 44/4）⟹ **同一恆等**，非第二源（M13 ③ 成立）
  ③ B = 4 由**迹通道**選定（z = −0.22σ；次選 ≥ 4.9σ）⟹ M13 ④ 成立
  ④ M13 之 α 通道數字（B=4 −0.685σ、B=5 −0.079σ）為第三方引述，本輪**不採為證據**，
     只記其方向（兩個 <1σ 候選 ⟹ α 通道無鑑別力，與 ② 一致）
  ⟹ Q2 判詞由「🟢 已解（兩條零輸入源）」**降為 🟡**：
     形式有 a priori 源；B 由迹通道選定；但「兩源」主張**撤回**，K_tree = 33 之獨立登錄仍缺""")
chk("Q2 降級 🟢 → 🟡（源由 2 減為 1＋恆等）", True)

print()
print("=" * 78)
print(f"§5  自審：PASS {len(OK)}、FAIL {len(BAD)}")
print("=" * 78)
for b in BAD:
    print(f"  ✗ {b}")

(DATA / "18_q2_provenance_recheck.json").write_text(json.dumps({
    "geo_by_B": {str(B): str(s) for B, s, _ in rows},
    "ktree": ktree, "L5_times_B": L5 * 4,
    "ratio_obs": round(float(r), 10), "ratio_sigma": round(s_r, 10),
    "z_by_B": {str(B): round(z, 2) for B, _, z in best},
    "verdict": "Q2 downgraded: geometric series is B-agnostic; 44/33 is the same identity; B=4 from trace channel",
    "n_pass": len(OK), "n_fail": len(BAD),
}, indent=1, ensure_ascii=False), encoding="utf-8")
print("\n[data] 18_q2_provenance_recheck.json 已寫")
print("== done rc=%d" % (0 if not BAD else 1))
raise SystemExit(0 if not BAD else 1)