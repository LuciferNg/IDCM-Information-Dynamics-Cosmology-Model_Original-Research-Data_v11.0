"""M5 docs/51 配套 — isospin 通道在 v3 模型內之定量檢定（律 B 之算子層）
================================================================================
`docs/50` 立律 B：`Δ0` 之符號由輕夸克味決定（含 down-type 者較重）。本輪問：
**v3 模型之 isospin 通道能否給出 ~4.45 MeV 之量級？**

模型結構（`docs/41`／`docs/43`）：`Δ0^model = A_model − 3·c_model`，
`c ∝ 1/(m_h m_l)`（`docs/41` v2：`C_q = (8π/3)αq/(16 m_h m_l)`）⟹ 令 `m_l = m_u, m_d` 可得 isospin 分裂：
    `Δ0^iso = −3[c(m_l^N) − c(m_l^C)]`，`c(m_l) = c_cal·(m_l^cal/m_l)`，`m_l^cal = 0.30 GeV`

判準（事前寫死）：
  I1 用模型**已校準之 contact 振幅**（不得重擬合）：`|c_D| = 0.3525`（校準）、`c_B = 0.0775`、`c_K = 0.0135`。
  I2 逐族反解「所需 `δ ≡ m_d − m_u`」：`δ_req = |Δ0^data|·m_l/(3|c|)`。
  I3 若 `δ_req / δ_phys > 50`（任一族）⟹ isospin 通道**定量上死**（`δ_phys ≈ 2.5 MeV`，FLAG）。
  I4 另報「即使 `δ/m_l = 100%`」之最大可達 `Δ0 = 3|c|`，與所需比較 ⟹ 上界判（不依賴 `δ` 之數值）。
  I5 不得以增大 `c` 補救（`c` 已由 hyperfine 通道釘死，`docs/35`／`docs/41`）。
  I6 結論只作「夠／唔夠」，不作新參數引入。
輸出：`data/53_isospin_channel.json`
"""
import json

MU = 0.300          # GeV，模型之輕夸克質量（v3 校準）
DELTA_PHYS = 0.0025  # GeV = 2.5 MeV（FLAG 2024 m_d − m_u，±0.5）

# 模型 contact 振幅（MeV）：D 為校準目標，其餘為模型預言（docs/41 §V2）
C = {"π": None, "K": 0.0135, "D": 0.3525, "B": 0.0775}
D0 = {"π": 4.5936, "K": 3.9340, "D": 4.8200, "B": 0.3100}

print("=" * 100)
print("M5 docs/51 — v3 模型之 isospin 通道：定量檢定")
print("=" * 100)
print(f"\n[I1] 模型 contact 振幅（MeV）：D {C['D']}（校準）、B {C['B']}、K {C['K']}；m_l^cal = {MU} GeV")

print("\n[I2] 逐族反解所需 δ ≡ m_d − m_u")
rows = []
for f in ("K", "D", "B"):
    c = C[f]
    d_req = D0[f] * MU / (3 * c)          # GeV
    rows.append({"fam": f, "c": c, "d0": D0[f], "delta_req_GeV": d_req,
                 "ratio_to_phys": d_req / DELTA_PHYS})
    print(f"  {f}: |c| = {c:.4f} MeV ⟹ δ_req = {D0[f]:.4f}×{MU}/(3×{c:.4f}) = "
          f"{d_req*1000:9.2f} MeV ⟹ 相對 δ_phys({DELTA_PHYS*1000:.1f} MeV) 之 "
          f"{d_req/DELTA_PHYS:9.1f}×")

print("\n[I4] 上界判（不依賴 δ）：即使 δ/m_l = 100% ⟹ Δ0^max = 3|c|")
for f in ("K", "D", "B"):
    mx = 3 * C[f]
    print(f"  {f}: Δ0^max = 3×{C[f]:.4f} = {mx:.4f} MeV vs 所需 {D0[f]:.4f} ⟹ "
          f"{'🔴 上界已不足（' + format(D0[f]/mx, '.1f') + '×）' if mx < D0[f] else '界內'}")

print("\n[I3] 判")
worst = max(r["ratio_to_phys"] for r in rows)
print(f"  max(δ_req/δ_phys) = {worst:.1f} ⟹ "
      f"{'🔴 isospin 通道定量上死（>50×）' if worst > 50 else '未達判死門檻'}")
print("  註：δ_phys 之不確定度（FLAG ±0.5 MeV，20%）遠不足以補 10²–10⁴ 之缺口。")

print("\n[I5] 為何不可靠增大 c 補救")
print("  c 已由 hyperfine（ΔI=1）通道釘死：c_B 預言 0.07960 vs 數據 0.0775（+2.72%）、")
print("  c_D = 0.3525 為校準目標 ⟹ 加大 c 會同時破壞 ΔI=1 通道（該通道現時 1/6 命中、4/6 無預言）。")

print("\n[I6] 結論")
print("  ① 律 B（輕夸克味）之**符號**成立（K/D/B 3/3），但 v3 模型之 isospin 通道")
print("     比所需弱 1.6×10²–1.2×10⁴ ⟹ **量級上不可能**由其 contact 項承載。")
print("  ② 加上 `docs/49` 之 EM 扇區判死 ⟹ v3 模型對 `Δ0` **兩個候選機制皆無力**。")
print("  ③ ⟹ `Δ0` 之 4.45 MeV 標度必需由**新扇區**承載；該扇區之量子數已定：")
print("     **單體、odd under 輕夸克味（isospin）交換、標度 ≈ 3.70×A_D（= 2.51×t）**。")

out = {"rows": rows, "delta_phys": DELTA_PHYS, "max_ratio": worst,
       "verdict": "isospin_channel_dead" if worst > 50 else "inconclusive"}
with open("data/53_isospin_channel.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("\n-> data/53_isospin_channel.json")