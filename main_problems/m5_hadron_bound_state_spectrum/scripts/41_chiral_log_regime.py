"""M5 docs/39 配套 — K 側手徵對數（χlog）路線之零參數否證
================================================================================
路線 E2。`docs/37` 判 K 側只需 **+0.66 MeV**（`A_K^data − A_K^model = 1.3335 − 0.6740`）。
本輪問：χPT 手徵對數項能否提供此殘差？改以**物理常數**直接算，不引入任何擬合係數。

判準（事前寫死；本檔 docstring 即預註冊）：
  C1 輸入常數一律外部給定（PDG／FLAG）：`f_π = 92.07 MeV`、`f_K = 110.1 MeV`，`m_π±`、`m_K±`、`m_K⁰` 由 `data/pdg_mass_width_2024.mcd` 直讀。
  C2 χlog 形式取標準單圈：`δA_P ∝ (3α/4π)·(m_P²/f_P)·ln(m_P²/μ²)`。
  C3 因 `μ` 未定，**只判比值**（`μ` 部分抵消），不得用 `μ` 之手選值制造命中。
  C4 判：若 `|χlog 比值預言 / 數據比值| > 2` 或 `< 1/2` ⟹ χlog 路線 **no-go**（不可作 K 殘差之解釋）。
  C5 另一支：報 `|ΔM|`（`M(neutral) − M(charged)`）之通用尺度，檢查係否與 `m_h` 無關（π／K／D 同值、B 抑制）。
  C6 不得改動模型任何數值；本輪只判生死。
輸出：`data/41_chiral_log_regime.json`
"""
import json
import math
from pathlib import Path

ALPHA = 1 / 137.035999084
F_PI, F_K = 92.07, 110.1          # PDG / FLAG 2024（外部輸入，C1）
A_K_DATA, A_K_MODEL = 1.33350, 0.6740
A_D_DATA, A_D_MODEL = 1.20250, 1.20250

M = {}
for line in Path("data/pdg_mass_width_2024.mcd").read_text().splitlines():
    if line.startswith("*") or len(line) < 60:
        continue
    ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
    try:
        mv = float(line[33:51])
    except ValueError:
        continue
    for k in ids:
        if k:
            M[int(k)] = mv * 1000.0

print("=" * 100)
print("M5 docs/39 — K 側 χlog 路線零參數否證")
print("=" * 100)

m_pi, m_K, m_D, m_B = M[211], M[321], M[411], M[521]
print(f"\n[C1] 輸入：m_π± = {m_pi:.4f}、m_K± = {m_K:.4f}、m_D± = {m_D:.4f}、m_B± = {m_B:.4f} MeV")
print(f"     f_π = {F_PI}、f_K = {F_K} MeV（外部）")

print("\n[C2/C3] χlog 比值預言（μ 抵消）")
r_chi = (m_K ** 2 / F_K) / (m_pi ** 2 / F_PI)
print(f"     (m_K²/f_K)/(m_π²/f_π) = {r_chi:.4f}")
for mu in (0.5, 1.0, 2.0):
    lg_pi, lg_K = math.log(m_pi ** 2 / mu), math.log(m_K ** 2 / mu)
    print(f"     μ = {mu:.1f} GeV：ln(m_π²/μ²) = {lg_pi:+.3f}、ln(m_K²/μ²) = {lg_K:+.3f}"
          f" ⟹ 比值 = {r_chi * lg_K / lg_pi:+.3f}")

print("\n[C4] 與數據比（同通道：電荷分裂 |ΔM|）")
dM = {"π": M[111] - M[211], "K": M[311] - M[321], "D": M[421] - M[411], "B": M[511] - M[521]}
for k, v in dM.items():
    print(f"     |ΔM({k})| = {abs(v):8.4f} MeV（帶號 {v:+8.4f}）")
r_data = abs(dM["K"] / dM["π"])
chi_mu1 = r_chi * math.log(m_K ** 2 / 1.0) / math.log(m_pi ** 2 / 1.0)
ratio = chi_mu1 / r_data
print(f"     數據 |ΔM_K/ΔM_π| = {r_data:.4f}；χlog（μ = 1 GeV）= {chi_mu1:.4f} ⟹ 預言/數據 = {ratio:.2f}×")
verdict_K = abs(ratio) > 2
print(f"     ⟹ K 側 χlog 路線：{'🔴 no-go（超 'f'{ratio:.1f}×）' if verdict_K else '🟡 未否證'}")

print("\n[C4b] 殘差直接比：K 需 +0.66 MeV，χlog 之自然尺度")
resid = A_K_DATA - A_K_MODEL
for P, m, f in (("π", m_pi, F_PI), ("K", m_K, F_K)):
    for mu in (1.0,):
        scale = (3 * ALPHA / (4 * math.pi)) * (m ** 2 / f) * abs(math.log(m ** 2 / mu))
        print(f"     χlog 尺度（{P}，μ = 1 GeV）= {scale:8.3f} MeV")
print(f"     K 所需殘差 = {resid:.4f} MeV ⟹ 若注入 χlog(K) 會超 "
      f"{(3*ALPHA/(4*math.pi))*(m_K**2/F_K)*abs(math.log(m_K**2/1.0))/resid:.1f}×")
verdict_resid = (3 * ALPHA / (4 * math.pi)) * (m_K ** 2 / F_K) * abs(math.log(m_K ** 2 / 1.0)) / resid > 2

print("\n[C5] 通用尺度：|ΔM| 對 m_h 之依賴")
for k in ("π", "K", "D"):
    print(f"     {k}: |ΔM| = {abs(dM[k]):.4f} MeV")
vals = [abs(dM["π"]), abs(dM["K"]), abs(dM["D"])]
mean = sum(vals) / 3
var = sum((v - mean) ** 2 for v in vals) / 2
print(f"     π/K/D 均值 = {mean:.4f} ± {math.sqrt(var):.4f} MeV（離散 {100*math.sqrt(var)/mean:.1f}%）")
print(f"     B: |ΔM| = {abs(dM['B']):.4f} MeV ⟹ 相對 π/K/D 抑制 "
      f"{mean/abs(dM['B']):.1f}×（m_h 由 1.27 → 4.18，m 比 3.29）")

out = {"inputs": {"f_pi": F_PI, "f_K": F_K, "m_pi": m_pi, "m_K": m_K, "m_D": m_D, "m_B": m_B},
       "chi_ratio": {"mass_factor": r_chi, "mu1GeV": chi_mu1, "data": r_data, "pred_over_data": ratio},
       "residual": {"needed_MeV": resid, "chi_scale_K_MeV": (3 * ALPHA / (4 * math.pi)) * (m_K ** 2 / F_K) * abs(math.log(m_K ** 2 / 1.0))},
       "universal": {"dM": dM, "mean_piK_D": mean, "scatter": math.sqrt(var), "B_suppression": mean / abs(dM["B"])},
       "verdict": {"K_chi_log_nogo": bool(verdict_K and verdict_resid)}}
with open("data/41_chiral_log_regime.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=float)
print(f"\n[總結] K 側 χlog 路線 = {'🔴 NO-GO（兩支皆超 2×）' if out['verdict']['K_chi_log_nogo'] else '🟡 未定'}")
print("-> data/41_chiral_log_regime.json")