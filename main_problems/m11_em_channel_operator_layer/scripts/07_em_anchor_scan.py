"""M11 scripts/07 — em(K⁺−K⁰) 外部錨之 ε 來源掃描（FLAG 2024 Eqs. 38–46）
================================================================================================
問題：M11／M4 之 K 通道錨 em(K⁺−K⁰) 由 ε 組裝。M4 `docs/36` 只取 **MILC 18 單篇**
      ε = 0.78(1)(+8/−11)；但 FLAG 2024（arXiv:2411.04268）§4.1.3 自己給出平均：
      Eq. (45) Nf=2+1+1  ε = 0.79(6)（平均 RM123 17、MILC 18、BMW 14 之 (ΔM²_K)_γ）
      Eq. (46) Nf=2+1     ε = 0.73(17)（BMW 16A）
      且 FLAG 明言：若 QED sea 為 (ΔM²_K)_γ 之 O(10%)，則 ε 有 **~0.15 之 quenching 誤差**。

本腳本：以 Fraction 精確算各 ε 來源之 em 與 M4 預言（docs/29 之 +0.3327）之 z；
並算 M11 C2 之 ⟨O⟩_K。全部輸入逐字取自已下載之 FLAG 2024 PDF（md5 `3eedc5a7a51d3cbda1c460e00f1e2a9c`）
與 M4 `docs/36`／`docs/40` 所用之同一組 PDG 質量。

判準（寫死）：
  E1  em = (1+ε)·ΔM²_π/(M_K⁺+M_K⁰)（FLAG Eq. (39)，LO isospin-breaking ⟹ 只作 LO 分解）
  E2  σ_em = σ_ε·ΔM²_π/(M_K⁺+M_K⁰)
  E3  z = (pred − em)/sqrt(σ_em² + σ_X²)，σ_X 由 X 之誤差傳播（M4 docs/36 之 0.156667）
  E4  判詞穩健 ⟺ 所有 ε 來源皆同向且 |z| ≥ 3
輸出：data/07_em_anchor_scan.json
"""
import json
import os
import math

from fractions import Fraction as F

DMPI2 = F("1261.157")        # ΔM²_π = M²_π⁺ − M²_π⁰ (MeV²)，同 M4 docs/36
MKSUM = F("991.288")         # M_K⁺ + M_K⁰ (MeV)，PDG
PRED = F("0.3327")           # M4 docs/29 單一形式預言 em(K⁺−K⁰) (MeV)
SIG_X = F("0.156667")        # X 對 em 之誤差傳播（M4 docs/36）
X = F("2.56")                # X = m_d − m_u（FLAG 2024 Eq. (43)：4.70(5) − 2.14(8)）
DM_K = F("-3.9340")          # Δm(K⁺−K⁰) 原始質量差（PDG，M11 data/01）
COEFF_K = F("-1")

SRC = [
    ("Dashen（ε=0，參考）", F("0"), F("0"), "—"),
    ("MILC 18（M4 docs/36 現用）", F("0.78"), F("0.10"), "0.78(1)(+8/−11) 對稱化"),
    ("FLAG 2024 Eq.(45) Nf=2+1+1", F("0.79"), F("0.06"), "平均 RM123 17 + MILC 18 + BMW 14"),
    ("FLAG 2024 + QED quenching 0.15", F("0.79"), F("0.15"), "FLAG 原文自述之保守系統"),
    ("FLAG 2024 Eq.(46) Nf=2+1", F("0.73"), F("0.17"), "BMW 16A"),
]


def main():
    out = {"inputs": {"DMPI2": float(DMPI2), "MKSUM": float(MKSUM), "pred": float(PRED),
                      "sigma_X": float(SIG_X), "X": float(X)},
           "rows": []}
    print("=" * 112)
    print("M11 scripts/07 — em(K⁺−K⁰) 之 ε 來源掃描（FLAG 2024 Eqs. 38–46）")
    print("=" * 112)
    print(f"公式：em = (1+ε)·ΔM²_π/(M_K⁺+M_K⁰)，ΔM²_π = {float(DMPI2)} MeV²、M_K⁺+M_K⁰ = {float(MKSUM)} MeV")
    print(f"預言：M4 docs/29 之 em = {float(PRED)} MeV｜σ_X = {float(SIG_X)} MeV")
    print()
    print(f"{'ε 來源':34s}{'ε':>7s}{'em (MeV)':>11s}{'σ_em':>9s}{'z(ε only)':>11s}{'z(ε+X)':>10s}{'⟨O⟩_K':>9s}")
    for lab, eps, sig, note in SRC:
        em = (1 + eps) * DMPI2 / MKSUM
        sem = sig * DMPI2 / MKSUM
        st = F(math.sqrt(float(sem ** 2 + SIG_X ** 2))) if sem > 0 else SIG_X
        z_eps = None if sem == 0 else (PRED - em) / sem
        z_tot = (PRED - em) / st
        oK = (DM_K - em) / (COEFF_K * X)
        print(f"{lab:34s}{float(eps):7.2f}{float(em):11.4f}{float(sem):9.4f}"
              f"{('—' if z_eps is None else f'{float(z_eps):+11.2f}')}{float(z_tot):+10.2f}{float(oK):9.4f}")
        out["rows"].append(dict(label=lab, note=note, eps=float(eps), sigma_eps=float(sig),
                                em_MeV=float(em), sigma_em=float(sem),
                                z_eps=None if z_eps is None else float(z_eps),
                                z_total=float(z_tot), O_K=float(oK)))
    ems = [r["em_MeV"] for r in out["rows"][1:]]
    zs = [r["z_total"] for r in out["rows"][1:]]
    print()
    print(f"⟹ em 範圍（除 Dashen 參考）：[{min(ems):.4f}, {max(ems):.4f}] MeV；預言 {float(PRED)} ⟹ 差 "
          f"{float(PRED)/max(ems):.3f}–{float(PRED)/min(ems):.3f} 倍")
    print(f"⟹ z(ε+X) 範圍：[{min(zs):+.2f}, {max(zs):+.2f}]σ ⟹ "
          f"{'✅ 全部 |z| ≥ 3 ⟹ KILL 判詞穩健' if all(abs(z) >= 3 for z in zs) else '⚠ 有 |z| < 3 ⟹ 判詞鬆動'}")
    print(f"⟹ M11 之 ⟨O⟩_K 變動：{float(out['rows'][1]['O_K']):.4f} → {float(out['rows'][2]['O_K']):.4f}"
          f"（{float(out['rows'][2]['O_K']/out['rows'][1]['O_K']-1)*100:+.3f}%）")
    print("\n[C] M11 C9 之 K 通道是否用 ε？")
    E_K = DM_K - COEFF_K * X
    print(f"    E_K = Δm − coeff·X = {float(DM_K)} − ({float(COEFF_K)})×{float(X)} = {float(E_K):.4f} MeV"
          f"（純 PDG ＋ FLAG Eq.(43) 夸克質量，**不含 ε**）")
    print(f"    t = E_K/k₀ = {float(E_K/K0):.4f}（k₀ = {float(K0):.4f} MeV）；箱 [0, +0.6667] ⟹ 符號封鎖")
    print("    ⟹ ε 之任何選擇都不改變 M11 之 K 判詞（亦與 κ 無關）")
    _out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "07_em_anchor_scan.json")
    with open(_out, "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n-> data/07_em_anchor_scan.json")


K0 = F("3.8152622441")     # data/01 正典 k₀（R6 重算後）
main()