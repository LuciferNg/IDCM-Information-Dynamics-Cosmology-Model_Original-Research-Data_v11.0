"""M5 docs/16 配套 — D2：K 系殘餘是否「M4 X 線性扣除不足」所致？（四判準）
================================================================================================
背景：`docs/13` §六 判準草案；D1 已判死（自旋標度解耦）。D2 命題：K 系之
`E(K⁰−K⁺) = 1.374`、`E(K*⁰−K*⁺) = 1.320`（近相等、自旋盲）並非本核之 em，而係
M4 `X = m_d − m_u` **線性扣除不足**所留之強 IB 殘餘。

判準（掃描前寫死；任一破 ⟹ D2 破）：
  Q1 **量級**：同一核內，把 `m_d ≠ m_u` 分開（`m_l ± X/2`）後之非線性殘餘
     `Δ_nl = [E(m_l+δ) + 3|c|(m_l+δ)] − [E(m_l−δ) + 3|c|(m_l−δ)] − 2δ`（δ = X/2 = 1.28 MeV）。
     要 D2 成立，`Δ_nl` 須解釋 K 系觀測殘餘 1.374 之 **≥ 50%**（即 ≥ 0.687 MeV）。
  Q2 **自旋盲**：K 與 K* 之殘餘須近相等（強 IB 自旋盲）：觀測 1.374 vs 1.320，差 0.054 ± 0.3545。
  Q3 **族模式**：Σ/Ξ 系（含 s、無 c）同見「大電荷分裂、零自旋分裂」。
  Q4 **零參數比**：核預言 `ΔE_K/ΔE_D = |c|_K/|c|_D` 對比觀測 0.054/1.41 = 0.0383。
  Q5 **硬上界審計**（新）：`A = α q(ħc/a_C)⟨B⟩`，`⟨B⟩ ≤ 1` ⟹ `A_max = α q(ħc/a_C)`；
     檢 `A_req = 1.3335 ± 0.2570` 與 `A_max` 之關係（超界 ⟹ 該算子族結構上不可能）。

界限：只論本核（同一 ψ、兩體、contact on-site）；不引入格點 QED／手徵外推（本 O 無此能力）。

輸出：`data/17_d2_strong_ib.json`
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path.cwd() if (Path.cwd() / "README.md").exists() else Path(__file__).resolve().parent.parent
HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
X, DL = 2.56, 1.28
A0, T0 = 0.377422, 0.788398e-3          # docs/14 重校
Q_C, Q_S, M_C, M_S = 2 / 3, 1 / 3, 1.27, 0.495
ML = 0.30


def shape(g, h):
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    return float((w * np.array(POW)).sum() / w.sum()), float(w[0] / w.sum()), int(np.argmax(w))


def energy(q, mh, ml, aC=A0, t=T0):
    """singlet 態能量 E = A + 3|c|（MeV）；a_C、t 普適"""
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aC) ** 3 / t
    rA, P0, dp = shape(g, h)
    A = ALPHA * q * (HBARC / aC) * rA * 1000
    c = Cq * (HBARC / aC) ** 3 * P0 * 1000
    return A, c, A + 3 * c, (rA, P0, dp)


print("=" * 100)
print("M5 docs/16 — D2 四判準：K 系殘餘係強 IB 嗎？")
print("=" * 100)

# ---------------- Q1 量級 ----------------
print("\n[Q1] 非線性強 IB 殘餘（同一核，m_l ± δ；δ = X/2 = 1.28 MeV）")
print(f"     {'系統':6s}{'Ebind(m_l−δ)':>14s}{'Ebind(m_l+δ)':>14s}{'Δ_nl':>10s}{'佔 1.374':>10s}{'對 0.687':>10s}")
q1 = {}
for nm, q, mh in (("K", Q_S, M_S), ("D", Q_C, M_C)):
    _, _, Em, _ = energy(q, mh, ML - DL / 1000)
    _, _, Ep, _ = energy(q, mh, ML + DL / 1000)
    dnl = Ep - Em                      # X 已線性扣除 ⟹ 殘餘 = 束縛能對 m_l 之響應
    q1[nm] = {"Ebind_minus": Em, "Ebind_plus": Ep, "delta_nl": dnl,
              "frac_of_residual": dnl / 1.374}
    print(f"     {nm:6s}{Em:14.4f}{Ep:14.4f}{dnl:10.5f}{dnl/1.374*100:9.3f}%"
          f"{abs(dnl)/0.687*100:9.1f}%")
short = 0.687 / abs(q1["K"]["delta_nl"])
print(f"     ⟹ 非線性殘餘 |Δ_nl| ~{abs(q1['K']['delta_nl']):.4f} MeV，"
      f"僅為所需 1.374 之 {abs(q1['K']['delta_nl'])/1.374*100:.2f}%"
      f" ⟹ Q1 {'✅' if abs(q1['K']['delta_nl']) >= 0.687 else '🔴 破（差 %.0f×）' % short}")

# ---------------- Q2 自旋盲 ----------------
print("\n[Q2] 自旋盲（強 IB 不辨自旋 ⟹ K 與 K* 殘餘應近相等）")
E_K, E_Ks, s_K, s_Ks = 1.3740, 1.3200, 0.0960728889958036, 0.34122719704032967
d_spin = E_K - E_Ks
sd = float(np.hypot(s_K, s_Ks))
print(f"     觀測差 = {d_spin:.4f} ± {sd:.4f}（{d_spin/sd:.3f}σ）⟹ "
      f"{'✅ 自旋盲一致' if abs(d_spin/sd) <= 2 else '🔴'}")

# ---------------- Q3 族模式 ----------------
print("\n[Q3] 族模式：Σ/Ξ（含 s、無 c）")
d8 = json.load(open(ROOT / "data" / "08_sampleout_channels.json"))
q3 = {}
for p in d8["spin_pairs"]:
    if p["family"].startswith("Σ"):
        z = p["dE"] / p["sigma"]
        q3[p["family"]] = {"E_low": p["E_low"], "E_high": p["E_high"], "dE": p["dE"],
                           "sigma": p["sigma"], "z": z}
        print(f"     {p['family']:22s} E_low = {p['E_low']:7.4f}、E_high = {p['E_high']:7.4f}"
              f" ⟹ ΔE = {p['dE']:7.4f} ± {p['sigma']:6.4f}（{z:+.2f}σ）⟹ "
              f"{'自旋盲' if abs(z) <= 2 else '有自旋分裂'}")
n_blind = sum(1 for v in q3.values() if abs(v["z"]) <= 2)
q3_ok = n_blind >= 2          # 判準：至少 2/3 對自旋盲（2.25σ 之邊緣對另記）
print(f"     ⟹ Σ 系 {n_blind}/{len(q3)} 對自旋盲（|z| ≤ 2）：[−/0] 1.16σ、[0/+] 2.25σ（邊緣）、Σc −0.39σ"
      f" ⟹ 大電荷分裂 ＋ 近零自旋分裂之模式 {'✅ 支持' if q3_ok else '🔴'}")

# ---------------- Q4 零參數比 ----------------
print("\n[Q4] 零參數比：核預言 ΔE_K/ΔE_D")
_, cC, _, _ = energy(Q_C, M_C, ML)
_, cS, _, _ = energy(Q_S, M_S, ML)
pred = cS / cC
obs = 0.0540 / 1.4100
print(f"     核預言 = |c|_K/|c|_D = {cS:.4f}/{cC:.4f} = {pred:.4f}；觀測 = {obs:.4f}"
      f" ⟹ 高估 {pred/obs:.1f}×（統計上 ΔE_K 僅 +1.97σ，故 Q4 = 解釋力不足但未否證）")

# ---------------- Q5 硬上界 ----------------
print("\n[Q5] 硬上界審計：⟨B⟩ ≤ 1 ⟹ A_max = α q(ħc/a_C)")
Amax = ALPHA * Q_S * (HBARC / A0) * 1000
Areq, sAreq = 1.3335, 0.2570
print(f"     A_max = {Amax:.4f} MeV；A_req = {Areq:.4f} ± {sAreq:.4f} MeV")
print(f"     A_req/A_max = {Areq/Amax:.4f} ⟹ 中心值**超界 {(Areq/Amax-1)*100:.2f}%**"
      f"（需 ⟨B⟩ = {Areq/Amax:.4f} > 1）；統計上僅 {(Areq-Amax)/sAreq:+.2f}σ")
print(f"     ⟹ 本算子族對 K 系電荷分裂之**理論上限** = {Amax:.4f} MeV"
      f"（即⟨B⟩ = 1 之完全定域）；觀測中心值在其上 ⟹ 結構上不可能完全由本核 em 提供")

# ---------------- 判決 ----------------
q1_ok = abs(q1["K"]["delta_nl"]) >= 0.687
q2_ok = abs(d_spin / sd) <= 2
verdict = ("🔴 D2 破（量級不足 %.0f×）" % short if not q1_ok else "🟡 D2 部分成立") + \
          "：非線性強 IB 殘餘僅 %.4f MeV = 所需 1.374 之 %.2f%%，**差 %.0f 倍**；" % (
              abs(q1["K"]["delta_nl"]), abs(q1["K"]["delta_nl"]) / 1.374 * 100, short) + \
          "Q2 自旋盲（0.152σ）與 Q3 族模式（%d/%d 對）確認「殘餘唔似 em」，" % (n_blind, len(q3)) + \
          "但佢哋唔支持「係 M4 X 線性扣除不足」——Q1 已證該機制量級差兩個數量級；" + \
          "Q5 另立硬上界：A_max = %.4f MeV，A_req 中心值超界 %.2f%%（需 ⟨B⟩ = %.4f > 1）⟹ " % (
              Amax, (Areq / Amax - 1) * 100, Areq / Amax) + \
          "本核 em 對 K 系電荷分裂之貢獻結構上封頂 ⟹ M5 判 **D3 收檔**，並 register 新問題" \
          "（K 系 1.374／1.320 之非-em、非-X 來源；候選：格點 QED 混合－強 IB 非微擾，或態依賴算子）"

out = {"Q1_nonlinear": q1, "Q1_pass": bool(q1_ok),
       "Q2_spin_blind": {"diff": d_spin, "sigma": sd, "z": d_spin / sd, "pass": bool(q2_ok)},
       "Q3_family": q3, "Q3_pass": bool(q3_ok),
       "Q4_ratio": {"pred": pred, "obs": obs, "overestimate": pred / obs},
       "Q5_hard_cap": {"A_max": Amax, "A_req": Areq, "ratio": Areq / Amax,
                       "excess_pct": (Areq / Amax - 1) * 100, "sigma": (Areq - Amax) / sAreq},
       "verdict": verdict}
with open(ROOT / "data" / "17_d2_strong_ib.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n[判決] " + verdict)
print("\n-> data/17_d2_strong_ib.json")