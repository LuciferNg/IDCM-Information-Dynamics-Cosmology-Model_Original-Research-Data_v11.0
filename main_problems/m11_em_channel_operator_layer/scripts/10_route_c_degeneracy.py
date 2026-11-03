"""M11 scripts/10 — 路線 C 之失效解剖：**模型空間退化**（結構性，非擬合不足）
================================================================================================
`scripts/09`（路線 C，判準事前寫死）判 **KILL**：2 標度模型 χ²/dof = 3154、max|z| = 128。
本腳本回答「為何死」——分兩層：

  E1｜**退化（結構性，零擬合）**：2 向量模型 (Φ, Ψ) 之下，若干通道具**完全相同**之 (Φ, Ψ)
      但 Δm 相差達 2.04 MeV ⟹ 任何標度 (a, b) 都預言佢哋相等 ⟹ 無參數可救。
      呢個係模型空間維度不足之證明，與擬合好壞無關。

  E2｜**最小第三向量**（後驗 exploratory，明確標記）：由 E1 之退化類別反推，第三向量必須能分辨
      「隱藏重味數」。取 n_s（隱藏 s 數）作候選，3 向量 (Φ, Ψ, n_s) 擬合（3 參數／8 數據 ⟹ 5 dof）。
  E3｜同 E2，但加理論誤差地板 0.05 MeV（後驗；因 π 之 σ = 0.0005 MeV 會獨佔 χ²，屬實驗精度而非模型可及）。

**判詞規則（事前寫死）**：
  E1 若存在 (Φ, Ψ) 完全相同而 Δm 差 > 3σ_合併 之通道對 ⟹ 判「模型空間退化，2 向量不足」。
  E2／E3 只作 exploratory 指引（後驗），**不得**作為 M11 之 pass／kill 判詞。

輸出：`data/10_route_c_degeneracy.json`
"""
import itertools
import json
import os

import numpy as np

X, K0 = 2.56, 3.8152622441
M = {"p": (938.27208943, 2.9e-7), "n": (939.56542052, 5.4e-7),
     "K+": (493.677, 0.016), "K0": (497.611, 0.013),
     "pi+": (139.57039, 1.8e-4), "pi0": (134.9768, 5e-4),
     "S-": (1197.449, 0.030), "S0": (1192.642, 0.024), "S+": (1189.37, 0.07),
     "Xi-": (1321.71, 0.07), "Xi0": (1314.86, 0.20),
     "D+": (1869.66, 0.05), "D0": (1864.84, 0.05)}
QB = {"u": 2 / 3, "d": -1 / 3, "s": -1 / 3, "c": 2 / 3}
CH = [("n-p", "n", "p", 1, 0, 0), ("K+-K0", "K+", "K0", -1, 1, 0), ("pi+-pi0", "pi+", "pi0", 1, 0, 0),
      ("S--S0", "S-", "S0", 1, 1, 0), ("S0-S+", "S0", "S+", 1, 1, 0), ("S--S+", "S-", "S+", 2, 1, 0),
      ("Xi--Xi0", "Xi-", "Xi0", 1, 2, 0), ("D+-D0", "D+", "D0", 1, 0, 1)]
CONTENT = {
    "n-p": ([("d", 2), ("u", 1)], [("u", 2), ("d", 1)]),
    "K+-K0": ([("u", 1), ("sbar", 1)], [("d", 1), ("sbar", 1)]),
    "pi+-pi0": ([("u", 1), ("dbar", 1)], "PI0"),
    "S--S0": ([("d", 2), ("s", 1)], [("u", 1), ("d", 1), ("s", 1)]),
    "S0-S+": ([("u", 1), ("d", 1), ("s", 1)], [("u", 2), ("s", 1)]),
    "S--S+": ([("d", 2), ("s", 1)], [("u", 2), ("s", 1)]),
    "Xi--Xi0": ([("d", 1), ("s", 2)], [("u", 1), ("s", 2)]),
    "D+-D0": ([("c", 1), ("dbar", 1)], [("c", 1), ("ubar", 1)]),
}


def q2(spec):
    if spec == "PI0":
        a, b = q2([("u", 1), ("ubar", 1)]), q2([("d", 1), ("dbar", 1)])
        return (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
    qs = [QB[k.rstrip("bar")] * (-1 if k.endswith("bar") else 1) for k, n in spec for _ in range(n)]
    return (sum(q * q for q in qs),
            sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs))))


def fit(A, y, s):
    W = 1 / s
    p, *_ = np.linalg.lstsq(A * W[:, None], y * W, rcond=None)
    r = (y - A @ p) / s
    return p, float((r ** 2).sum()), r


def main():
    rows = []
    for name, hi, lo, coeff, ns, nc in CH:
        dm = M[hi][0] - M[lo][0]
        sd = (M[hi][1] ** 2 + M[lo][1] ** 2) ** 0.5
        a2, ap = q2(CONTENT[name][0])
        b2, bp = q2(CONTENT[name][1])
        rows.append(dict(ch=name, dm=dm, sd=sd, phi=((a2 - b2) + (ap - bp)) / 2, coeff=coeff, ns=ns, nc=nc))
    Y = np.array([r["dm"] for r in rows]); S = np.array([r["sd"] for r in rows])

    print("=" * 100)
    print("M11 scripts/10 — 路線 C 之失效解剖：模型空間退化")
    print("=" * 100)
    print("[E1] 同 (Φ, Ψ) 之通道類（零擬合；同類必須預言相等）")
    print(f"{'類別 (Φ, Ψ)':18s}{'通道':10s}{'Δm':>9s}{'σ':>8s}")
    bad = []
    for key, grp in itertools.groupby(sorted(rows, key=lambda r: (r["phi"], r["coeff"])),
                                      key=lambda r: (round(r["phi"], 6), r["coeff"])):
        g = list(grp)
        if len(g) < 2:
            continue
        tag = f"({key[0]:.4f}, {key[1]:d})"
        for r in g:
            print(f"{tag:18s}{r['ch']:10s}{r['dm']:9.4f}{r['sd']:8.4f}")
        spread = max(x["dm"] for x in g) - min(x["dm"] for x in g)
        sig = (sum(x["sd"] ** 2 for x in g)) ** 0.5
        z = spread / max(sig, 1e-12)
        print(f"{'':18s}⟹ 類內差 {spread:.4f} MeV、合併 σ = {sig:.4f} ⟹ **{z:.1f}σ** ⟹ "
              f"{'退化（違反）' if z > 3 else 'OK'}")
        if z > 3:
            bad.append(dict(cls=tag, channels=[x["ch"] for x in g], spread=spread, z=float(z)))
    E1 = dict(degenerate_classes=bad, verdict="模型空間退化：2 向量 (Φ,Ψ) 不足" if bad else "OK")
    print(f"\n  ⟹ 判：**{E1['verdict']}**（退化類 {len(bad)} 個）")
    print("  ⟹ 結論：任何 (a,b) 都預言同類相等；類內差係**無參數可救**之結構缺陷 ⟹ 必須引入第三向量。")

    print()
    print("[E2] 後驗 exploratory：3 向量 (Φ, Ψ, n_s) —— 3 參數／8 數據 ⟹ 5 dof")
    ns = np.array([r["ns"] for r in rows])
    A3 = np.array([np.array([r["phi"] for r in rows]), np.array([r["coeff"] for r in rows]), ns]).T
    p3, c3, r3 = fit(A3, Y, S)
    print(f"  a = {p3[0]:.4f}、b = {p3[1]:.4f}（⟨O⟩_eff = {p3[1]/X:.3f}）、c_s = {p3[2]:.4f} MeV／s")
    print(f"  χ² = {c3:.1f}、χ²/dof = {c3/5:.1f}、max|z| = {max(abs(r3)):.1f}")
    for i, r in enumerate(rows):
        print(f"    {r['ch']:9s} obs={r['dm']:8.4f} pred={Y[i]-r3[i]*S[i]:8.4f} z={r3[i]:+8.1f}")

    print()
    print("[E3] 同 E2，加理論誤差地板 0.05 MeV（後驗；π 之 σ=0.0005 屬實驗精度）")
    S2 = (S ** 2 + 0.05 ** 2) ** 0.5
    p3b, c3b, r3b = fit(A3, Y, S2)
    print(f"  a = {p3b[0]:.4f}、b = {p3b[1]:.4f}（⟨O⟩_eff = {p3b[1]/X:.3f}）、c_s = {p3b[2]:.4f} MeV／s")
    print(f"  χ² = {c3b:.1f}、χ²/dof = {c3b/5:.2f}、max|z| = {max(abs(r3b)):.2f}")
    for i, r in enumerate(rows):
        print(f"    {r['ch']:9s} z={r3b[i]:+6.2f}")

    print()
    print("[E4] 4 向量 (Φ, Ψ, n_s, n_c) —— 4 參數／8 數據 ⟹ 4 dof（後驗）")
    nc = np.array([r["nc"] for r in rows])
    A4 = np.hstack([A3, nc[:, None]])
    p4, c4, r4 = fit(A4, Y, S2)
    print(f"  a = {p4[0]:.4f}、b = {p4[1]:.4f}（⟨O⟩_eff = {p4[1]/X:.3f}）、c_s = {p4[2]:.4f}、c_c = {p4[3]:.4f}")
    print(f"  χ² = {c4:.2f}、χ²/dof = {c4/4:.2f}、max|z| = {max(abs(r4)):.2f}")

    _out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "10_route_c_degeneracy.json")
    with open(_out, "w") as f:
        json.dump(dict(rows=rows, E1=E1, E2=dict(a=float(p3[0]), b=float(p3[1]), c_s=float(p3[2]),
                                                  chi2=c3, chi2_dof=c3 / 5, max_z=float(max(abs(r3)))),
                       E3=dict(a=float(p3b[0]), b=float(p3b[1]), c_s=float(p3b[2]),
                               chi2=c3b, chi2_dof=c3b / 5, max_z=float(max(abs(r3b)))),
                       E4=dict(a=float(p4[0]), b=float(p4[1]), c_s=float(p4[2]), c_c=float(p4[3]),
                               chi2=c4, chi2_dof=c4 / 4, max_z=float(max(abs(r4))))),
                  f, ensure_ascii=False, indent=2)
    print("\n-> data/10_route_c_degeneracy.json")


main()