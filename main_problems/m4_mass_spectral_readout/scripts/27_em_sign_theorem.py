"""M4 docs/27 配套 — em 符號定理（零參數、與波函數無關）：8 通道逐通道判定
=========================================================================================
定理（不假設任何標度、任何波函數形狀）：
  對任意正核（`V(d) > 0`，如 Coulomb `B^{−d}` 或 on-site 接觸）與任意歸一化態，
  `em(A) − em(B) = Δ(Σq²)·S + Δ(Σ_{i<j} q_iq_j)·P`，其中 `S, P ≥ 0`。
  ⟹ 若兩個 Δ 之**符號相同**，則 `sign(em(A) − em(B))` **被固定**，與態、標度、耦合強度全部無關。
  ⟹ 若符號相反，則無符號預言（本輪記「無判別力」）。

判據（事前寫死）：
  S1 用 Fraction 精確計 `ΔΣq²`、`ΔΣ_{i<j}q_iq_j`；只對同號通道出預言。
  S2 預言符號 vs 實測殘餘 `E`（`docs/18`：`E = Δm − coeff·X`，`X = 2.56`）之符號；
     `|E| ≥ 3σ` 才算「有判別力之檢定」。
  S3 通道唔准改；夸克電荷用 Fraction（`u = +2/3`、`d = s = −1/3`、`c = +2/3`；反夸克取負）。
  S4 唔准為救任何通道而改 `X` 或改 `E` 定義。
輸出：`data/27_em_sign.json`
"""
import json
from fractions import Fraction as F

QU, QD, QS, QC = F(2, 3), F(-1, 3), F(-1, 3), F(2, 3)
QMAP = {"u": QU, "d": QD, "s": QS, "c": QC}
X_STRONG, X_SIGMA = 2.56, 0.094

# (名稱, 第一態價夸克, 第二態價夸克, Δm = m(第一) − m(第二), coeff, E, σ_E)  ← docs/18
CH = [
    ("n−p", "udd", "uud", +1.2933, +1, -1.2667, 0.0940),
    ("K⁺−K⁰", "u sbar", "d sbar", -3.9340, -1, -1.3740, 0.0963),
    ("Σ⁻−Σ⁰", "dds", "uds", +4.8070, +1, +2.2470, 0.1010),
    ("Σ⁰−Σ⁺", "uds", "uus", +3.2720, +1, +0.7120, 0.1196),
    ("Σ⁻−Σ⁺", "dds", "uus", +8.0790, +2, +2.9590, 0.2028),
    ("Ξ⁻−Ξ⁰", "dss", "uss", +6.8500, +1, +4.2900, 0.2301),
    ("D⁺−D⁰", "c dbar", "c ubar", +4.8200, +1, +2.2600, 0.1172),
]


def parse(spec):
    toks = []
    for tok in spec.split():
        if tok in QMAP:
            toks.append(tok)
        elif tok.endswith("bar"):
            toks.append("-" + tok[:-3])
        else:
            toks.extend(list(tok))
    out = []
    for t in toks:
        neg = t.startswith("-")
        base = t.lstrip("-")
        q = QMAP[base]
        out.append(-q if neg else q)
    return out


def algebra(qs):
    return (sum(q * q for q in qs),
            sum(qs[i] * qs[j] for i in range(len(qs)) for j in range(i + 1, len(qs))))


print("=" * 104)
print("em 符號定理：對任意正核、任意波函數，同號通道之 sign(em(A)−em(B)) 被固定")
print(f"（殘餘 E ≡ em(A) − em(B)，A = 通道名第一態；X = {X_STRONG} ± {X_SIGMA}；預言方向 = −sign(Δs2+Δpair)）")
print("=" * 104)
print(f"{'通道':<10}{'ΔΣq²':>9}{'ΔΣpair':>10}{'同號?':>7}{'預言':>7}{'實測 E ± σ':>18}{'實測符':>7}{'|E|/σ':>7}{'判':>8}")
rows, verdicts = [], []
for name, a, b, dm, coeff, E, sE in CH:
    qa, qb = parse(a), parse(b)
    s2a, pra = algebra(qa)
    s2b, prb = algebra(qb)
    ds2, dpr = s2b - s2a, prb - pra          # 係 em(B) − em(A) 之代數係數
    same = (ds2 > 0 and dpr > 0) or (ds2 < 0 and dpr < 0)
    # 方向校正（必守）：em(A) − em(B) = −[Δ s2·S + Δ pair·P]，A = 通道名第一態、B = 第二態。
    # E 之定義同為 em(A) − em(B)（由 n−p 核對：Δm=+1.2933、coeff=+1 ⟹ E=+1.2933−2.56=−1.2667 = em(n)−em(p)）。
    pred = (-1 if ds2 + dpr > 0 else 1) if same else 0     # 預言 sign(em(A) − em(B))
    obs = 1 if E > 0 else -1                                 # 實測 sign(m(B) − m(A)) 之殘餘
    z = abs(E) / sE
    if not same:
        judge = "無判別"
        verdicts.append((name, "ambiguous", None, obs, z))
    elif obs == pred and z >= 3:
        judge = "PASS"
        verdicts.append((name, "pass", pred, obs, z))
    elif obs != pred and z >= 3:
        judge = "FAIL"
        verdicts.append((name, "fail", pred, obs, z))
    else:
        judge = "弱"
        verdicts.append((name, "weak", pred, obs, z))
    rows.append({"channel": name, "dSigmaQ2": str(ds2), "dSigmaPair": str(dpr), "same_sign": same,
                 "predicted_sign": pred, "E": E, "sigma": sE, "observed_sign": obs,
                 "z": z, "verdict": judge})
    print(f"{name:<10}{str(ds2):>9}{str(dpr):>10}{('是' if same else '否'):>7}"
          f"{({1:'+',-1:'−',0:'—'}[pred]):>7}{E:>12.4f} ± {sE:<5.4}{({1:'+',-1:'−'}[obs]):>7}"
          f"{z:>7.1f}{judge:>8}")

tested = [v for v in verdicts if v[1] in ("pass", "fail")]
print("-" * 104)
print(f"有判別力之通道：{len(tested)}（PASS {sum(1 for v in tested if v[1]=='pass')}／FAIL {sum(1 for v in tested if v[1]=='fail')}）；"
      f"其餘 {len(verdicts)-len(tested)} 條因兩項反號而無符號預言")

print("\n【定理之邏輯內容】")
print("  1. 同號通道之預言**不依賴**耦合強度、節點間距 a、波函數形狀、動態 P₀ —— 只要核為正、態歸一。")
print("  2. 故 FAIL 之通道**不可能**由任何 em 正核解釋（含 Coulomb、contact、及其任意混合）。")
print("  3. PASS 之通道只證「符號相容」，不證幅度（幅度見 docs/26 之承載窗）。")

print("\n【X 不等式化】E(X) = Δm − coeff·X；要求 sign(E(X)) = 預言 ⟹ 對 X 之單邊界")
lo, hi, cons = [], [], []
for (name, a, b, dm, coeff, E, sE), v in zip(CH, verdicts):
    if v[1] not in ("pass", "fail"):
        continue
    pred = v[2]
    thr = dm / coeff if coeff > 0 else -dm / abs(coeff)
    if (coeff > 0 and pred < 0) or (coeff < 0 and pred > 0):
        lo.append((name, thr))
        cons.append(f"{name}: X > {thr:.4f}")
    else:
        hi.append((name, thr))
        cons.append(f"{name}: X < {thr:.4f}")
print("  " + " ；".join(cons))
Xmin = max([t for _, t in lo]) if lo else float("-inf")
Xmax = min([t for _, t in hi]) if hi else float("inf")
print(f"  ⟹ 三通道聯合：{Xmin:.4f} < X < {Xmax if Xmax != float('inf') else '∞'}")
print(f"  框架值 X = {X_STRONG} ± {X_SIGMA} ⟹ 距下界 {Xmin:.4f} 有 {(Xmin - X_STRONG)/X_SIGMA:.1f}σ"
      f" ⟹ {'🔴 不一致' if Xmin > X_STRONG else '✅ 相容'}")
print("  註：X 之物理身份 = m_d − m_u 標度（FLAG24 Δm_ud = 2.56 ± 0.094，見 docs/18）。")

print("\n【物理註】")
print("  K⁰ 比 K⁺ 重 3.934 MeV、Σ⁰ 比 Σ⁺ 重 3.272 MeV —— 兩者皆為文獻上『帶電態反而較輕』之異常；")
print("  本輪證：即使扣除 u↔d 質量差（X = 2.56），該異常喺殘餘層仍然存在（14σ / 5.9σ），")
print("  而 em（正核）之符號預言與之相**反** ⟹ 殘餘不可能係純 em。")

json.dump({"X_strong": X_STRONG, "rows": rows,
           "X_inequalities": {"lower": lo, "upper": hi,
                              "Xmin": Xmin, "Xmax": None if Xmax == float("inf") else Xmax,
                              "sigma_distance": (Xmin - X_STRONG) / X_SIGMA},
           "summary": {"tested": len(tested),
                       "pass": sum(1 for v in tested if v[1] == "pass"),
                       "fail": sum(1 for v in tested if v[1] == "fail"),
                       "ambiguous": len(verdicts) - len(tested)}},
          open("data/27_em_sign.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/27_em_sign.json")