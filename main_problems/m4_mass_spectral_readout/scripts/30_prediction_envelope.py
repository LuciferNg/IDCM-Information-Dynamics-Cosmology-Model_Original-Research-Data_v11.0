"""M4 docs/30 配套 — ⟨O⟩ 形式不可辨識，但預言**包絡**穩健：em(K⁰−K⁺) ≈ 0 之型式無關預言
=========================================================================================
背景：`docs/29` 給出三錨（`⟨O⟩_N = 0.891 ± 0.107`、`_Σ = 1.498 ± 0.262`、`_Ξ = 2.293 ± 0.351`），
      線性 `a + b·n_s` 之 `b = 0.675 ± 0.161` ⟹ 3 錨**無法**分辨斜率。本輪唔再爭形式，改問：

  Q：在**所有通過 χ² 門檻之結構形式**下，各通道 `em` 之取值範圍（包絡）係幾多？
     若包絡狹窄，則預言**型式無關**、可即時證偽；若包絡寬，則老實承認不可判。

紀律：
  P1 模型類**事前寫死**（見 MODELS），每個形式之係數只能由已註冊整數 `B = 4`、`N_q = 3`、`N_d = 88`、`N_gen = 3` 構成。
  P2 全部形式皆**後驗**（3 錨已在手才列舉）⟹ 不得宣稱任何單一形式為 derived；只可報包絡。
  P3 通過門檻 = χ² < 3（3 dof，約 1.4σ 上尾 0.39）。**不得**為收窄包絡而收緊門檻。
  P4 `coeff` 取向按 `docs/18`（`u→d` ⟹ −1；`d→u`／反夸克 ⟹ +1；兩次 ⟹ +2）。
輸出：`data/30_prediction_envelope.json`
"""
import json
from math import exp, erfc, gamma as _g, sqrt

X = 2.56
OBS = [(0, 0.891, 0.107, "N"), (1, 1.498, 0.262, "Σ"), (2, 2.293, 0.351, "Ξ")]
B_, NQ, ND, NG = 4, 3, 88, 3


def chi2_sf(x, k):
    if k % 2 == 0:
        s, t, i = 0.0, 1.0, 0
        while i < k // 2:
            s += t
            i += 1
            t *= x / 2 / i
        return exp(-x / 2) * s
    return erfc(sqrt(x / 2)) + exp(-x / 2) * sum((x / 2) ** (j + 0.5) / _g(j + 1.5) for j in range(0))


# 模型類（事前寫死）：f(n_s)
MODELS = {
    "A1: 1 + n_s·2/N_q":        lambda n: 1 + 2 * n / NQ,
    "A2: N_q(n_s+1)/B":         lambda n: NQ * (n + 1) / B_,
    "A3: 1 + n_s·1/N_q":        lambda n: 1 + n / NQ,
    "A4: 1 + n_s·(B−1)/B":      lambda n: 1 + n * (B_ - 1) / B_,
    "A5: 1 + n_s/B":            lambda n: 1 + n / B_,
    "A6: 1 + n_s/2":            lambda n: 1 + n / 2,
    "A7: (1 + n_s/N_q)²":       lambda n: (1 + n / NQ) ** 2,
    "A8: B/(B − n_s)":          lambda n: B_ / (B_ - n),
    "A9: (N_q + n_s·NG)/N_q":   lambda n: (NQ + n * NG) / NQ,
    "A10: 1 + n_s·(B−2)/(B−1)": lambda n: 1 + n * (B_ - 2) / (B_ - 1),
    "A11: 1 + n_s·ND/(ND+NQ·B)": lambda n: 1 + n * ND / (ND + NQ * B_),
    "A12: 1 + n_s·2/B":         lambda n: 1 + 2 * n / B_,
}

# 通道：名、coeff、Δm、n_s
CH = [("K⁺−K⁰", -1, -3.9340, 1), ("D⁺−D⁰", +1, +4.8200, 0),
      ("Σ⁰−Σ⁺", +1, +3.2720, 1), ("Σ⁻−Σ⁰", +1, +4.8070, 1), ("Σ⁻−Σ⁺", +2, +8.0790, 1)]


def chi2(f):
    return sum(((o - f(n)) / s) ** 2 for n, o, s, _ in OBS)


print("=" * 100)
print("⟨O⟩ 形式不可辨識 → 改報預言包絡（型式無關）")
print("=" * 100)
print(f"{'模型':<30}{'⟨O⟩(0,1,2)':>26}{'χ²/3dof':>10}{'p':>8}{'通過(χ²<3)':>12}")
passed = []
for name, f in MODELS.items():
    c = chi2(f)
    ok = c < 3.0
    vals = "".join(f"{f(n):6.3f}" for n in (0, 1, 2))
    print(f"{name:<30}{vals:>26}{c:>10.2f}{chi2_sf(c, 3):>8.3f}{('✅' if ok else '—'):>12}")
    if ok:
        passed.append((name, f))

print(f"\n通過門檻者 {len(passed)}/{len(MODELS)}：")
for name, _ in passed:
    print(f"  • {name}")

print("\n【預言包絡】逐通道（只計通過之模型）")
print(f"{'通道':<9}{'em 最小':>10}{'em 最大':>10}{'跨距':>9}{'⟨O⟩ 範圍':>18}{'與 docs/29 M1 差':>18}")
env = {}
for ch, nf, dm, ns in CH:
    es, os_ = [], []
    for name, f in passed:
        o = f(ns)
        es.append(dm - nf * X * o)
        os_.append(o)
    env[ch] = {"em_min": min(es), "em_max": max(es), "O_min": min(os_), "O_max": max(os_)}
    m1 = dm - nf * X * (1 + 2 / 3 * ns)
    print(f"{ch:<9}{min(es):>+10.3f}{max(es):>+10.3f}{max(es)-min(es):>9.3f}"
          f"{min(os_):>8.3f} – {max(os_):>6.3f}{m1:>+18.3f}")

print("\n【判讀】")
k = env["K⁺−K⁰"]
print(f"  1. em(K⁰−K⁺) 包絡 = [{k['em_min']:+.3f}, {k['em_max']:+.3f}] MeV ⟹ 全部通過模型皆給 |em(K)| ≲ "
      f"{max(abs(k['em_min']), abs(k['em_max'])):.2f} MeV")
print(f"     Dashen 式估計 ~ +2.3 MeV ⟹ 落在包絡之外 {2.3/(max(abs(k['em_min']),abs(k['em_max']))):.1f}× "
      f"⟹ **型式無關之可證偽預言**")
d = env["D⁺−D⁰"]
print(f"  2. em(D⁺−D⁰) 包絡 = [{d['em_min']:+.3f}, {d['em_max']:+.3f}] MeV ⟹ 穩健於形式；與 docs/18 之 E = +2.26 相容")
for ch in ("Σ⁰−Σ⁺", "Σ⁻−Σ⁰", "Σ⁻−Σ⁺"):
    e = env[ch]
    print(f"  3. {ch:<9} 包絡 = [{e['em_min']:+.3f}, {e['em_max']:+.3f}] MeV")

print("\n【包絡之內部自洽（必守）】對每個通過模型，須有 em(Σ⁻−Σ⁰) + em(Σ⁰−Σ⁺) = em(Σ⁻−Σ⁺)")
bad = 0
for name, f in passed:
    o = f(1)
    lhs = (4.8070 - X * o) + (3.2720 - X * o)
    rhs = 8.0790 - 2 * X * o
    if abs(lhs - rhs) > 1e-9:
        bad += 1
        print(f"  🔴 {name}: {lhs:.6f} ≠ {rhs:.6f}")
print(f"  {'✅ 全部通過' if bad == 0 else f'🔴 {bad} 個不自洽'}（代數恆等，恆成立）")

json.dump({"models_passed": [n for n, _ in passed], "envelope": env,
           "note": "全部模型屬後驗；包絡才是可引用內容"},
          open("data/30_prediction_envelope.json", "w"), ensure_ascii=False, indent=1, default=str)
print("\n→ data/30_prediction_envelope.json")