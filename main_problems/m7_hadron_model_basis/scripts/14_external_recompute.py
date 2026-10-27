"""M7 scripts/14 — 外部結果獨立重算與判決：em(K⁺−K⁰) 之 lattice QED 值
==========================================================================================================
背景：M4 `docs/36`（R5）判本線 em 通道 KILL（−15.98σ）。M7 `docs/21` 只**轉錄**、未重算 ⟹ 本檔補此缺口。
依 skill `cross-record-result-import-verification`：§1 由 primary source 重建、§2 約定、§3 判詞。

Primary source（已取原文，非二手）：
  · FLAG Review 2024 (arXiv:2411.04268v3) §4.1.3 Eq. (38)：
        ϵ = (ΔM²_K − ΔM²_π)_γ / ΔM²_π ,  ΔM²_π = M²_{π⁺} − M²_{π⁰} ,  ΔM²_K = M²_{K⁺} − M²_{K⁰}
    同節明文：G-parity ⟹ ΔM²_π 無 O(m_u−m_d) 同位旋破缺修正 ⟹ (ΔM²_π)_γ = ΔM²_π
  · MILC 18 (PRD 99, 034503 (2019)；arXiv:1807.05556) Eq. (1)（**同一式，原文**）：
        ϵ ≡ [(M²_{K±} − M²_{K⁰})_γ − (M²_{π±} − M²_{π⁰})_γ] / (M²_{π±} − M²_{π⁰})_expt
    值：ϵ = 0.78(1)_stat (+8/−11)_syst ⟹ 帶 [0.67, 0.86]
  · PDG 2024 質量（MeV）：m_{K±} 493.677±0.016；m_{K⁰} 497.611±0.013；
        m_{π±} 139.57039±0.00018；m_{π⁰} 134.9768±0.0005
被檢預言（M4 `docs/29` §三，2026-10-25 寫死；取向 M(K⁺)−M(K⁰)）：em = +0.3327 MeV
輸出：data/14_external_recompute.json
"""
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "14_external_recompute.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<48} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 104 + f"\n{s}\n" + "=" * 104)


def dec(s):                                    # 字串 → 精確 Fraction
    return F(s)


def show(x, n=8):
    return nstr(mpf(x.numerator) / mpf(x.denominator), n)


# ── 輸入（純精確有理數）────────────────────────────────────────────
mKc, sKc = dec("493.677"), dec("0.016")
mK0, sK0 = dec("497.611"), dec("0.013")
mPic, sPic = dec("139.57039"), dec("0.00018")
mPi0, sPi0 = dec("134.9768"), dec("0.0005")
EPS, EPS_LO, EPS_HI = dec("0.78"), dec("0.67"), dec("0.86")
PRED = dec("0.3327")                           # M4 docs/29 §三（取向 K⁺−K⁰）
SUM_M = mKc + mK0

print("=" * 104)
print("M7 scripts/14 — lattice QED em(K⁺−K⁰) 獨立重算（FLAG Eq.38 ＋ MILC Eq.1 ＋ PDG）")
print("=" * 104)

# ══════════ §1 平方質量分裂 ══════════
hdr("§1 由 PDG 質量重算 ΔM²_π、ΔM²_K（精確有理數）")
dM2_pi = mPic ** 2 - mPi0 ** 2
dM2_K = mKc ** 2 - mK0 ** 2
print(f"    ΔM²_π = {show(dM2_pi, 10)} MeV²   （M4 表：1261.157）")
print(f"    ΔM²_K = {show(dM2_K, 10)} MeV²   （M4 表：−3899.73，負 ⟹ K⁰ 較重：PDG 497.611 vs 493.677）")
ok("ΔM²_π 重算 == 1261.157（M4 表值）", abs(mpf(dM2_pi.numerator) / dM2_pi.denominator - mpf("1261.157")) < mpf("0.001"),
   show(dM2_pi, 10), "1261.157")
ok("ΔM²_K 重算 == −3899.73（M4 表值）", abs(mpf(dM2_K.numerator) / dM2_K.denominator + mpf("3899.73")) < mpf("0.01"),
   show(dM2_K, 10), "-3899.73")

# ══════════ §2 ε 定義 → (ΔM²_K)_γ ══════════
hdr("§2 由 FLAG Eq.(38) ＋ G-parity 推 (ΔM²_K)_γ = (1+ϵ)ΔM²_π（代數，非引用結論）")
print("    ϵ = [(ΔM²_K)_γ − (ΔM²_π)_γ]/ΔM²_π ；G-parity ⟹ (ΔM²_π)_γ = ΔM²_π")
print("    ⟹ ϵ = (ΔM²_K)_γ/ΔM²_π − 1 ⟹ (ΔM²_K)_γ = (1+ϵ)·ΔM²_π")
dM2_K_gamma = (1 + EPS) * dM2_pi
dM2_K_gamma_lo = (1 + EPS_LO) * dM2_pi
dM2_K_gamma_hi = (1 + EPS_HI) * dM2_pi
em = dM2_K_gamma / SUM_M
em_lo = dM2_K_gamma_lo / SUM_M
em_hi = dM2_K_gamma_hi / SUM_M
sigma = (em_hi - em_lo) / 2
print(f"    ΣM = m_K⁺+m_K⁰ = {show(SUM_M, 10)} MeV")
print(f"    em(K⁺−K⁰) = (ΔM²_K)_γ/ΣM = {show(em, 8)} MeV；68% 帶 [{show(em_lo, 7)}, {show(em_hi, 7)}]；σ = {show(sigma, 6)}")
ok("em 重算 == +2.2646（M4 值）", abs(mpf(em.numerator) / em.denominator - mpf("2.2646")) < mpf("5e-5"),
   show(em, 8), "2.2646")
ok("68% 帶 == [2.1246, 2.3664]（M4 值）",
   abs(mpf(em_lo.numerator) / em_lo.denominator - mpf("2.1246")) < mpf("5e-4")
   and abs(mpf(em_hi.numerator) / em_hi.denominator - mpf("2.3664")) < mpf("5e-4"),
   f"[{show(em_lo,7)}, {show(em_hi,7)}]", "[2.1246, 2.3664]")
ok("σ == 0.1209（M4 用值）", abs(mpf(sigma.numerator) / sigma.denominator - mpf("0.1209")) < mpf("5e-5"),
   show(sigma, 6), "0.1209")

# ══════════ §3 判決 ══════════
hdr("§3 對 M4 寫死預言 em = +0.3327 之判決")
z = (PRED - em) / sigma
z_band = (PRED - em_lo) / sigma                     # 帶內最小 |z|
print(f"    z = (0.3327 − {show(em, 7)})/{show(sigma, 6)} = {float(z):+.2f}")
print(f"    帶內最保守（ϵ=0.67 ⟹ em 最小）= {float(z_band):+.2f} ⟹ 整個 ϵ 帶皆 KILL")
ok("判決 = KILL（|z| ≥ 3），且 ϵ 帶內全 KILL", abs(z) >= 3 and abs(z_band) >= 3,
   f"z = {float(z):+.2f}；帶內 {float(z_band):+.2f}", "|z|≥3",
   "本 M 獨立重算，非轉錄")
ok("判決對 M4 之 z 一致（−15.98）", abs(float(z) + 15.98) < 0.05, f"{float(z):+.2f}", "-15.98",
   "同一 σ 約定（帶半寬）下逐位一致")

# ══════════ §4 新診斷：預言與「(1−ϵ) 符號錯」退化 ══════════
hdr("§4 新診斷（本 M 提出）：預言 +0.3327 與『(1−ϵ) 讀法』數值退化")
em_wrong = (1 - EPS) * dM2_pi / SUM_M
print(f"    若誤讀為 (1−ϵ)：em_wrong = {show(em_wrong, 8)} MeV")
print(f"    預言 0.3327 距 em_wrong = {float((PRED-em_wrong)/sigma):+.2f}σ ；距正確值 = {float(z):+.2f}σ")
ok("預言落在『(1−ϵ) 錯誤讀法』之 1σ 內 ⟹ 預言 ≡ 符號錯之產物",
   abs((PRED - em_wrong) / sigma) < 1, f"{float((PRED-em_wrong)/sigma):+.2f}σ", "<1σ",
   "⟹ 本線 em 通道之『命中』只係複製咗 M4 v1 之符號錯，非結構巧合")

# ══════════ §5 反推：外部值要求之 ⟨O⟩_K ══════════
hdr("§5 反推：外部值要求之 ⟨O⟩_K，與 M4 docs/29 之 ⟨O⟩ = 1+(B−2)n_s/3 比較")
DM_K, X_M7 = dec("-3.9340"), dec("2.56")            # Δm（M4 docs/18 取向）、X（M7 產權）
O_pred = (PRED - DM_K) / X_M7
O_req = (em - DM_K) / X_M7
O_sig = dec("0.262")                                # BMW Σ 錨 1.498±0.262
print(f"    預言用 ⟨O⟩_K = (0.3327 − (−3.9340))/2.56 = {show(O_pred, 6)}（M4 記 1.667）")
print(f"    外部值要求 ⟨O⟩_K = ({show(em,6)} + 3.9340)/2.56 = {show(O_req, 6)}")
print(f"    與 Σ 錨 1.498±0.262 之距離 = {float((O_req - dec('1.498'))/O_sig):+.2f}σ")
ok("預言用 ⟨O⟩_K 重算 == 1.667", abs(mpf(O_pred.numerator) / O_pred.denominator - mpf("1.667")) < mpf("5e-4"),
   show(O_pred, 6), "1.667")
ok("外部值要求 ⟨O⟩_K 遠離 Σ 錨（>3σ）⟹ ⟨O⟩=a+b·n_s 線亦被 K 通道排除",
   abs((O_req - dec("1.498")) / O_sig) > 3, f"{float((O_req-dec('1.498'))/O_sig):+.2f}σ", ">3σ",
   "M4 docs/29 C2（χ²=0.09/1 相容三錨）在 K 通道失效")

# ══════════ §6 穩健性 ══════════
hdr("§6 穩健性：三個假設逐個鬆動，判決是否改變")
alt = {}
for lbl, gp, fac in [("(ΔM²_π)_γ = 0.9ΔM²_π", dec("0.9"), None),
                     ("(ΔM²_π)_γ = 0.8ΔM²_π", dec("0.8"), None),
                     ("σ = 0.05（只用 ϵ 統計+對稱化系統）", None, dec("0.0634"))]:
    if gp is not None:
        v = (EPS * dM2_pi + gp * dM2_pi) / SUM_M
        zz = (PRED - v) / sigma
    else:
        v, zz = em, (PRED - em) / fac
    alt[lbl] = (float(v), float(zz))
    print(f"    {lbl:<34} ⟹ em = {float(v):.4f} MeV，z = {float(zz):+.2f}")
ok("三種鬆動下判決皆 KILL（|z| ≥ 3）", all(abs(zz) >= 3 for _, zz in alt.values()),
   ", ".join(f"{float(zz):+.1f}" for _, zz in alt.values()), "全部 |z|≥3",
   "⟹ 判決唔依賴單一約定選擇")

print("\n" + "=" * 104)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(dM2_pi=float(mpf(dM2_pi.numerator) / dM2_pi.denominator),
                        dM2_K=float(mpf(dM2_K.numerator) / dM2_K.denominator),
                        sum_M=float(mpf(SUM_M.numerator) / SUM_M.denominator),
                        em=float(mpf(em.numerator) / em.denominator),
                        em_lo=float(mpf(em_lo.numerator) / em_lo.denominator),
                        em_hi=float(mpf(em_hi.numerator) / em_hi.denominator),
                        sigma=float(mpf(sigma.numerator) / sigma.denominator),
                        z=float(z), z_band=float(z_band), pred=0.3327,
                        em_wrong_1meps=float(mpf(em_wrong.numerator) / em_wrong.denominator),
                        O_pred=float(mpf(O_pred.numerator) / O_pred.denominator),
                        O_req=float(mpf(O_req.numerator) / O_req.denominator),
                        O_req_sigma_off=float((O_req - dec("1.498")) / O_sig)),
           verdict="KILL (independent recompute)")
json.dump(OUT, open(os.path.join(M7, "data", "14_external_recompute.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/14_external_recompute.json")