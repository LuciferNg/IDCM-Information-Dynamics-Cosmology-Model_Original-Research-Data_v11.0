"""M5 docs/10 配套（第二段）— A2b：以**結構標度**固定節點間距 a，再由 A 定耦合，盲測 P₀
================================================================================
A2a（scripts/10）已證：(A_D, K_D) 兩條方程對 (a, t) **秩虧**（a→0, t→∞ 整條射線都符合）
  ⟹ 絕對標度**不可由數據定**，必須由結構供。
A2b：取結構標度候選（全部已登記，非擬合），固定 a 後只剩 t 一個未知：
  ① a = ħc/Λ_IDCM（O142 v16：Λ_IDCM = 45.2849 MeV，IDCM 自己導出之 QCD 標度）
  ② a = 1/(2μ)，μ = m_c m_l/(m_c+m_l)（約化質量倒數）
  ③ a = 1/(m_c+m_l)
  ④ a = 1/(2m_c)
對每個候選：由 A_D = 1.2025 MeV 解 t ⟹ 得 P₀（**預言**）⟹ 對比 P₀^req = K_D/(C(ħc/a)³)。
判準：只有①係「IDCM 自己嘅標度」；②③④屬常規量子力學標度，只作對照，唔可以揀嚟凑。
輸出：`data/11_scale_fixed_test.json`
"""
import json
import mpmath as mp

mp.mp.dps = 30
HBARC = mp.mpf("0.1973269804")
ALPHA = 1 / mp.mpf("137.035999084")
B, L = 4, 40
A_D, K_D = mp.mpf("1.2025"), mp.mpf("0.5287")
qh, mh, ml = mp.mpf(2) / 3, mp.mpf("1.27"), mp.mpf("0.30")


def shape(g, h):
    N = [mp.mpf(1)] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
    H = mp.zeros(L + 1, L + 1)
    for d in range(L + 1):
        H[d, d] = -g * mp.power(B, -d) - (h if d == 0 else 0)
    off1, off = -mp.sqrt(B), -mp.sqrt(B - 1)
    for d in range(L):
        H[d, d + 1] = off1 if d == 0 else off
        H[d + 1, d] = H[d, d + 1]
    ev, evec = mp.eigsy(H)
    i = 0    # ⚠ 基態（升序第一）—— 原 min|E| 誤揀激發態，已修
    v = [evec[d, i] for d in range(L + 1)]
    norm = sum(N[d] * v[d] ** 2 for d in range(L + 1))
    return (sum(N[d] * v[d] ** 2 * mp.power(B, -d) for d in range(L + 1)) / norm,
            v[0] ** 2 / norm)


def predict(a_fm):
    inv_a = HBARC / a_fm                                    # GeV
    C = (8 * mp.pi / 3) * ALPHA * qh / (16 * mh * ml)        # GeV^-2
    P0req = K_D / 1000 / (C * inv_a ** 3)                    # 由 K 反推所需 P₀
    def rA_of(logt):
        t = mp.e ** logt
        g = ALPHA * qh * inv_a / t
        h = C * inv_a ** 3 / t
        rA, P0 = shape(g, h)
        return ALPHA * qh * inv_a * rA * 1000 - A_D, P0     # MeV
    # 由 A 解 t（掃描 + brent）
    lo, hi = mp.log(1e-6), mp.log(1e6)
    flo = rA_of(lo)[0]
    fhi = rA_of(hi)[0]
    if flo * fhi > 0:
        return None
    logt = mp.findroot(lambda x: rA_of(x)[0], (lo, hi), solver="bisect", tol=mp.mpf("1e-18"))
    t = mp.e ** logt
    rA, P0 = shape(ALPHA * qh * inv_a / t, C * inv_a ** 3 / t)
    return {"a_fm": mp.nstr(a_fm, 8), "t_GeV": mp.nstr(t, 8), "rA": mp.nstr(rA, 8),
            "P0_pred": mp.nstr(P0, 8), "P0_req": mp.nstr(P0req, 8),
            "ratio": mp.nstr(P0 / P0req, 6)}


LAM = mp.mpf("0.04528487")                                   # GeV (O142 v16)
mu = mh * ml / (mh + ml)
CANDS = [("Lam_IDCM", HBARC / LAM),
         ("1/(2mu)", HBARC / (2 * mu)),
         ("1/(mh+ml)", HBARC / (mh + ml)),
         ("1/(2mh)", HBARC / (2 * mh))]

print("=" * 96)
print("M5 docs/10 — A2b: 結構標度固定 a ⟹ 由 A 定 t ⟹ 盲測 P₀（vs 由 K 反推所需值）")
print("=" * 96)
out = {}
for name, a in CANDS:
    r = predict(a)
    out[name] = r
    if r is None:
        print(f"\n[{name:10s}] a = {mp.nstr(a,6)} fm ⟹ A 方程無解（掃描範圍內）")
        continue
    print(f"\n[{name:10s}] a = {r['a_fm']} fm ⟹ t = {r['t_GeV']} GeV")
    print(f"   rA = {r['rA']}  |  P₀(預言) = {r['P0_pred']}  vs  P₀(由 K 需要) = {r['P0_req']}"
          f"  ⟹ ratio = {r['ratio']}")

print("""
[判讀原則]
 • 只有 `Lam_IDCM` 係 IDCM 自己導出嘅標度（結構）；其餘三個係常規量子力學標度，只作對照。
 • 若 `Lam_IDCM` 下 P₀(預言) 與 P₀(需要) 同量級 ⟹ 接觸項與 Coulomb 項由**同一 ψ** 自洽，
   即 docs/09 之「分母」得解；若差幾個數量級 ⟹ 兩體樹層模型不足以承載，A2 誠實失敗。
""")
with open("data/11_scale_fixed_test.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print("-> data/11_scale_fixed_test.json")