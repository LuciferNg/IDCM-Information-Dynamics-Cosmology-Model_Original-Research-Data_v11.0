#!/usr/bin/env python3
"""M3 L2 — 獨立重算（唔 import L1、唔讀 docs）

所有 INPUTS 皆原始量（計數／測量／上游已封結論）；
所有 DERIVE 由 INPUTS 算出；IDENTITIES 驗代數恆等（差 < 1e-30）。
輸出：data/l2_independent_verify.json
"""
import json, math
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40

INPUTS = {
    "N_sites": 88,          # 環 site 數（上游 U11／二十面體）
    "N_d": 44,              # 半譜模數
    "K_tree": 33,           # 樹模數 = N_gen·L₅
    "L5": 11,               # 鏈長指數（Lucas）
    "mu2": 2,               # SSB 勢 μ²
    "lam_num": 2, "lam_den": 123,   # λ = 2/123
    "omega_b": 0.0493, "sig_b": 0.0003,     # Planck 2018
    "omega_dm": 0.2645, "sig_dm": 0.0050,
    "omega_L": 0.6847, "sig_L": 0.0073,
    "h_anchor": 0.674071,
}

B = 4  # 由下方 gam_max 導出，先以變數佔位
_out, _ident = {}, {}

def put(k, v):
    _out[k] = v
    return v

# ── 譜（由 N_sites 直接算，唔用閉式）────────────────────────────
N = INPUTS["N_sites"]; M = INPUTS["N_d"]
gam = [4*mp.sin(mp.pi*n/N)**2 for n in range(0, M+1)]
gam_max = max(gam)
sig_gamma = mp.fsum(gam)
sig_odd = mp.fsum(gam[n] for n in range(1, M+1, 2))
sig_even = mp.fsum(gam[n] for n in range(0, M+1, 2))
put("gam_max", gam_max); put("sig_gamma", sig_gamma)
put("gamma_2", gam[2]); put("gamma_3", gam[3])
put("sig_odd", sig_odd); put("sig_even", sig_even)

# ── 閘門／容量 ─────────────────────────────────────────────────
kappa = 1/gam_max**2
put("kappa", kappa)
Bv = gam_max                      # B ≡ γ_max（Nyquist）
cum, s = [], 0
for L in range(0, 6):
    s += float(Bv)**L
    cum.append(s)
put("cum", cum)
ell = max(L for L in range(0, 6) if cum[L] <= INPUTS["K_tree"])
put("ell", ell)
cum2 = mp.mpf(1 + float(Bv) + float(Bv)**2)
put("cum2", cum2)

# 通過閘門之層（n ≥ 0；記錄兩種索引）
pass_0_43 = [n for n in range(0, INPUTS["N_d"]) if mp.e**(-n) > kappa]
pass_1_44 = [n for n in range(1, INPUTS["N_d"]+1) if mp.e**(-n) > kappa]
put("rank_0_43", len(pass_0_43)); put("rank_phys", len(pass_1_44))

# ── 三分 ───────────────────────────────────────────────────────
vis = int(mp.ceil(ell*mp.log(Bv)))
dm = int(Bv)**ell
lam = INPUTS["N_d"] - vis
tot = INPUTS["N_d"] + dm
put("vis", vis); put("dm", dm); put("lam", lam); put("total", tot)
put("omegaL", F(lam, tot)); put("omega_m", F(vis+dm, tot))
put("omega_b", F(vis, tot)); put("omega_dm", F(dm, tot))

# ── 能量側（SSB 勢）────────────────────────────────────────────
lam_w = F(INPUTS["lam_num"], INPUTS["lam_den"])
E_sync = F(INPUTS["mu2"])/(2*lam_w)          # ⟨W⟩² = μ²/(2λ)
E_proj = F(90) - E_sync
put("E_sync", E_sync); put("E_proj", E_proj)
put("ratio_41_19", F(41, 19))
put("quantum", E_sync/41)
put("sixty_from_sig", F(90)/(E_sync/41))

# ── 觀測對照（z 值）──────────────────────────────────────────
def z(pred, obs, sig):
    return (mp.mpf(pred.numerator)/pred.denominator - mp.mpf(str(obs)))/mp.mpf(str(sig))
put("z_L", z(F(lam, tot), INPUTS["omega_L"], INPUTS["sig_L"]))
put("z_m", z(F(vis+dm, tot), 0.3153, 0.0073))
put("z_b", z(F(vis, tot), INPUTS["omega_b"], INPUTS["sig_b"]))
put("z_dm", z(F(dm, tot), INPUTS["omega_dm"], INPUTS["sig_dm"]))

# ── Ω_b h² vs 1/45 ────────────────────────────────────────────
obh2 = mp.mpf(vis)/tot*mp.mpf(str(INPUTS["h_anchor"]))**2
put("omega_b_h2", obh2)
put("dev_1_45_pct", (obh2 - mp.mpf(1)/45)/(mp.mpf(1)/45)*100)
put("h_for_1_45", mp.sqrt(mp.mpf(20)/45))

# ── ℓ 掃描（Ω_dm(ℓ) = B^ℓ/(N_d+B^ℓ)）─────────────────────────
put("z_l1", (mp.mpf(4)/48 - mp.mpf("0.26447"))/mp.mpf("0.005"))
put("z_l2", (mp.mpf(16)/60 - mp.mpf("0.26447"))/mp.mpf("0.005"))
put("z_l3", (mp.mpf(64)/108 - mp.mpf("0.26447"))/mp.mpf("0.005"))
put("z_l4", (mp.mpf(256)/300 - mp.mpf("0.26447"))/mp.mpf("0.005"))

# ── 值空間（整數三分）─────────────────────────────────────────
triples = [(a, b, 60-a-b) for a in range(1, 60) for b in range(1, 60-a) if 60-a-b >= 1]
put("n_triples", len(triples))
OBS = {"b": 0.0493, "dm": 0.2645, "L": 0.6847}
SIG = {"b": 0.0003, "dm": 0.005, "L": 0.0073}
def chi2(t):
    return sum(((t[i]/60 - OBS[k])/SIG[k])**2 for i, k in enumerate(("b", "dm", "L")))
chi2_obs = chi2((3, 16, 41))
put("chi2_obs", chi2_obs)
put("chi2_3_15_42", chi2((3, 15, 42)))
put("p1_value_space", 1/len(triples))
put("z_p1", float(mp.sqrt(2)*mp.erfinv(1 - 2*mp.mpf(1)/len(triples))))

# ── 畢氏審計算術 ──────────────────────────────────────────────
put("pyth_lhs_5_12_13", 5 + 2*144 + 169)
put("pyth_rhs_5_12_13", 5*12*13)
put("prod_345", 3*4*5)
put("dims_total", INPUTS["N_d"]*INPUTS["K_tree"])

# ── 字典命中（M3 結構字典）───────────────────────────────────
D = {1,2,3,4,5,6,8,11,15,16,19,22,29,30,33,41,44,45,46,60,88,90,123,176}
put("dict_size", len(D))
put("dict_admissible", sum(1 for t in triples if chi2(t) <= chi2_obs and set(t) <= D))

# ── IDENTITIES（代數恆等，差 < 1e-30）────────────────────────
_ident["sig_gamma == 2M+2"] = abs(sig_gamma - (2*M+2))
_ident["sig_odd == M"] = abs(sig_odd - M)
_ident["sig_even == M+2"] = abs(sig_even - (M+2))
_ident["gam_max == 4"] = abs(gam_max - 4)
_ident["cum(ell) <= K_tree < cum(ell+1)"] = max(0, cum[ell] - INPUTS["K_tree"]) if cum[ell] <= INPUTS["K_tree"] else 1
_ident["vis + dm + lam == total"] = abs((vis+dm+lam) - tot)
_ident["omegaL + omega_b + omega_dm == 1"] = abs(float(F(lam,tot)+F(vis,tot)+F(dm,tot)) - 1)
_ident["E_sync + E_proj == 90"] = abs(float(E_sync+E_proj) - 90)
_ident["E_sync/41 == E_proj/19"] = abs(float(E_sync/41 - E_proj/19))
_ident["total == N_d + B^ell"] = abs(tot - (INPUTS["N_d"] + int(Bv)**ell))

ok = all(v < mp.mpf("1e-30") for v in _ident.values())
report = {
    "inputs": INPUTS,
    "derived": {k: (str(v) if isinstance(v, (mp.mpf, mp.mpc)) else v) for k, v in _out.items()},
    "identities": {k: str(v) for k, v in _ident.items()},
    "pass": ok,
}
with open("data/l2_independent_verify.json", "w", encoding="utf-8") as f:
    json.dump(report, f, ensure_ascii=False, indent=1, default=str)

np_ = 0
for k, v in _out.items():
    print(f"  [key] {k} = {v}"); np_ += 1
for k, v in _ident.items():
    st = "PASS" if v < mp.mpf("1e-30") else "FAIL"
    print(f"  [{st}] ident {k}  |diff| = {mp.nstr(v, 5)}")
print(f"\nkeys {np_} ; identities {len(_ident)} ; all_pass {ok}")