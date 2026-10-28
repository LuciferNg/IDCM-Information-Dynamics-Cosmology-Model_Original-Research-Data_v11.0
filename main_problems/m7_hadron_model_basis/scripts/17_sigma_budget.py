"""M7 scripts/17 — 回核 M4 docs/40 之誤差預算更正（σ_X 漏計）
==========================================================================================================
M4 docs/40 ④ 自主更正：原 z = −15.98σ 漏計 σ_X = 0.094（X = m_d − m_u 之外部誤差）。
  σ_pred = (5/3)·σ_X ；σ_tot = √(σ_em² + σ_pred²) ；z = (em_pred − em)/σ_tot ⇒ 宣稱 −9.76
本檔獨立重算並判定：① 新 z；② 是否仍 KILL；③ 對本 M 既有結論（(1−ε) 退化、帶內最保守）之影響；
  ④ 為何 −15.98 仍可保留為「X 視作精確」之約定值（並說明本 M 之 γ 載體推導已 NO-GO ⟹ 必須用外部 X）。
輸出 data/17_sigma_budget.json
"""
import json
import os
from fractions import Fraction as F

from mpmath import mp, mpf, sqrt, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "17_sigma_budget.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


EM = mpf("2.264589")          # scripts/14（外部 em）
PRED = mpf("0.332667")        # M1 @ n_s=1：Δm + X·5/3（X = 2.56 精確）
S_EM = mpf("0.120863")        # ε 帶半寬（scripts/14）
X = mpf("2.56")
S_X = mpf("0.094")            # FLAG24：X = 2.56 ± 0.094（M4 docs/18 §2）
DM_K = mpf("-3.9340")

print("=" * 96)
print("M7 scripts/17 — 誤差預算回核（M4 docs/40 ④）")
print("=" * 96)

s_pred = mpf(5) / 3 * S_X
s_tot = sqrt(S_EM ** 2 + s_pred ** 2)
z_tot = (PRED - EM) / s_tot
z_econv = (PRED - EM) / S_EM
print(f"    σ_pred = (5/3)·σ_X = {nstr(s_pred, 8)}")
print(f"    σ_tot  = √(σ_em² + σ_pred²) = √({nstr(S_EM**2,8)} + {nstr(s_pred**2,8)}) = {nstr(s_tot, 8)}")
print(f"    z（含 σ_X）= {float(z_tot):+.3f}；z（X 視作精確，約定值）= {float(z_econv):+.2f}")

ok("σ_pred = 0.156667（M4 值）", abs(s_pred - mpf("0.156667")) < mpf("5e-7"), nstr(s_pred, 8), "0.156667")
ok("σ_tot = 0.197869（M4 值）", abs(s_tot - mpf("0.197869")) < mpf("5e-7"), nstr(s_tot, 9), "0.197869")
ok("z 含 σ_X = −9.76（M4 宣稱）", abs(z_tot + mpf("9.76")) < mpf("0.01"), f"{float(z_tot):+.3f}", "-9.76")
ok("判決仍 KILL（|z| ≫ 3）", abs(z_tot) > 3, f"{float(z_tot):+.2f}", "|z|>3", "誤差預算不救預言")

# 帶內最保守（ε = 0.67 ⟹ em 最小 2.124642）
EM_LO = mpf("2.124642")
z_band = (PRED - EM_LO) / s_tot
print(f"    ε 帶內最保守（em = 2.124642）⟹ z = {float(z_band):+.2f}")
ok("帶內最保守仍 KILL", abs(z_band) > 3, f"{float(z_band):+.2f}", "|z|>3")

# (1−ε) 退退化在新区間下之距離
EM_WRONG = mpf("0.279893")
d_wrong = (PRED - EM_WRONG) / s_tot
print(f"    預言 vs (1−ε) 誤讀值：{float(d_wrong):+.2f}σ（舊 σ_em 下 +0.44σ）")
ok("(1−ε) 退化診斷在含 σ_X 後仍成立（<1σ）", abs(d_wrong) < 1, f"{float(d_wrong):+.2f}σ", "<1σ",
   "⟹ 符號層診斷不受誤差預算影響")

# 為何 −15.98 可保留：γ 載體精確 X（本 M）已 NO-GO
X_gamma = mpf("2.5461523")
z_gamma = (DM_K + X_gamma * mpf(5) / 3 - EM) / S_EM
print(f"    若用本 M γ 載體之 X = {nstr(X_gamma,8)}（無外部 σ）：z = {float(z_gamma):+.2f}")
ok("γ 載體 X 分支給 z ≈ −16（但該分支已 NO-GO ⟹ 不得引用）",
   abs(z_gamma + mpf("16")) < mpf("0.2"), f"{float(z_gamma):+.2f}", "≈ −16",
   "docs/20 §七：u/d 指派 NO-GO ⟹ 外部 X 為唯一可用 ⟹ canonical z = −9.76")

print("\n" + "=" * 96)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(s_pred=float(s_pred), s_tot=float(s_tot), z_tot=float(z_tot),
                        z_econv=float(z_econv), z_band_tot=float(z_band),
                        d_wrong_tot=float(d_wrong), z_gamma=float(z_gamma)),
           verdict="accept M4 ④: canonical z=-9.76 (with sigma_X); -15.98 kept as X-exact convention")
json.dump(OUT, open(os.path.join(M7, "data", "17_sigma_budget.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/17_sigma_budget.json")