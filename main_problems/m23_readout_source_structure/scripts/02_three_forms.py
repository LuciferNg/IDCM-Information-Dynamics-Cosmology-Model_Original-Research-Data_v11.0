#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M23 scripts/02 — R3／R4：三型判定 ＋ 零參數／一參數分界 ＋ 格點閘

預註冊判準（跑前寫死，不得改；一律依 `01_problem_record §七`）：
  J0  格點閘（§七⑥）：任何候選之 ρ 若落 ℚ/41 格點 ⟹ 依 M12 R25 判**零資訊**，只作披露。
  J1  零參數（§七①）：候選式之整數須由已註冊結構固定；存在另一同樣合法之值 ⟹ 該整數自由。
  J2  盲測（§七③）：形式已於 `docs/01` 寫死；本檔只算值。
  J3  ✅ 升格條件（三項同時）：自由整數 = 0 ∧ |z_d| ≤ 3 ∧ |z_u| ≤ 3 ∧ MC null p < 0.05（`scripts/01` 之 10⁷）
  J4  「一參數」分界：即使准許 1 個自由整數，若 min χ 仍 > 3 ⟹ 一參數亦不閉合。
  J5  T2 之參數計數：4 個整數 − 1 個共同標度 = 有效 3；再減「分子分母同泛函」之要求 ⟹ 報告實際有效自由數。

分段：
  H. T3-a／T3-b 精確值（Gram 幾何，零參數）
  I. T2 全族枚舉（946² 精確）：密度、最佳、參數計數
  J. 一參數族 G1–G4（各一個自由整數）之 min χ
  K. T1 已註冊源對（引用 scripts/01 之 JSON，不重算）
  L. §七 六項逐項報

探針項（預期 FAIL，判詞不依賴）：
  I1 —— 曾預期大族 3σ 密度 ≥ 5%；實測 0.017%（族大但值域寬）⟹ 預期被否。
  J1 —— 曾預期「1 參數最好仍 > 3σ」；實測 2.169σ ⟹ 預期被否（可到 2.2σ，但 0 個 1σ 命中）。
"""
from fractions import Fraction as F
import json
import os
import sys

import mpmath as mp

mp.mp.dps = 50

N = 44
M = N + 1
SIG_D = mp.mpf('0.0351')
SIG_U = mp.mpf('0.03292')
TGT_D = mp.mpf('0.40350365')
TGT_U = mp.mpf('-0.02588837')

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
results = []


def chk(name, ok, detail=""):
    results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[FAIL] ") + name + ("  " + detail if detail else ""))

# ── 探針機制（預期假設／預註冊閘被否之項）：記錄事實，不影響 EXIT ──
probe_results = []


def chkp(name, ok, detail=""):
    probe_results.append((name, bool(ok)))
    print(("[PASS] " if ok else "[REC!] ") + name + ("  " + detail if detail else ""))


def zchi(rd, ru):
    zd = (mp.mpf(F(rd).numerator) / F(rd).denominator - TGT_D) / SIG_D if not isinstance(rd, mp.mpf) else (rd - TGT_D) / SIG_D
    zu = (mp.mpf(F(ru).numerator) / F(ru).denominator - TGT_U) / SIG_U if not isinstance(ru, mp.mpf) else (ru - TGT_U) / SIG_U
    return zd, zu, mp.sqrt(zd ** 2 + zu ** 2)


def delta(a, b, k):
    L = b - a
    if k <= a:
        return F(-k * L, M)
    elif k <= b:
        return F(45 * (k - a) - k * L, M)
    else:
        return F(L * (45 - k), M)


print("=== M23 scripts/02 — R3／R4：三型判定 ===")

# ── H. T3 關係型（Gram，零參數）────────────────────────────────────
print("── H. T3 關係型 ──")
# v_x(n) = φ_n(k_x) Δ_n ; 核 Gram 量由 scripts/00 G4 已精確錨定：
#   ‖v_u‖²=‖v_d‖²=1/18、‖v_s‖²=2/45、⟨v_u,v_s⟩=1/90、⟨v_d,v_s⟩=−1/45
IP_uu, IP_dd, IP_ss = F(1, 18), F(1, 18), F(2, 45)
IP_us, IP_ds, IP_ud = F(1, 90), F(-1, 45), F(-4, 5) * F(1, 18)
t3a_d, t3a_u = IP_ds / IP_ss, IP_us / IP_ss
nv = mp.sqrt(mp.mpf(1) / 18) * mp.sqrt(mp.mpf(2) / 45)      # ‖v_d‖‖v_s‖ = ‖v_u‖‖v_s‖
def to_mp(q):
    return mp.mpf(q.numerator) / mp.mpf(q.denominator)


t3b_d, t3b_u = to_mp(IP_ds) / nv, to_mp(IP_us) / nv
zd, zu, chi = zchi(t3a_d, t3a_u)
print("      T3-a 投影型：ρ_d = %s (=%s)  ρ_u = %s  ⟹ z=(%s, %s) χ=%s"
      % (t3a_d, float(t3a_d), t3a_u, mp.nstr(zd, 8), mp.nstr(zu, 8), mp.nstr(chi, 8)))
chk("H1  T3-a 之 ρ_d = −1/2 精確", t3a_d == F(-1, 2), "⟨v_d,v_s⟩/‖v_s‖² = %s / %s" % (IP_ds, IP_ss))
chk("H2  T3-a 之 ρ_u = +1/4 精確", t3a_u == F(1, 4), "⟨v_u,v_s⟩/‖v_s‖² = %s / %s" % (IP_us, IP_ss))
print("      T3-b 餘弦型：ρ_d = %s  ρ_u = %s" % (mp.nstr(t3b_d, 10), mp.nstr(t3b_u, 10)))
chk("H3  T3-b 之 ρ_d = −1/√5 且 ρ_u = 1/(2√5)（= cos 值，振幅被除掉）",
    abs(t3b_d + 1 / mp.sqrt(5)) < mp.mpf('1e-40') and abs(t3b_u - 1 / (2 * mp.sqrt(5))) < mp.mpf('1e-40'),
    "z=(%s, %s)" % (mp.nstr((t3b_d - TGT_D) / SIG_D, 8), mp.nstr((t3b_u - TGT_U) / SIG_U, 8)))
chi_t3 = min(mp.mpf(chi), mp.sqrt(((t3b_d - TGT_D) / SIG_D) ** 2 + ((t3b_u - TGT_U) / SIG_U) ** 2))
print("      ⟹ T3 最佳 χ = %s（零參數、形式唯一）" % mp.nstr(chi_t3, 8))

# ── I. T2 雙源全族 ───────────────────────────────────────────────
print("── I. T2 雙源（946²）──")
vals_dd = {}
vals_ds = {}
vals_du = {}
for a in range(1, N + 1):
    for b in range(a + 1, N + 1):
        vals_dd[(a, b)] = delta(a, b, 3)
        vals_ds[(a, b)] = delta(a, b, 4)
        vals_du[(a, b)] = delta(a, b, 2)
keys = list(vals_ds)
print("      δ 值集大小：δ_d %d ｜ δ_s %d ｜ δ_u %d"
      % (len(set(vals_dd.values())), len(set(vals_ds.values())), len(set(vals_du.values()))))

best = None
cnt1 = cnt3 = 0
tot = 0
chi_min = mp.inf
for k1 in keys:
    d1, u1 = vals_dd[k1], vals_du[k1]
    for k2 in keys:
        s2 = vals_ds[k2]
        if s2 == 0:
            continue
        tot += 1
        rd, ru = d1 / s2, u1 / s2
        zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
        zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
        c = mp.sqrt(zd ** 2 + zu ** 2)
        if abs(zd) < 1 and abs(zu) < 1:
            cnt1 += 1
        if abs(zd) < 3 and abs(zu) < 3:
            cnt3 += 1
        if c < chi_min:
            chi_min = c
            best = (k1, k2, rd, ru, zd, zu)
print("      T2 成員 = %d ｜ 最佳：(分子 %s, 分母 %s) ρ_d=%s (z=%s) ρ_u=%s (z=%s) χ=%s"
      % (tot, best[0], best[1], best[2], mp.nstr(best[4], 6), best[3], mp.nstr(best[5], 6), mp.nstr(chi_min, 6)))
dens2_3 = mp.mpf(cnt3) / tot
dens2_1 = mp.mpf(cnt1) / tot
print("      T2 密度：|z|<3 ∧ |z_u|<3 = %d/%d = %s ｜ |z|<1 = %d/%d = %s"
      % (cnt3, tot, mp.nstr(dens2_3 * 100, 6) + "%", cnt1, tot, mp.nstr(dens2_1 * 100, 6) + "%"))
chkp("I1  [探針·預期假設] T2 族密度 ≥ 5%？", dens2_3 >= mp.mpf('0.05'),
    "密度 = %s" % (mp.nstr(dens2_3 * 100, 6) + "%"))
chk("I2  T2 參數計數：4 整數 − 1 共同標度 = 有效 3 ≥ 1 ⟹ 零參數違反（§七①）", True,
    "自由有效整數 = 3（(a₁,b₁) 兩整數 ⊕ (a₂,b₂) 兩整數 − 標度 1）")

# ── J. 一參數分界（各一個自由整數）──────────────────────────────
print("── J. 一參數族（放寬至 1 個自由整數）──")
one_par = {
    "G1 (2,b) b∈3..44（錨 k_u=2）": [(2, b) for b in range(3, N + 1)],
    "G2 (a,3) a∈1..2（錨 k_d=3 為上源）": [(a, 3) for a in range(1, 3)],
    "G3 (3,b) b∈4..44（錨 k_d=3 為下源）": [(3, b) for b in range(4, N + 1)],
    "G4 (a,4) a∈1..3（錨 k_s=4 為上源）": [(a, 4) for a in range(1, 4)],
}
jp = {}
for name, pairs in one_par.items():
    bb = None
    for (a, b) in pairs:
        dd, du, dsn = delta(a, b, 3), delta(a, b, 2), delta(a, b, 4)
        rd, ru = dd / dsn, du / dsn
        zd = (mp.mpf(rd.numerator) / rd.denominator - TGT_D) / SIG_D
        zu = (mp.mpf(ru.numerator) / ru.denominator - TGT_U) / SIG_U
        c = mp.sqrt(zd ** 2 + zu ** 2)
        if bb is None or c < bb[2]:
            bb = (a, b, c, rd, ru)
    jp[name] = {"n": len(pairs), "best": (bb[0], bb[1]), "chi": str(bb[2]), "rho_d": str(bb[3]), "rho_u": str(bb[4])}
    print("      %-32s 成員 %2d ｜ 最佳 (a=%d,b=%d) ρ_d=%s ρ_u=%s χ=%s"
          % (name, len(pairs), bb[0], bb[1], bb[3], bb[4], mp.nstr(bb[2], 8)))
all_one = min(mp.mpf(v["chi"]) for v in jp.values())
chkp("J1  [探針·預期假設] 1 個自由整數下 min χ > 3？", all_one > 3,
    "min χ = %s" % mp.nstr(all_one, 8))

# ── K. T1 已註冊源對（引 scripts/01 JSON）──────────────────────────
print("── K. T1 已註冊源對（scripts/01 結果，不重算）──")
with open(os.path.join(DATA, "01_family_enum.json")) as fh:
    fam01 = json.load(fh)
for r in fam01["registered"]:
    _zd, _zu, _c = zchi(F(r["rho_d"]), F(r["rho_u"]))
    print("      (a=%d,b=%d)  ρ_d=%s (z=%s) ρ_u=%s (z=%s) χ=%s"
          % (r["a"], r["b"], r["rho_d"], mp.nstr(_zd, 8), r["rho_u"], mp.nstr(_zu, 8), r["chi"]))
print("      F-A 族：min χ = %s @ (a=%d,b=%d) ｜ 946 對中 |z|<1 = %d"
      % (fam01["min_chi"], fam01["min_chi_pair"]["a"], fam01["min_chi_pair"]["b"], fam01["hits_z1"]))
chk("K1  T1 全族（946 對）|z_d|<1 ∧ |z_u|<1 之成員 = 0（零參數下無閉合點）",
    fam01["hits_z1"] == 0, "hits(z<1) = %d/946" % fam01["hits_z1"])

# ── L. §七 六項逐項報 ────────────────────────────────────────────
print("── L. 驗證標準（§七）六項 ──")
rows = [
    ("① 零參數", "T1: 2 整數（(2,3) 由 M14 軸固定 ⟹ 0 自由；但 P3 混段位）｜ T2: 3 自由 🔴 ｜ T3: 0 自由 ✅"),
    ("② MC null 10⁷", "T1 p = %s（`scripts/01`）｜ T2 密度 %s ｜ T3 零參數無 null 需求（值唯一）"
     % (fam01["mc_null"]["p"], mp.nstr(dens2_3 * 100, 4) + "%")),
    ("③ 盲測", "docs/01 先寫死三型語法（零數值）；本檔才算值 ⟹ 次序合規"),
    ("④ 剛性", "N=43／45 之族級重算見 `scripts/03`"),
    ("⑤ 防重複", "T1 之可達集 vs M20 單鍵 ±1 類：見 `scripts/03`（族級 min χ 對照 ＋ 結構層理由）"),
    ("⑥ 格點閘", "T1/T2 之 ρ 全部有理；落 ℚ/41 者按 M12 R25 判零資訊（只披露，不作證據）"),
]
for a, b in rows:
    print("      %-14s %s" % (a, b))
# 格點狀態
grid_hits = [r for r in fam01["registered"] if (F(r["rho_d"]) * 41).denominator == 1]
print("      F-B 落 ℚ/41 格點之成員：%d/3（%s）⟹ 依 §七⑥ 判零資訊"
      % (len(grid_hits), [("(%d,%d)" % (r["a"], r["b"])) for r in grid_hits]))
chk("L1  格點閘已報（不作證據）", True, "F-B 落格點 = %d/3" % len(grid_hits))

out = {
    "t3a": {"rho_d": str(t3a_d), "rho_u": str(t3a_u), "chi": str(chi)},
    "t3b": {"rho_d": mp.nstr(t3b_d, 20), "rho_u": mp.nstr(t3b_u, 20),
            "chi": mp.nstr(mp.sqrt(((t3b_d - TGT_D) / SIG_D) ** 2 + ((t3b_u - TGT_U) / SIG_U) ** 2), 10)},
    "t3_best_chi": str(chi_t3),
    "t2": {"members": tot, "best": {"pair_num": list(best[0]), "pair_den": list(best[1]),
                                    "rho_d": str(best[2]), "rho_u": str(best[3]),
                                    "z_d": mp.nstr(best[4], 10), "z_u": mp.nstr(best[5], 10),
                                    "chi": str(chi_min)},
           "density_z3": str(dens2_3), "density_z1": str(dens2_1), "hits_z1": cnt1, "hits_z3": cnt3},
    "one_param": jp,
}
with open(os.path.join(DATA, "02_three_forms.json"), "w") as fh:
    json.dump(out, fh, indent=2, ensure_ascii=False)

npass = sum(1 for _, ok in results if ok)
print("")
print("──────────────────────────────────────────")
nrf = sum(1 for _, ok in probe_results if not ok)
print("結果: %d/%d PASS ＋ %d 記錄 FAIL（探針／預期假設被否：判詞不依賴）" % (npass, len(results), nrf))
print("──────────────────────────────────────────")
sys.exit(0 if npass == len(results) else 1)
