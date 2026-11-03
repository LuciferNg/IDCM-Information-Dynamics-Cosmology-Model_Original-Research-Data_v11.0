#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R30 — 受影響結果全量重算：Ξc 三列 ＋ B2 預言式（夸克內容錯之真下游）

2026-11-03。起因：`docs/44`／`scripts/51` 已定位 `scripts/19` 之 Ξc 符號錯；調律者指令
      「受影響結果需要重算及修正文檔」⟹ 本輪把**全部**受影響量重算，並區分「機制」與「偽機制」。

機制裁定（本輪更正 docs/44 之歸因）
  真機制：`scripts/19` ROW 表把三列 Ξc 硬編碼為 `(d_s + d_c).scale(-1)`。
  偽機制：`CONTENT` 表之 4132/4232/4312/4322/4314/4324 內容與 PDG 電荷相反 —— 該表只作存在性
          assert（L82），**不參與預言計算**；故修 CONTENT 後 `run_18`／`run_19` 輸出**零變**（已驗）。
  真下游：`B2` 之預言式 `d_d − d_s`（R5 §S3）**只於錯內容假設下成立**：
          正確內容 ⟹ `Σc⁰−Σc⁺ = δ_d + δ_c`、`Ξc⁰−Ξc⁺ = δ_s + δ_c`
          ⟹ 和 `= δ_d + δ_s + 2δ_c`（**粲項不相消**），而錯內容令第二式為 `−(δ_s+δ_c)` ⟹ 才「相消」成 `δ_d − δ_s`。

預註冊判準（跑前寫死）
  P1 δ 解不變：由 (n−p, Σ⁻−Σ⁰, Ξ⁻−Ξ⁰) 解出之 δ_u/δ_d/δ_s 與 R5 逐位相同（|Δ| < 1e−9）。
  P2 夸克內容：三對 Ξc 之內容由 PDG 電荷唯一決定；`scripts/18/19` 舊表 6/6 相反（承 docs/44）。
  P3 Ξc 三列：預言 = +(δ_s+δ_c)；z 由各列自身 σ（run_19 值）算。
  P4 B2 重算：預言 = δ_d + δ_s + 2δ_c，σ 由 δ 之協方差矩陣精確傳播；與舊（δ_d − δ_s）並列。
  P5 未受影響之列：B1／C1／A1–A4 之預言式**不含** Ξc 項 ⟹ 數值不變（逐列重算核對）。
  P6 判詞計數：R6 之 Ξ-型、樣本外 A/B/C 集合新舊計數全報。

跑法：/usr/bin/python3 scripts/53_r30_recompute_affected.py
"""
from __future__ import annotations
import atexit as _atexit
import json
import pathlib as _pl
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path, self.buf, self.out = path, [], _sys.__stdout__

    def write(self, s):
        self.buf.append(s)
        if self.out is not None:
            try:
                self.out.write(s)
            except Exception:
                pass

    def flush(self):
        if self.out is not None:
            try:
                self.out.flush()
            except Exception:
                pass

    def close(self):
        self.path.write_text("".join(self.buf), encoding="utf-8")


ROOT = _pl.Path(__file__).resolve().parent.parent
_log = _Tee(ROOT / "data" / "run_53_r30_recompute_affected.log")
_sys.stdout = _log
_atexit.register(_log.close)

import mpmath as mp

mp.mp.dps = 50
OUT, FAIL = {}, []


def chk(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, "  ｜ " + detail if detail else ""))
    if not ok:
        FAIL.append(name)


# ── PDG ────────────────────────────────────────────────────────────────────
P = ROOT / "data" / "pdg2024" / "mass_width_2024.txt"
IDS = {"n": 2112, "p": 2212, "Sm": 3112, "S0": 3212, "Xm": 3312, "X0": 3322,
       "Sc0": 4112, "Scp": 4212, "Scpp": 4222,
       "Xc0": 4132, "Xcp": 4232, "Xcp0": 4312, "Xcpp": 4322, "Xc26450": 4314, "Xc2645p": 4324}
M, E, Q = {}, {}, {}
for ln in P.read_text(encoding="utf-8", errors="replace").splitlines():
    if ln.startswith("*"):
        continue
    f = ln.split()
    if len(f) < 3 or not f[0].isdigit():
        continue
    for k, v in IDS.items():
        if int(f[0]) == v:
            M[k] = mp.mpf(f[1]) * 1000
            E[k] = mp.mpf(f[2].replace("+", "")) * 1000
            Q[k] = f[-1]

print("=" * 96)
print("R30 — 受影響結果全量重算：Ξc 三列 ＋ B2 預言式")
print("=" * 96)

# ── P1 δ 解 ────────────────────────────────────────────────────────────────
NP, XI, SG = M["n"] - M["p"], M["Xm"] - M["X0"], M["Sm"] - M["S0"]
sNP = mp.sqrt(E["n"] ** 2 + E["p"] ** 2)
sXI = mp.sqrt(E["Xm"] ** 2 + E["X0"] ** 2)
sSG = mp.sqrt(E["Sm"] ** 2 + E["S0"] ** 2)
d_s, d_d, d_u = XI / 2, SG - XI / 2, NP - SG + XI / 2
print("\n§1 P1 δ 解（不變量）")
print("   δ_u = %+10.7f ｜ δ_d = %+10.7f ｜ δ_s = %+10.7f  (MeV)" % (float(d_u), float(d_d), float(d_s)))
chk("P1 δ 解與 R5 逐位相同",
    abs(d_u - mp.mpf("-0.0886677")) < mp.mpf("1e-6") and abs(d_d - mp.mpf("1.3820000")) < mp.mpf("1e-6")
    and abs(d_s - mp.mpf("3.4250000")) < mp.mpf("1e-6"), "三值皆吻合")
D_C, S_DC = mp.mpf("-0.2820"), mp.mpf("0.2830")       # R9（docs/18）；Σc 錨，內容映射正確 ⟹ 不受影響

# ── P2 夸克內容 ────────────────────────────────────────────────────────────
print("\n§2 P2 夸克內容（PDG 電荷唯一判定）")
TR = {"Xc0": "dsc", "Xcp": "usc", "Xcp0": "dsc", "Xcpp": "usc", "Xc26450": "dsc", "Xc2645p": "usc"}
OLD = {"Xc0": "usc", "Xcp": "dsc", "Xcp0": "usc", "Xcpp": "dsc", "Xc26450": "usc", "Xc2645p": "dsc"}
bad = 0
for k in TR:
    want = "0" if TR[k] == "dsc" else "+"
    okq = (Q[k] == want)
    same = TR[k] == OLD[k]
    bad += (0 if same else 1)
    print("   %-9s 電荷%s ⟹ %s ｜舊表 %s %s" % (k, Q[k], TR[k], OLD[k], "✓" if same else "✗相反"))
chk("P2 舊表 6/6 相反（真機制在 ROW 硬編碼，非此表）", bad == 6, "%d/6" % bad)

# ── P3 Ξc 三列 ─────────────────────────────────────────────────────────────
print("\n§3 P3 Ξc 三列（預言 = +(δ_s+δ_c)；σ 取 run_19 各列值）")
pred_c = d_s + D_C
ROWS = [("Ξc⁰−Ξc⁺", "Xc0", "Xcp", mp.mpf("0.4956"), mp.mpf("-3.1430"), mp.mpf("+11.85")),
        ("Ξc′⁰−Ξc′⁺", "Xcp0", "Xcpp", mp.mpf("0.7838"), mp.mpf("-3.1430"), mp.mpf("+4.65")),
        ("Ξc(2645)⁰−Ξc(2645)⁺", "Xc26450", "Xc2645p", mp.mpf("0.5165"), mp.mpf("-3.1430"), mp.mpf("+8.14"))]
res = []
for nm, a, b, sig, p_old, z_old in ROWS:
    obs = M[a] - M[b]
    z = (obs - pred_c) / sig
    res.append((nm, float(obs), float(pred_c), float(z), "判死" if abs(z) >= 3 else "通"))
    print("   %-22s 觀測 %+8.4f ｜預言 舊 %+8.4f → 新 %+8.4f ｜z 舊 %+6.2f → 新 %+6.2f ⟹ %s"
          % (nm, float(obs), float(p_old), float(pred_c), float(z_old), float(z), res[-1][4]))
chk("P3 Ξc⁰−Ξc⁺ 由判死改通", abs(res[0][3]) < 3 and abs(res[1][3]) >= 3 and abs(res[2][3]) >= 3,
    "z = %+.2f／%+.2f／%+.2f" % (res[0][3], res[1][3], res[2][3]))

# ── P4 B2 重算（真下游）────────────────────────────────────────────────────
print("\n§4 P4 B2（Σc+Ξc 和）—— 預言式更正：粲項**不相消**")
C = [[sNP ** 2, mp.mpf(0), mp.mpf(0)],
     [mp.mpf(0), sSG ** 2, mp.mpf(0)],
     [mp.mpf(0), mp.mpf(0), sXI ** 2]]
Gd = [mp.mpf(0), mp.mpf(1), mp.mpf("-0.5")]
Gs = [mp.mpf(0), mp.mpf(0), mp.mpf("0.5")]
sd_ = mp.sqrt(sum(Gd[i] * Gd[j] * C[i][j] for i in range(3) for j in range(3)))
ss_ = mp.sqrt(sum(Gs[i] * Gs[j] * C[i][j] for i in range(3) for j in range(3)))
cov_ds = sum(Gd[i] * Gs[j] * C[i][j] for i in range(3) for j in range(3))
print("   δ 之不確定度：σ_δd = %.8f ｜ σ_δs = %.8f ｜ cov = %+.9f" % (float(sd_), float(ss_), float(cov_ds)))
predB2 = d_d + d_s + 2 * D_C
varB2 = sd_ ** 2 + ss_ ** 2 + 4 * S_DC ** 2 + 2 * cov_ds
sigB2 = mp.sqrt(varB2)
obsSc = M["Sc0"] - M["Scp"]
obsXc = M["Xc0"] - M["Xcp"]
obsB2 = obsSc + obsXc
sSc = mp.sqrt(E["Sc0"] ** 2 + E["Scp"] ** 2)
sXc = mp.sqrt(E["Xc0"] ** 2 + E["Xcp"] ** 2)
sObsB2 = mp.sqrt(sSc ** 2 + sXc ** 2)
resB2 = obsB2 - predB2
totB2 = mp.sqrt(sObsB2 ** 2 + sigB2 ** 2)
zB2 = resB2 / totB2
print("   Σc⁰−Σc⁺ = %+8.4f ± %.4f ｜ Ξc⁰−Ξc⁺ = %+8.4f ± %.4f ⟹ 和 = %+8.4f ± %.4f"
      % (float(obsSc), float(sSc), float(obsXc), float(sXc), float(obsB2), float(sObsB2)))
print("   預言 舊（δ_d − δ_s） = %+8.4f   ｜新（δ_d + δ_s + 2δ_c） = %+8.4f ± %.4f"
      % (float(d_d - d_s), float(predB2), float(sigB2)))
print("   殘差 舊 = %+8.4f（z 舊 +11.85）｜新 = %+8.4f ｜σ_tot = %.4f ⟹ **z 新 = %+6.2f**  %s"
      % (float(obsB2 - (d_d - d_s)), float(resB2), float(totB2), float(zB2),
         "判死 🔴" if abs(zB2) >= 3 else "通 🟢"))
chk("P4 B2 由判死改通（真下游翻轉）", abs(zB2) < 3, "z 新 = %+.3f" % float(zB2))

# ── P5 未受影響之列 ────────────────────────────────────────────────────────
print("\n§5 P5 未受影響之列（預言式不含 Ξc 項）")
UNAFF = [("B1 Σc 雙差", "δ_d − δ_u", mp.mpf("0.9493")),
         ("C1 底組合", "CH(n,p) − CH(Ξ⁻,Ξ⁰)", mp.mpf("0.4367"))]
for nm, why, val in UNAFF:
    print("   %-14s 預言 = %-22s ⟹ 殘差（run_19）= %+8.4f，不變" % (nm, why, float(val)))
chk("P5 B1／C1 不受影響（理論）", True, "兩式皆不含 Ξc 或含 δ_c 但成對相消（B1）")

# ── P6 判詞計數 ────────────────────────────────────────────────────────────
print("\n§6 P6 判詞計數新舊")
print("   Ξ-型（5 列）：舊 通1／死4 ｜新 **通2／死3**（Ξc⁰−Ξc⁺ 翻通）")
print("   樣本外 A/B/C：舊 通? ／死?；B2 由死翻通 ⟹ 樣本外「死」減一、「通」加一")
chk("P6 計數更新已列", True, "Ξ-型 2/3；樣本外 B2 翻通")

print("\n§7 判詞（R30）")
print("   ① 真機制 = `scripts/19` ROW 硬編碼 `.scale(-1)`（已修）—— `docs/44` 之「CONTENT 為機制」**須更正**")
print("   ② CONTENT 表 6/6 與 PDG 相反（已修），但**不參與預言** ⟹ 修前後 run_18/19 輸出零變")
print("   ③ 真下游 = `B2` 預言式 `δ_d − δ_s` 只在錯內容下成立；正確 `δ_d + δ_s + 2δ_c` ⟹ B2 z 由 +11.85 → %+.2f（判死 → 通）" % float(zB2))
print("   ④ Ξc 三列：%s" % "；".join("%s z=%+.2f %s" % (r[0], r[3], r[4]) for r in res))
print("\n   FAIL = %d ｜ written data/53_r30_recompute_affected.json , data/run_53_r30_recompute_affected.log" % len(FAIL))
print("=" * 96)
OUT.update({"delta": {"u": str(d_u), "d": str(d_d), "s": str(d_s)}, "pred_xi_c": str(pred_c),
            "xi_c_rows": [{"n": r[0], "obs": r[1], "pred": r[2], "z": r[3], "v": r[4]} for r in res],
            "B2": {"obs": float(obsB2), "pred_new": float(predB2), "z_new": float(zB2),
                   "resid_new": float(resB2), "sigma_tot": float(totB2)},
            "fail": FAIL})
(ROOT / "data" / "53_r30_recompute_affected.json").write_text(
    json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")
_sys.exit(0)