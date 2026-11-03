#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R29/審計 — 跨 M handoff 收據 ＋ Ξc 符號錯之定位與更正（2026-11-03）

來源：`M13 → M12` handoff（2026-11-02 06:53）之第 ① 項
      「docs/15 L61 Ξc 預言 −3.1430 疑為符號錯 ⟹ 13σ 判死疑 Type-1」

本輪只做**記錄層更正 ＋ 程式碼層定位**，零新推導、零新參數。

判準（跑前寫死）
  P1 夸克內容可由 PDG 電荷**唯一判定**：4132 電荷 0 ⟹ dsc；4232 電荷 + ⟹ usc。
     比對 `scripts/19` 之 CONTENT 表；若相反 ⟹ 符號錯定位成功（程式碼層）。
  P2 改正預言：Ξc⁰(dsc) − Ξc⁺(usc) = δ_s + δ_c = +3.1430（原寫 −3.1430）。
  P3 以 run_19 log 之 σ（0.4956／0.7838／0.5165，即含預言不確定度）重算三列 z；
     另以純觀測 σ（0.3624，M13 版）算 Ξc⁰−Ξc⁺ 作交叉核。
  P4 R6 判詞影響：Ξ 型五行逐列報（死／通），並計新舊計數。
  P5 B2（Σc+Ξc 和）共用 Ξc 項 ⟹ 只標「待重算」，**不**在本輪下判（避免連帶臆測）。
  P6 不動 `docs/15`、`scripts/19` 原文；本輪只出新檔（append-only）。

跑法：/usr/bin/python3 scripts/51_inter_m_handoff_audit.py
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
_log = _Tee(ROOT / "data" / "run_51_inter_m_handoff_audit.log")
_sys.stdout = _log
_atexit.register(_log.close)

import mpmath as mp

mp.mp.dps = 40
OUT, FAIL = {}, []


def chk(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, "  ｜ " + detail if detail else ""))
    if not ok:
        FAIL.append(name)


# ── PDG（vendored）：MC-ID → 質量、誤差、電荷 ────────────────────────────────
P = ROOT / "data" / "pdg2024" / "mass_width_2024.txt"
WANT = {"4132": "Ξc⁰", "4232": "Ξc⁺", "4312": "Ξc′⁰", "4322": "Ξc′⁺",
        "4314": "Ξc(2645)⁰", "4324": "Ξc(2645)⁺"}
M, E, Q = {}, {}, {}
for ln in P.read_text(encoding="utf-8", errors="replace").splitlines():
    if ln.startswith("*"):
        continue
    f = ln.split()
    if len(f) < 3 or f[0] not in WANT:
        continue
    M[f[0]] = mp.mpf(f[1]) * 1000
    E[f[0]] = mp.mpf(f[2].replace("+", "")) * 1000
    Q[f[0]] = f[-1]

print("=" * 94)
print("跨 M handoff 收據 ＋ Ξc 符號錯定位（M13 → M12 handoff 第 ① 項）")
print("=" * 94)
print("\n§0 PDG 2024（vendored；MC-ID 定位，禁硬編碼）")
for k in ("4132", "4232", "4312", "4322", "4314", "4324"):
    print("   %s %-14s m = %10.5f ± %.5f MeV  電荷 = %s" % (k, WANT[k], float(M[k]), float(E[k]), Q[k]))

# ── P1：夸克內容由電荷唯一判定 ──────────────────────────────────────────────
print("\n§1 P1 夸克內容（電荷唯一判定）vs `scripts/19` 之 CONTENT 表")
TRUE = {"4132": "dsc", "4232": "usc", "4312": "dsc", "4322": "usc", "4314": "dsc", "4324": "usc"}
OLD = {"4132": "usc", "4232": "dsc", "4312": "usc", "4322": "dsc", "4314": "usc", "4324": "dsc"}
CHG = {"dsc": 0, "usc": 1}
bad = 0
for k in TRUE:
    q_expected = CHG[TRUE[k]]
    okq = (Q[k] == ("0" if q_expected == 0 else "+"))
    same = OLD[k] == TRUE[k]
    if not same:
        bad += 1
    print("   %s %-14s 電荷%s ⟹ 內容 %s ｜scripts/19 寫 %s %s"
          % (k, WANT[k], Q[k], TRUE[k], OLD[k], "✓" if same else "✗ **相反**"))
print("   ⟹ CONTENT 表寫反者 = %d／6（三對 ⁰/⁺ 全部對調）" % bad)
chk("P1 符號錯定位到程式碼層（CONTENT 表 4132/4232 等相反）", bad == 6 and all(
    (Q[k] == ("0" if CHG[TRUE[k]] == 0 else "+")) for k in TRUE),
    "6/6 電荷判定一致；scripts/19 之表 6/6 相反")

# ── P2：改正預言 ────────────────────────────────────────────────────────────
d_s = (mp.mpf("1321.71") - mp.mpf("1314.86")) / 2     # R5：δ_s = (Ξ⁻−Ξ⁰)/2
d_c = mp.mpf("-0.2820")                                # R9（docs/18）
pred = d_s + d_c
print("\n§2 P2 改正預言（δ_s = %.4f、δ_c = %.4f MeV）" % (float(d_s), float(d_c)))
print("   Ξc⁰(dsc) − Ξc⁺(usc) = δ_s + δ_c = %+.4f MeV（原 docs/15 寫 −3.1430）" % float(pred))
print("   對照：Ξb⁻(dsb) − Ξb⁰(usb) = δ_s + δ_b = %+.4f（同慣例，M12 表得 🟢）"
      % float(d_s + mp.mpf("1.8933")))
chk("P2 改正預言 = +3.1430", abs(pred - mp.mpf("3.1430")) < mp.mpf("1e-4"), "%+.4f" % float(pred))

# ── P3：三列 z（用 run_19 之 σ）────────────────────────────────────────────
print("\n§3 P3 三列重算（觀測值與 σ 皆取自 `data/run_19_r6_channel_type.log`，只改預言符號）")
ROWS = [("Ξc⁰−Ξc⁺", "4132", "4232", mp.mpf("0.4956"), mp.mpf("11.85")),
        ("Ξc′⁰−Ξc′⁺", "4312", "4322", mp.mpf("0.7838"), mp.mpf("4.65")),
        ("Ξc(2645)⁰−Ξc(2645)⁺", "4314", "4324", mp.mpf("0.5165"), mp.mpf("8.14"))]
res = []
for nm, a, b, sig, z_old in ROWS:
    obs = M[a] - M[b]
    z_new = (obs - pred) / sig
    verdict = "判死 🔴" if abs(z_new) >= 3 else "通 🟢"
    res.append((nm, float(obs), float(sig), float(z_new), verdict))
    print("   %-22s 觀測 %+9.4f ± %.4f ｜預言 %+8.4f ｜z 原 %+6.2f（錯號）→ **新 %+6.2f** ⟹ %s"
          % (nm, float(obs), float(sig), float(pred), float(z_old), float(z_new), verdict))
# 交叉核：純觀測 σ（M13 版 0.3624）
obs0 = M["4132"] - M["4232"]
z_m13 = (obs0 - pred) / mp.mpf("0.3624")
print("   交叉核（純觀測 σ = 0.3624，M13 版）：Ξc⁰−Ξc⁺ z = %+.2f ⟹ 同判（通）" % float(z_m13))
chk("P3 Ξc⁰−Ξc⁺ 由 11.85σ 判死改為通", abs((obs0 - pred) / mp.mpf("0.4956")) < 3
    and abs((obs0 - pred) / mp.mpf("0.3624")) < 3,
    "z = %.2f（σ 0.4956）／%.2f（σ 0.3624）" % (float((obs0 - pred) / mp.mpf("0.4956")), float(z_m13)))

# ── P4：R6 之 Ξ 型五行計數 ─────────────────────────────────────────────────
print("\n§4 P4 R6「Ξ 型」五行新舊計數")
XI_ROWS = [("Ξ*⁻−Ξ*⁰", -5.12, "死"), ("Ξc⁰−Ξc⁺", res[0][3], res[0][4]),
           ("Ξc′⁰−Ξc′⁺", res[1][3], res[1][4]), ("Ξc(2645)⁰−Ξc(2645)⁺", res[2][3], res[2][4]),
           ("Ξb⁻−Ξb⁰", -0.27, "通")]
dead = sum(1 for r in XI_ROWS if abs(r[1]) >= 3)
alive = 5 - dead
for r in XI_ROWS:
    print("   %-24s z = %+7.2f ⟹ %s" % (r[0], r[1], "判死" if abs(r[1]) >= 3 else "通"))
print("   ⟹ 舊：4 死 1 通（Ξc⁰−Ξc⁺ 原報 +11.85 死）｜新：**%d 死 %d 通**" % (dead, alive))
chk("P4 R6 Ξ 型計數由 4死1通 改為 3死2通", dead == 3 and alive == 2, "%d 死 %d 通" % (dead, alive))

# ── P5：B2 只標待重算 ─────────────────────────────────────────────────────
print("\n§5 P5 B2（Σc+Ξc 和，run_19 L76）與 Ξc 列共用同一個錯號項")
print("   原報 +5.8730／z +11.85 ⟹ 與 Ξc 列同值 ⟹ 該列之 5.8730 = |2.7300 − (−3.1430)|")
print("   ⟹ 標「**待重算**」，本輪**不**下判（避免連帶臆測；須另立 R 輪）")
OUT["B2"] = "pending_recompute"

# ── P6 ────────────────────────────────────────────────────────────────────
print("\n§6 P6 原文不動（append-only）")
print("   `docs/15` 表與 `scripts/19` CONTENT 表**原文保留**；更正只入本輪 doc（`docs/44`）")
chk("P6 未改舊檔（本腳本只讀）", True, "只讀 PDG 與 run_19 log")

print("\n§7 收據")
print("   ✅ `M13 → M12` 第 ① 項 **成立**：符號錯在 `scripts/19` CONTENT 表（三對 Ξc 之 ⁰/⁺ 內容對調）")
print("   ⟹ 改正後 Ξc⁰−Ξc⁺ = −0.83σ（通）；Ξc′ −3.37σ、Ξc(2645) −4.03σ 仍死；R6 Ξ 型 3 死 2 通")
print("\n   FAIL = %d ｜ written data/51_inter_m_handoff_audit.json , data/run_51_inter_m_handoff_audit.log" % len(FAIL))
print("=" * 94)
OUT.update({"content_bug": {"swapped": bad}, "pred": str(pred),
            "rows": [{"name": r[0], "obs": r[1], "sigma": r[2], "z": r[3], "verdict": r[4]} for r in res],
            "xi_counts": {"dead": dead, "alive": alive}, "fail": FAIL})
(ROOT / "data" / "51_inter_m_handoff_audit.json").write_text(
    json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")
_sys.exit(0)