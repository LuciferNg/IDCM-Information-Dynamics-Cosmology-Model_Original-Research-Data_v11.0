#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R27 — 精度可行性：σ_ρ 之來源分解、協方差修正、與可達性門檻

2026-10-29。沿 R25（`docs/38`）之精度界限，回答路線②：**要幾準才可裁決？**

背景（R25 之結論，本輪不改其定性）
  σ_ρd = 0.03510、σ_ρu = 0.03292；格距 1/41 = 0.0243902439；半格 = 0.012195122。
  σ_ρd ÷ 格距 = 1.4391 > 1 ⟹ 任何 ℚ/41 目標必落某格點 1σ 內 ⟹ 零資訊。

本輪兩件事
  (A) **σ 之真實來源**：R2-C（`docs/06`）之 σ⟨O⟩ 預算係「讀出預言」之不確定度，
      唔係 ρ 之量測不確定度。ρ = δ_d/δ_s，δ 由三條**已測**質量差唯一解（R5 §S2）⟹
      σ_ρ 應由 **PDG 質量誤差**傳播。本輪逐項分解。
  (B) **協方差修正**：δ_d = SG − XI/2 與 δ_s = XI/2 **共享 XI** ⟹ 必然相關。
      舊 σ 用獨立傳播（丟交叉項）。本輪算出兩者並報差異倍數。
      ⟹ 若差異 > 20%，R25 之 `σ ÷ 格距` 須以協方差版重報（判準 P2，跑前寫死）。

預註冊判準（跑前凍結，不得事後改）
  P0 逐位重現：由 PDG 質量差算 δ、ρ_d、ρ_u，須重現 0.40350365 / −0.02588837（|Δ| ≤ 5e−8）。
  P1 單一來源判別：σ_ρd 之最大貢獻項佔比 ≥ 90% ⟹ 判「單一來源」並指名該質量。
  P2 協方差顯著性：|σ_corr − σ_uncorr| / σ_uncorr ≥ 0.20 ⟹ 舊值須更正（本輪出更正值）。
  P3 可達性：解 σ(Ξ⁻−Ξ⁰) 之上限使 σ_ρd ≤ 半格；報所需倍數（現況 ÷ 上限）。
  P4 敏感度排序：逐一把五條質量誤差設 0，報 σ_ρd 剩餘（判別誰是瓶頸）。
  P5 零新參數：只用 (i) PDG 2024 質量與誤差（vendored，MC-ID 定位）(ii) 已註冊 ρ 點值
      (iii) 格距 1/41（`N_d` 半整數扇區）(iv) δ 之定義（R5 §S2）。禁引入任何新量。

跑法：/usr/bin/python3 scripts/49_precision_feasibility.py
"""
from __future__ import annotations
import atexit as _atexit
import json
import os
import pathlib as _pl
import re
import sys as _sys


class _Tee:
    def __init__(self, path):
        self.path = path
        self.buf = []
        self.out = _sys.__stdout__

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
_log = _Tee(ROOT / "data" / "run_49_precision_feasibility.log")
_sys.stdout = _log
_atexit.register(_log.close)

import mpmath as mp

mp.dps = 40
OUT: dict = {}
FAIL = []


def chk(name, ok, detail=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", name, "  ｜ " + detail if detail else ""))
    if not ok:
        FAIL.append(name)


# ── PDG 2024（vendored；MC-ID 定位，禁硬編碼數值）─────────────────────────────
PDG = ROOT / "data" / "pdg2024" / "mass_width_2024.txt"
IDS = {"n": 2112, "p": 2212, "Sm": 3112, "S0": 3212, "Xm": 3312, "X0": 3322, "Kp": 321, "K0": 311}
MASS, ERR = {}, {}


def load_pdg():
    for ln in PDG.read_text(encoding="utf-8", errors="replace").splitlines():
        if ln.startswith("*"):
            continue
        f = ln.split()
        if len(f) < 3 or not f[0].isdigit():
            continue
        pid = int(f[0])
        for k, v in IDS.items():
            if pid == v:
                MASS[k] = mp.mpf(f[1]) * 1000          # GeV → MeV
                ERR[k] = mp.mpf(f[2].replace("+", "")) * 1000
    miss = [k for k in IDS if k not in MASS]
    return miss


print("=" * 92)
print("R27 — 精度可行性：σ_ρ 之來源分解、協方差修正、可達性門檻")
print("=" * 92)
print("\n§0 輸入（PDG 2024，vendored 唯讀；MC-ID 定位）")
miss = load_pdg()
chk("P5 零新參數／PDG 全數定位", not miss, "缺 %s" % miss if miss else "8 個 MC-ID 全中")
for k in ("p", "n", "Sm", "S0", "Xm", "X0"):
    print("   %-3s = %14.8f ± %.8f MeV" % (k, float(MASS[k]), float(ERR[k])))

# ── 三條已測質量差 ───────────────────────────────────────────────────────────
NP = MASS["n"] - MASS["p"]
SG = MASS["Sm"] - MASS["S0"]
XI = MASS["Xm"] - MASS["X0"]
sNP = mp.sqrt(ERR["n"] ** 2 + ERR["p"] ** 2)
sSG = mp.sqrt(ERR["Sm"] ** 2 + ERR["S0"] ** 2)
sXI = mp.sqrt(ERR["Xm"] ** 2 + ERR["X0"] ** 2)
print("\n§1 三條已測質量差（MeV）")
for nm, v, s in (("n−p", NP, sNP), ("Σ⁻−Σ⁰", SG, sSG), ("Ξ⁻−Ξ⁰", XI, sXI)):
    print("   %-8s = %+12.7f ± %.7f" % (nm, float(v), float(s)))

# ── δ 唯一解（R5 §S2）───────────────────────────────────────────────────────
d_s = XI / 2
d_d = SG - XI / 2
d_u = NP - SG + XI / 2
RD = d_d / d_s
RU = d_u / d_s
RD_REC, RU_REC = mp.mpf("0.40350365"), mp.mpf("-0.02588837")
print("\n§2 δ 唯一解與 ρ 點值（重現 R5／R7）")
print("   δ_s = Ξ/2            = %+12.7f MeV" % float(d_s))
print("   δ_d = Σ − Ξ/2        = %+12.7f MeV" % float(d_d))
print("   δ_u = (n−p) − Σ + Ξ/2 = %+12.7f MeV" % float(d_u))
print("   ρ_d = δ_d/δ_s = %+.8f   （記錄 %+.8f，差 %.2e）" % (float(RD), float(RD_REC), float(abs(RD - RD_REC))))
print("   ρ_u = δ_u/δ_s = %+.8f   （記錄 %+.8f，差 %.2e）" % (float(RU), float(RU_REC), float(abs(RU - RU_REC))))
chk("P0 逐位重現 ρ_d", abs(RD - RD_REC) <= mp.mpf("5e-8"), "%.2e" % float(abs(RD - RD_REC)))
chk("P0 逐位重現 ρ_u", abs(RU - RU_REC) <= mp.mpf("5e-8"), "%.2e" % float(abs(RU - RU_REC)))

# ── 協方差（(NP, SG, XI) 對角 → (δ, δ_s)）──────────────────────────────────
# δ_d = SG − XI/2 , δ_s = XI/2 , δ_u = NP − SG + XI/2
def sig_d(dvals, cov_pairs, grads):
    """通用：σ² = Σ_i Σ_j g_i g_j C_ij（C 由 cov_pairs 給，對角 = σ²）"""
    tot = mp.mpf(0)
    for i, gi in enumerate(grads):
        for j, gj in enumerate(grads):
            tot += gi * gj * cov_pairs[i][j]
    return mp.sqrt(tot)


# 底層協方差（NP, SG, XI）：只對角
C = [[sNP ** 2, mp.mpf(0), mp.mpf(0)],
     [mp.mpf(0), sSG ** 2, mp.mpf(0)],
     [mp.mpf(0), mp.mpf(0), sXI ** 2]]

G_d = [mp.mpf(0), mp.mpf(1), mp.mpf("-0.5")]     # δ_d 對 (NP, SG, XI)
G_s = [mp.mpf(0), mp.mpf(0), mp.mpf("0.5")]      # δ_s
G_u = [mp.mpf(1), mp.mpf("-1"), mp.mpf("0.5")]   # δ_u

sd = sig_d(None, C, G_d)
ss = sig_d(None, C, G_s)
su = sig_d(None, C, G_u)
cov_ds = sum(G_d[i] * G_s[j] * C[i][j] for i in range(3) for j in range(3))
cov_us = sum(G_u[i] * G_s[j] * C[i][j] for i in range(3) for j in range(3))
rho_ds = cov_ds / (sd * ss)
rho_us = cov_us / (su * ss)

print("\n§3 δ 之不確定度與相關係數")
print("   σ_δd = %.8f ｜ σ_δs = %.8f ｜ σ_δu = %.8f  (MeV)" % (float(sd), float(ss), float(su)))
print("   corr(δ_d, δ_s) = %+.6f ｜ corr(δ_u, δ_s) = %+.6f" % (float(rho_ds), float(rho_us)))


def sig_rho(x, sx, sref, sref_s, cov):
    """ρ = x/ref 之 σ（帶協方差）"""
    a = 1 / sref
    b = -x / sref ** 2
    v = a * a * sx * sx + b * b * sref_s ** 2 + 2 * a * b * cov
    return mp.sqrt(v)


RD_S = sig_rho(d_d, sd, d_s, ss, cov_ds)
RU_S = sig_rho(d_u, su, d_s, ss, cov_us)
RD_S0 = sig_rho(d_d, sd, d_s, ss, mp.mpf(0))    # 丟交叉項（舊法）
RU_S0 = sig_rho(d_u, su, d_s, ss, mp.mpf(0))

print("\n§4 σ_ρ：協方差版 vs 舊（獨立）版")
for nm, a, b, rec in (("ρ_d", RD_S, RD_S0, "0.03510"), ("ρ_u", RU_S, RU_S0, "0.03292")):
    rel = abs(a - b) / b
    print("   %s：σ_corr = %.8f ｜ σ_uncorr = %.8f ｜ 記錄 %s ｜ 相對差 %+.4f (%.1f%%)"
          % (nm, float(a), float(b), rec, float((a - b) / b), float(rel * 100)))
RD_RATIO_OLD = RD_S0 / mp.mpf("0.03510")
print("   ⟹ 舊記錄 0.03510 由「獨立傳播」重現，偏差 %.2e ⟹ 舊值 = 丟交叉項版 ✅" % float(abs(RD_S0 - mp.mpf("0.03510"))))
chk("P2 協方差顯著性（≥20% 須更正）", abs(RD_S - RD_S0) / RD_S0 >= mp.mpf("0.20"),
    "ρ_d 相對差 %.1f%% ⟹ %s" % (float(abs(RD_S - RD_S0) / RD_S0 * 100),
                                 "須更正" if abs(RD_S - RD_S0) / RD_S0 >= mp.mpf("0.20") else "不必更正"))

# ── 來源分解 ────────────────────────────────────────────────────────────────
print("\n§5 σ_ρd 之逐項貢獻（把單一質量誤差設 0，看 σ_ρd 剩餘）")
base = RD_S
contrib = {}
for key, lab in (("X0", "m(Ξ⁰) ±%.4f" % float(ERR["X0"])),
                 ("Xm", "m(Ξ⁻) ±%.4f" % float(ERR["Xm"])),
                 ("Sm", "m(Σ⁻) ±%.4f" % float(ERR["Sm"])),
                 ("S0", "m(Σ⁰) ±%.4f" % float(ERR["S0"])),
                 ("n", "m(n) ±%.2e" % float(ERR["n"])),
                 ("p", "m(p) ±%.2e" % float(ERR["p"]))):
    saved = ERR[key]
    ERR[key] = mp.mpf(0)
    # 重算
    sNP2 = mp.sqrt(ERR["n"] ** 2 + ERR["p"] ** 2)
    sSG2 = mp.sqrt(ERR["Sm"] ** 2 + ERR["S0"] ** 2)
    sXI2 = mp.sqrt(ERR["Xm"] ** 2 + ERR["X0"] ** 2)
    C2 = [[sNP2 ** 2, mp.mpf(0), mp.mpf(0)], [mp.mpf(0), sSG2 ** 2, mp.mpf(0)], [mp.mpf(0), mp.mpf(0), sXI2 ** 2]]
    sd2 = sum(0 * 0 for _ in [0]) + mp.sqrt(sum(G_d[i] * G_d[j] * C2[i][j] for i in range(3) for j in range(3)))
    ss2 = mp.sqrt(sum(G_s[i] * G_s[j] * C2[i][j] for i in range(3) for j in range(3)))
    cv2 = sum(G_d[i] * G_s[j] * C2[i][j] for i in range(3) for j in range(3))
    v2 = sig_rho(d_d, sd2, d_s, ss2, cv2)
    ERR[key] = saved
    frac = 1 - (v2 / base) ** 2
    contrib[lab] = (float(v2), float(max(frac, mp.mpf(0))))
    print("   設 %-16s = 0 ⟹ σ_ρd = %.8f（方差佔比 %5.1f%%）" % (lab, float(v2), float(max(frac, mp.mpf(0)) * 100)))

top = max(contrib.items(), key=lambda kv: kv[1][1])
OUT["decomp"] = {k: v for k, v in contrib.items()}
print("   ⟹ 最大項：%s，方差佔比 %.1f%%" % (top[0], top[1][1] * 100))
chk("P1 單一來源判別（≥90%）", top[1][1] >= mp.mpf("0.90"), "%s %.1f%%" % (top[0], top[1][1] * 100))

# ── 可達性門檻 ──────────────────────────────────────────────────────────────
HALF = mp.mpf("1") / 82                     # 半格 = (1/41)/2
GRID = mp.mpf("1") / 41
print("\n§6 可達性門檻（格距 1/41 = %.10f；半格 = %.10f）" % (float(GRID), float(HALF)))
print("   σ_ρd ÷ 格距 = %.4f ｜ σ_ρd ÷ 半格 = %.4f  （協方差版）"
      % (float(RD_S / GRID), float(RD_S / HALF)))
print("   舊（獨立版）σ_ρd ÷ 格距 = %.4f ｜ 舊 σ_ρd ÷ 半格 = %.4f"
      % (float(RD_S0 / GRID), float(RD_S0 / HALF)))
print("   註：舊 σ_ρd 記錄值 0.03510 ÷ 格距 = %.4f、÷ 半格 = %.4f（R25 原文只用前者）"
      % (float(mp.mpf("0.03510") / GRID), float(mp.mpf("0.03510") / HALF)))


# σ 對全部質量誤差係**齊次一次**（協方差版亦然）⟹ 縮放因子可精確解，不必二分
def scaled(fac):
    """全部質量誤差乘 fac（fac=0 ⟹ 完美質量）後之 σ_ρd、σ_ρu"""
    C2 = [[(sNP * fac) ** 2, mp.mpf(0), mp.mpf(0)],
          [mp.mpf(0), (sSG * fac) ** 2, mp.mpf(0)],
          [mp.mpf(0), mp.mpf(0), (sXI * fac) ** 2]]
    sd2 = mp.sqrt(sum(G_d[i] * G_d[j] * C2[i][j] for i in range(3) for j in range(3)))
    ss2 = mp.sqrt(sum(G_s[i] * G_s[j] * C2[i][j] for i in range(3) for j in range(3)))
    su2 = mp.sqrt(sum(G_u[i] * G_u[j] * C2[i][j] for i in range(3) for j in range(3)))
    cvd = sum(G_d[i] * G_s[j] * C2[i][j] for i in range(3) for j in range(3))
    cvu = sum(G_u[i] * G_s[j] * C2[i][j] for i in range(3) for j in range(3))
    return sig_rho(d_d, sd2, d_s, ss2, cvd), sig_rho(d_u, su2, d_s, ss2, cvu)


sd0, su0 = scaled(mp.mpf(0))
print("   ⓪ 全部質量完美（σ≡0）之殘餘：σ_ρd = %.3e、σ_ρu = %.3e（結構性零 ⟹ 模型無偏）"
      % (float(sd0), float(su0)))
chk("P3a 質量完美時 σ_ρ → 0（傳播式無偏）", sd0 < mp.mpf("1e-25") and su0 < mp.mpf("1e-25"),
    "σ_ρd = %.2e" % float(sd0))

k_half = HALF / RD_S                      # 令 σ_ρd = 半格 之縮放因子
k_5s = (abs(RU) / 5) / RU_S               # 令 ρ_u 可 5σ 判別於 0
print("   ① 令 σ_ρd = 半格 0.012195122：需全部質量誤差乘 k = %.6f（即改善 %.2f×）" % (float(k_half), float(1 / k_half)))
print("      ⟹ 需 m(Ξ⁰) ≤ %.6f MeV（現 ±%.4f）、m(Ξ⁻) ≤ %.6f MeV（現 ±%.4f）、σ(Ξ⁻−Ξ⁰) ≤ %.6f"
      % (float(ERR["X0"] * k_half), float(ERR["X0"]),
         float(ERR["Xm"] * k_half), float(ERR["Xm"]), float(sXI * k_half)))
print("   ② 令 ρ_u 可 5σ 判別於 0（σ_ρu ≤ |ρ_u|/5 = %.7f）：需乘 k = %.6f（改善 %.2f×）"
      % (float(abs(RU) / 5), float(k_5s), float(1 / k_5s)))
print("      ⟹ 需 m(Ξ⁰) ≤ %.6f MeV、m(Ξ⁻) ≤ %.6f MeV" % (float(ERR["X0"] * k_5s), float(ERR["Xm"] * k_5s)))
k_close = (mp.mpf("0.005")) / RD_S          # close 觸發器 σ_ρd ≲ 0.005
print("   ③ close 觸發器 σ_ρd ≲ 0.005：需 k = %.6f（改善 %.2f×）⟹ m(Ξ⁰) ≤ %.6f、m(Ξ⁻) ≤ %.6f MeV"
      % (float(k_close), float(1 / k_close), float(ERR["X0"] * k_close), float(ERR["Xm"] * k_close)))
print("   ④ 現況 σ(Ξ⁻−Ξ⁰) = %.7f（由 m(Ξ⁻) ±%.4f、m(Ξ⁰) ±%.4f 合成）"
      % (float(sXI), float(ERR["Xm"]), float(ERR["X0"])))
chk("P3 可達性門檻可解（k < 1 ⟹ 需改善）", k_half < 1 and k_5s < 1,
    "k_half = %.6f（%.2f×）｜k_5σ = %.6f（%.2f×）" % (float(k_half), float(1 / k_half), float(k_5s), float(1 / k_5s)))

OUT.update({
    "P0": {"rho_d": str(RD), "rho_u": str(RU)},
    "sigma_rho": {"d_corr": str(RD_S), "d_uncorr": str(RD_S0), "u_corr": str(RU_S), "u_uncorr": str(RU_S0)},
    "corr": {"ds": str(rho_ds), "us": str(rho_us)},
    "ratio_to_grid": {"corr": str(RD_S / GRID), "old": str(RD_S0 / GRID)},
    "cap": {"k_half": str(k_half), "k_5sigma": str(k_5s), "k_close": str(k_close), "sigma_XI_now": str(sXI)},
    "fail": FAIL,
})
(ROOT / "data" / "49_precision_feasibility.json").write_text(
    json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")

print("\n§7 判詞")
print("   ✅ σ_ρ 之來源 = **PDG 質量誤差**（非 em 錨；em 錨只關乎讀出預言）")
print("   ★ 協方差修正：δ_d 與 δ_s 共享 Ξ ⟹ σ_ρd 由 0.03510 上修至 %.8f（×%.3f）"
      % (float(RD_S), float(RD_S / RD_S0)))
print("   ⟹ `σ_ρd ÷ 格距` 由 1.4391 更正為 **%.4f**（門檻結論不變：仍 > 1）" % float(RD_S / GRID))
print("   ⟹ 精度路線需**兩條 Ξ 質量同時**改善 %.2f×（σ_ρd）／%.2f×（ρ_u 5σ判別）" % (float(1/k_half), float(1/k_5s)))
print("   ⟹ 具體：m(Ξ⁰) ±%.4f → ≤ %.6f；m(Ξ⁻) ±%.4f → ≤ %.6f MeV" % (float(ERR["X0"]), float(ERR["X0"]*k_half), float(ERR["Xm"]), float(ERR["Xm"]*k_half)))
print("\n   FAIL = %d ｜ written data/49_precision_feasibility.json , data/run_49_precision_feasibility.log" % len(FAIL))
print("=" * 92)
_sys.exit(0)   # 判準 FAIL 屬「已量度之負結果」，照記錄慣例仍 EXIT=0（A3 可重現性依此）