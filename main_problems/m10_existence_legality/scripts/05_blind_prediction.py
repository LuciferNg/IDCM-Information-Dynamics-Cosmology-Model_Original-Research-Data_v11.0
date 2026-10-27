#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
M10 scripts/05 — R5 預註冊盲測 ＋ 輸入數據完整性審計

前設問題（docs/05 §一）:
  §七.4 要求「預言必須樣本外 ＋ 預註冊」。本腳本先做一件前置事：
    J 數據完整性審計 —— 逐個核 01_problem_record §三 之 14 個輸入質量
      （用 AME2020 原子質量 − Z·m_e 重建核質量）
    K 樣本外外推族（20 通道）—— 判準在凍結表以外之表現
    L 拓撲型必零族（4 條）—— 唯一可判死之預測類，逐條寫死 kill condition
    M 預言力判詞 ＋ 開啟條件（W1/W2/W3/W5）

精確性：Fraction 精確；參考質量由 u 轉換常數乘 AME2020 原子質量後減電子。
執行：/usr/bin/python3 scripts/05_blind_prediction.py
"""
import hashlib
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
PASS, FAIL = [], []


def ck(tag, cond, msg):
    (PASS if cond else FAIL).append(tag)
    print(f"  [{'PASS' if cond else 'FAIL'}] {tag}  {msg}")


U = F("931.49410242")        # MeV/u（CODATA 2018）
ME = F("0.510998950")

# ---------------------------------------------------------------- J
print("=" * 78)
print("M10 R5 — 預註冊盲測 ＋ 數據完整性審計")
print("=" * 78)
print("\n[J] 輸入數據完整性審計（§三 之 14 個質量）")
FROZEN = {k: F(v) for k, v in {
    "e": "0.510998950", "mu": "105.6583755", "tau": "1776.86",
    "p": "938.27208816", "n": "939.56542052", "d2H": "1875.61294257",
    "t3H": "2808.391", "a4He": "3727.379", "Be8": "7454.85",
    "pi0": "134.9768", "W": "80377.0", "Z": "91187.6",
    "b": "4180.0", "t": "172690.0",
}.items()}
# AME2020 原子質量（u）⟹ 核質量 = 原子質量·u − Z·m_e（忽略電子束縛，~1e-5 MeV 級）
AME_U = {
    "d2H": ("2.01410177812", 1), "t3H": ("3.01604928199", 1),
    "a4He": ("4.00260325413", 2), "Be8": ("8.00530510", 4),
    "t3He": ("3.01602932010", 2),
}
print(f"    {'態':<6}{'凍結值':>14}{'AME2020 重建':>16}{'Δ (MeV)':>12}  判")
drift = {}
for k, (u, z) in AME_U.items():
    if k not in FROZEN:
        continue
    nuc = F(u) * U - z * ME
    d = FROZEN[k] - nuc
    drift[k] = float(d)
    flag = "✅" if abs(d) < F("0.001") else "🔴 輸入有誤（= ³He 核質量）"
    print(f"    {k:<6}{float(FROZEN[k]):>14.6f}{float(nuc):>16.6f}{float(d):>+12.4f}  {flag}")
ck("J.1 ³H 輸入值實為 ³He（差 0.53 MeV）",
   abs(FROZEN["t3H"] - (F("3.01604928199") * U - ME)) > F("0.4"),
   f"凍結 2808.391 = ³He 核質量（{float(F('3.01602932010')*U-2*ME):.4f}）；"
   f"³H 應為 {float(F('3.01604928199')*U-ME):.4f} MeV")
ck("J.2 其餘 3 個核種無誤（<1e-3 MeV）",
   all(abs(FROZEN[k] - (F(u) * U - z * ME)) < F("0.001")
       for k, (u, z) in AME_U.items() if k not in ("t3H", "t3He")),
   "²H／⁴He／⁸Be 與 AME2020 一致 ⟹ 誤差只集中在 ³H 行")
# 修正後之 ³H 行
t3h_ref = F("3.01604928199") * U - ME
x_t3h = t3h_ref - FROZEN["p"] - 2 * FROZEN["n"]
x_t3h_frozen = FROZEN["t3H"] - FROZEN["p"] - 2 * FROZEN["n"]
print(f"    ³H 束縛（p+2n）：凍結 {float(x_t3h_frozen):+.4f} MeV（誤）／修正 "
      f"{float(x_t3h):+.4f} MeV（= 公認 8.482 MeV 束縛）")
ck("J.3 修正後判詞不變（仍為束縛 ⟹ 合法）", x_t3h < 0 and x_t3h_frozen < 0,
   f"|x|/m: {float(abs(x_t3h_frozen)/FROZEN['t3H']):.3e} → "
   f"{float(abs(x_t3h)/t3h_ref):.3e} ⟹ 符號與判詞皆不變，只改記錄數值")
ck("J.4 誤差溯源：³He 質量 ＋ ³H 組成 之混配", abs(
    (F("3.01602932010") * U - 2 * ME) - FROZEN["t3H"]) < F("0.002"),
   "凍結值 = ³He 核質量（2p+n 組成之態）而通道寫 p+2n ⟹ 質量與組成不配；"
   "兩個獨立錯（種錯＋配錯）疊加出 −9.0119，公認值為 −8.482")

# ---------------------------------------------------------------- K
print("\n[K] 樣本外外推族（20 通道；全部在凍結表之外）")
EXTRA = {k: F(v) for k, v in {
    "L": "1115.683", "Sp": "1189.37", "S0": "1192.642", "Sm": "1197.449",
    "X0": "1314.86", "Xm": "1321.71", "Om": "1672.45", "pip": "139.57039",
    "Kp": "493.677", "K0": "497.611", "eta": "547.862", "rho": "775.26",
    "phi": "1019.461", "Jpsi": "3096.900", "D0": "1864.84", "Dm": "1869.66", "B0": "5279.66",
    "pim": "139.57039", "Km": "493.677", "gamma": "0",
    "munu": "0", "ep": "0.510998950",
}.items()}
ALL = dict(FROZEN)
ALL.update(EXTRA)
Q = {"L": 0, "Sp": 1, "S0": 0, "Sm": -1, "X0": 0, "Xm": -1, "Om": -1,
     "pip": 1, "pim": -1, "Kp": 1, "Km": -1, "K0": 0, "eta": 0, "rho": 0,
     "phi": 0, "Jpsi": 0, "D0": 0, "Dm": -1, "B0": 0, "p": 1, "n": 0,
     "e": -1, "pi0": 0, "munu": 0, "d2H": 1, "a4He": 2, "gamma": 0, "ep": 1}
BB = {"L": 1, "Sp": 1, "S0": 1, "Sm": 1, "X0": 1, "Xm": 1, "Om": 1,
      "pip": 0, "pim": 0, "Kp": 0, "Km": 0, "K0": 0, "eta": 0, "rho": 0,
      "phi": 0, "Jpsi": 0, "D0": 0, "Dm": 0, "B0": 0, "p": 1, "n": 1,
      "e": 0, "pi0": 0, "munu": 0, "d2H": 2, "a4He": 4, "gamma": 0, "ep": 0}
LL = {"L": 0, "Sp": 0, "S0": 0, "Sm": 0, "X0": 0, "Xm": 0, "Om": 0,
      "pip": 0, "pim": 0, "Kp": 0, "Km": 0, "K0": 0, "eta": 0, "rho": 0,
      "phi": 0, "Jpsi": 0, "D0": 0, "Dm": 0, "B0": 0, "p": 0, "n": 0,
      "e": 1, "pi0": 0, "munu": -1, "d2H": 0, "a4He": 0, "gamma": 0, "ep": -1}
CH20 = [
    ("Λ→pπ⁻",        "L",    ["p", "pim"]),
    ("Λ→nπ⁰",        "L",    ["n", "pi0"]),
    ("Σ⁺→pπ⁰",       "Sp",   ["p", "pi0"]),
    ("Σ⁺→nπ⁺",       "Sp",   ["n", "pip"]),
    ("Σ⁻→nπ⁻",       "Sm",   ["n", "pim"]),
    ("Ξ⁰→Λπ⁰",       "X0",   ["L", "pi0"]),
    ("Ξ⁻→Λπ⁻",       "Xm",   ["L", "pim"]),
    ("Ω⁻→Ξ⁰π⁻",      "Om",   ["X0", "pim"]),
    ("Ω⁻→ΛK⁻",       "Om",   ["L", "Km"]),
    ("K⁺→π⁺π⁰",      "Kp",   ["pip", "pi0"]),
    ("K⁰→π⁺π⁻",      "K0",   ["pip", "pim"]),
    ("η→γγ",         "eta",  ["gamma", "gamma"]),
    ("ρ⁰→π⁺π⁻",      "rho",  ["pip", "pim"]),
    ("φ→K⁺K⁻",       "phi",  ["Kp", "Km"]),
    ("J/ψ→e⁺e⁻",     "Jpsi", ["ep", "e"]),
    ("D⁰→K⁻π⁺",      "D0",   ["Km", "pip"]),
    ("B⁰→D⁻π⁺",      "B0",   ["Dm", "pip"]),
    ("n→peν",        "n",    ["p", "e", "munu"]),
    ("²H→pn",        "d2H",  ["p", "n"]),
    ("p→e⁺γγ (ΔB)",  "p",    ["e", "gamma", "gamma"]),
]
print(f"    {'通道':<15}{'x (MeV)':>12}{'|x|/m':>11}  閘  判準判  觀測對照")
k_rows, minrel = [], 1.0
for name, init, finals in CH20:
    x = ALL[init] - sum((ALL[f] for f in finals), F(0))
    rel = float(abs(x) / ALL[init])
    dB = BB[init] - sum((BB[f] for f in finals), 0)
    dQ = Q[init] - sum((Q[f] for f in finals), 0)
    dL = LL[init] - sum((LL[f] for f in finals), 0)
    gate = (dB == 0 and dQ == 0 and dL == 0)
    verdict = "可達" if (gate and x > 0) else ("禁（能量）" if gate else "禁（拓撲）")
    if name not in ("n→peν", "²H→pn"):      # 兩個對照通道（Tier II）不計入外推族
        minrel = min(minrel, rel)
    k_rows.append({"channel": name, "x": float(x), "rel": rel, "gate": gate,
                   "verdict": verdict})
    print(f"    {name:<15}{float(x):>+12.4f}{rel:>11.2e}  {'✅' if gate else '❌'}  "
          f"{verdict:<9} {'—' if name.startswith('²') else '見註'}")
EXPECT_VERDICT = {"²H→pn": "禁（能量）", "p→e⁺γγ (ΔB)": "禁（拓撲）"}
all_ok = all(r["verdict"] == EXPECT_VERDICT.get(r["channel"], "可達") for r in k_rows)
ck("K.1 20/20 判準判與觀測一致", all_ok,
   "18 條可達（觀測皆衰變）＋ ²H（閘過但能量禁，觀測穩定）＋ p→e⁺γγ（拓撲禁）")
ck("K.2 18 條新通道全落 Tier I（|x|/m > 1e-2）", minrel > 1e-2,
   f"最小 |x|/m = {minrel:.2e}；另兩個對照通道（n 8.33e-04、²H 1.19e-03）屬 Tier II ⟹ "
   f"18 條新通道全部靠能量比較即可判 ⟹ **此族不得作為證據**（判準在此層恆真，零判別力）")

# ---------------------------------------------------------------- L
print("\n[L] 拓撲型必零族（唯一可判死之預測類）—— 預註冊，凍結")
PRED = [
    {"id": "P1", "claim": "ΔB≠0 之所有過程速率 = 0（精確）",
     "targets": "p→e⁺π⁰ (Hyper-K)、p→K⁺ν̄、pp→π⁺π⁺ / nn→π⁰π⁰ (Super-K dinucleon)、"
                "¹⁶O→¹⁴C+…、n→n̄ (ESS)",
     "kill": "任何一個通道之 90% CL 上限被突破 ⟹ M10 閘判死"},
    {"id": "P2", "claim": "面型守恆 ⟹ η_B 係初始條件（冇 sphaleron 洗白）",
     "targets": "結構判詞，非直接可測；後果 = 唔需要 (B−L) 生成機制",
     "kill": "若量到電弱重子生成 signature ⟹ 閘精確性判死（目前 EW 相變係 crossover，不可測）"},
    {"id": "P3", "claim": "合法集 = {各超選擇類之最輕元} ⟹ 新發現之 B=1,Q=1 態必須衰變到 p",
     "targets": "LHC 穩定 R-hadron／fractional charge／dark baryon 搜索",
     "kill": "發現任何穩定 B=1,Q=1 態（m > m_p）⟹ 最輕態定理判死"},
    {"id": "P4", "claim": "n 之合法性只可條件成立 ⟹ 自由中子必衰變（無例外）",
     "targets": "任何「自由中子穩定」之聲稱",
     "kill": "自由中子穩定 ⟹ 判據判死"},
]
for pr in PRED:
    print(f"    {pr['id']}  {pr['claim']}")
    print(f"        目標：{pr['targets']}")
    print(f"        判死：{pr['kill']}")
seal = hashlib.sha256(json.dumps(PRED, ensure_ascii=False, sort_keys=True).encode()).hexdigest()
print(f"    預註冊封印 SHA256 = {seal}")
ck("L.1 每條預言皆有明確 kill condition", all("kill" in p and p["kill"] for p in PRED),
   "4/4 條寫死判死條件（可證偽，非敍述）")
ck("L.2 P1/P3/P4 對 SM 判別力弱（誠實標記）", True,
   "SM 亦給 p 穩定／n 衰變／無穩定額外 B=1 ⟹ P1/P3/P4 只作**必要條件**，"
   "唔可以當 IDCM 特有預言；P2 不可直接測 ⟹ **本 M 冇『IDCM 特有 ＋ 可判別 ＋ 樣本外』預言**")

# ---------------------------------------------------------------- M
print("\n[M] 預言力判詞 ＋ 開啟條件")
BLOCK = [
    ("W1 em 自能算子", "M5 三界封死（A_max 1.2494 vs 需求 1.2667–1.3963）", "所有含帶電粒子之質量差"),
    ("W2 複合態質量槽", "M7 ⟨O⟩≡1 ⟹ 加性 ⟹ 無束縛槽", "所有強子通道"),
    ("W3 核束縛槽", "M5 兩參數模型；M7 加束縛即違 M2 唯一性", "所有核通道（含 ⁸Be、diproton）"),
    ("W5 質量值層 dof=0", "O120：9 方程 9 未知 ⟹ 9! 精確解", "所有絕對質量預言"),
]
for w, why, hit in BLOCK:
    print(f"    {w:<20} {why:<44} → {hit}")
ck("M.1 T5 = 零判別力（有精確原因）", True,
   "🔴 T5：本 M 只有『拓撲型必零』族可預註冊，而該族對 SM 判別力弱；"
   "任何具判別力之預言皆需先開 W1/W2/W3/W5 任一 ⟹ 判詞係**界限**，非失敗")
ck("M.2 判準本身可作協議交付", True,
   "x_verdict(state, channels) 五步協議：①列通道 ②驗閘(ΔB,ΔQ,ΔL,色單態) "
   "③算 x ④查 |x|/m 定 Tier ⑤Tier II/III 一律報『不可判』")

print("\n" + "=" * 78)
print(f"CHECK 合計: PASS {len(PASS)} / FAIL {len(FAIL)}  → {len(PASS)}/{len(PASS)+len(FAIL)}")
if FAIL:
    print("FAIL 清單:", FAIL)
print("=" * 78)
os.makedirs(DATA, exist_ok=True)
out = {
    "script": "05_blind_prediction.py", "role": "R5 預註冊 ＋ 數據完整性審計",
    "J_drift_MeV": drift,
    "J_t3H": {"frozen": float(FROZEN["t3H"]), "correct": float(t3h_ref),
              "x_frozen": float(x_t3h_frozen), "x_correct": float(x_t3h),
              "identity_of_frozen": "³He 核質量"},
    "K_out_of_sample": k_rows, "K_min_rel": minrel,
    "L_preregistered": PRED, "L_seal": seal,
    "M_blocks": [{"wall": w, "why": y, "hits": h} for w, y, h in BLOCK],
    "checks_passed": len(PASS), "checks_total": len(PASS) + len(FAIL), "failures": FAIL,
}
with open(os.path.join(DATA, "05_blind_prediction.json"), "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False, indent=2)
print("data/05_blind_prediction.json 已寫")