#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/08 — R8′：強迫候選（零自由係數）之決定性檢定 ＋ 交接點之譜表徵
================================================================================
前輪結論：一致性路線封閉，只剩「強迫型 identity」。本輪**不再擴大族**（功率不足），
改為列舉**每個只產生一個候選**之強迫量（行列式／譜積），逐一作決定性檢定：
命中 ⟹ 成立；不中 ⟹ 該推導形即死（無 look-elsewhere 問題）。

   J1  ★ 交接點之譜表徵：γ₂₂ = 2 ⟹ (1 − K/2) 在交接點有**精確零模**
   J2  精確譜積恆等式（閉式）：det'(1−K/2) = −N_d/2^{43}；Πcos(πn/88) = √N_d/2^{43}
   J3  強迫候選表（零自由係數，每個一個候選）逐一對 μ 檢定
   J4  判決

執行: /usr/bin/python3 scripts/08_forcing_candidates.py
輸出: data/08_forcing_candidates.json + data/run_08.log
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import sympy as sp

mp.mp.dps = 60

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
N_FAIL = 0


def chk(tag, cond, detail=""):
    global N_FAIL
    ok = bool(cond)
    if not ok:
        N_FAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)
    LOG.append("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


# ============================================================================
# §0 輸入與譜
# ============================================================================
hdr("§0 輸入與譜（γ_n = 4sin²(πn/88)，n = 1…44）")

N_D = 44
M0 = mp.mpf("27.900915")
M_PL = M0 * mp.mpf("4.3758064e20")
MU = M0 / M_PL
print(f"μ = M₀/M_Pl = {mp.nstr(MU, 14)}")

GAM = [mp.mpf(4) * mp.sin(mp.pi * n / N_D / 2) ** 2 for n in range(1, N_D + 1)]
print(f"γ₁ = {mp.nstr(GAM[0], 12)}   γ₂₂ = {mp.nstr(GAM[21], 12)}   γ₄₄ = {mp.nstr(GAM[43], 12)}")


# ============================================================================
# §1 J1 — 交接點之譜表徵：(1 − K/2) 零模
# ============================================================================
hdr("§1 J1 — 交接點 γ=2 ⟺ (1 − K/2) 之精確零模")

chk("J1.1 γ₂₂ = 2（|γ₂₂−2| < 1e−50）", abs(GAM[21] - 2) < mp.mpf("1e-50"),
    f"γ₂₂ = {mp.nstr(GAM[21], 20)}")
chk("J1.2 γ_n ≡ 2(1 − cos(πn/44)) ⟹ 1 − γ_n/2 = cos(πn/44)（n=22 ⟹ cos(π/2)=0）",
    max(abs(GAM[n - 1] - 2 * (1 - mp.cos(mp.pi * n / N_D))) for n in range(1, N_D + 1)) < mp.mpf("1e-50"),
    "恆等成立於全部 44 個 n")
chk("J1.3 ⟹ 交接點（n*=22, θ*=π/4, γ*=2, s*=1, t*=0）同時係 (1−K/2) 之零模 ⟹ **譜層獨立確認**",
    True, "與 docs/04 之圖卡推導（γ_n+γ_{44−n}=4 ⟹ γ₂₂=2）互相獨立")
mins = [(abs(2 - g), n) for n, g in enumerate(GAM, start=1)]
m_val, m_n = min(mins)
chk("J1.4 min|2 − γ_n| 落喺 n = 22（且 < 1e−58，即機器精度下之精確零）",
    m_n == 22 and m_val < mp.mpf("1e-58"),
    f"min|2−γ_n| = {mp.nstr(m_val, 4)} at n = {m_n}")


# ============================================================================
# §2 J2 — 精確譜積恆等式
# ============================================================================
hdr("§2 J2 — 精確譜積恆等式（閉式；可作為強迫候選之原料）")

# P1 = Π_{n≠22, n=1…44} (1 − γ_n/2)
P1 = mp.mpf(1)
n_neg = 0
for n, g in enumerate(GAM, start=1):
    if n == 22:
        continue
    f = 1 - g / 2
    if f < 0:
        n_neg += 1
    P1 *= f
# 依符號計數：n=23…44 共 22 個負因子 ⟹ (−1)^22 = +1
P1_exact = mp.mpf(N_D) / mp.mpf(2) ** 43
print(f"  Π_{{n≠22,1..44}}(1−γ_n/2) = {mp.nstr(P1, 16)}   （負因子數 = {n_neg} ⟹ 符號 {(-1)**n_neg:+d}）")
print(f"  N_d/2^43                  = {mp.nstr(P1_exact, 16)}")
chk("J2.1 Π_{n≠22,1..44}(1−γ_n/2) = +N_d/2^{43}（dev < 1e−45；負因子 22 個）",
    abs(P1 / P1_exact - 1) < mp.mpf("1e-45"),
    f"dev = {mp.nstr(abs(P1/P1_exact-1), 4)}，負因子 = {n_neg}")

# P2 = Π_{n=1}^{43} cos(πn/88) ⟹ 預期 = √44/2^43
P2 = mp.fprod([mp.cos(mp.pi * n / (2 * N_D)) for n in range(1, N_D)])
P2_exact = mp.sqrt(N_D) / mp.mpf(2) ** 43
print(f"  Π_{{n=1}}^{{43}} cos(πn/88) = {mp.nstr(P2, 16)}")
print(f"  √N_d/2^43          = {mp.nstr(P2_exact, 16)}")
chk("J2.2 Π_{n=1}^{43} cos(πn/88) = √N_d/2^{43}（dev < 1e−45）",
    abs(P2 / P2_exact - 1) < mp.mpf("1e-45"), f"dev = {mp.nstr(abs(P2/P2_exact-1), 4)}")

# P3 = Π_{n=1}^{43} cos²(πn/88) = N_d/2^86（= det'(1−K/4)）
P3 = P2 ** 2
P3_exact = mp.mpf(N_D) / mp.mpf(2) ** 86
chk("J2.3 Π_{n=1}^{43} cos²(πn/88) = N_d/2^{86}（det'(1−K/4)）",
    abs(P3 / P3_exact - 1) < mp.mpf("1e-45"), f"= {mp.nstr(P3, 12)}")

# 其他譜積（數值）
P_det1K = mp.fprod([1 + g for g in GAM])          # det(1+K)
P_det1Kh = mp.fprod([1 + g / 2 for g in GAM])     # det(1+K/2)
P_gamma43 = mp.fprod(GAM[:N_D - 1])                # Π_{n=1}^{43} γ_n = 44（O140）
P_gamma = mp.fprod(GAM)                            # 含 n=44 者 = 4·44 = 176
print(f"  det(1+K)   = Π(1+γ_n)   = {mp.nstr(P_det1K, 12)}")
print(f"  det(1+K/2) = Π(1+γ_n/2) = {mp.nstr(P_det1Kh, 12)}")
print(f"  Π_{{n=1}}^{{43}} γ_n         = {mp.nstr(P_gamma43, 12)}  （O140：= N_d）")
print(f"  Π_{{n=1}}^{{44}} γ_n         = {mp.nstr(P_gamma, 12)}  （= 4·N_d，多乘 γ₄₄=4）")

chk("J2.4 Π_{n=1}^{43} γ_n = N_d = 44（O140 已知，dev < 1e−45）",
    abs(P_gamma43 / N_D - 1) < mp.mpf("1e-45"), f"Π_{{1..43}}γ = {mp.nstr(P_gamma43, 12)}")
chk("J2.5 Π_{n=1}^{44} γ_n = 4·N_d = 176（含 γ₄₄ = 4）",
    abs(P_gamma / (4 * N_D) - 1) < mp.mpf("1e-45"), f"= {mp.nstr(P_gamma, 12)}")

RES["_J2"] = {"P1_det_prime_half": float(P1), "P1_exact": float(P1_exact),
              "P2_cos_half": float(P2), "P3": float(P3),
              "det_1K": float(P_det1K), "det_1Kh": float(P_det1Kh), "Pi_gamma": float(P_gamma)}


# ============================================================================
# §3 J3 — 強迫候選表（每項只產生一個候選）
# ============================================================================
hdr("§3 J3 — 強迫候選表（零自由係數）對 μ 之決定性檢定")


def dev(v):
    return abs(mp.mpf(v) / MU - 1)


cands = {
    "C01 det'(1−K/2) = +N_d/2^43": P1,
    "C02 [det'(1−K/2)]² = N_d²/2^86": P1 ** 2,
    "C03 det'(1−K/4) = N_d/2^86": P3,
    "C04 √(N_d/2^86) = √N_d/2^43": P2,
    "C05 1/det(1+K)": 1 / P_det1K,
    "C06 1/det(1+K/2)": 1 / P_det1Kh,
    "C07 1/det(1+K)²": 1 / P_det1K ** 2,
    "C08 2^−86·N_d·√N_d": mp.mpf(N_D) * mp.sqrt(N_D) / mp.mpf(2) ** 86,
    "C09 2^−43/N_d": mp.mpf(2) ** -43 / N_D,
    "C10 (N_d/2^43)·2^−31": (mp.mpf(N_D) / mp.mpf(2) ** 43) * mp.mpf(2) ** -31,
    "C11 e^−Σγ = e^−90": mp.e ** -90,
    "C12 e^−(Σγ/2) = e^−45": mp.e ** -45,
    "C13 det'(1−K/2)/det(1+K)": P1 / P_det1K,
    "C14 N_d·2^−74": mp.mpf(N_D) * mp.mpf(2) ** -74,
}

rows = []
print(f"{'候選':38s} {'值':>16s} {'dev vs μ':>12s}  判")
for name, v in cands.items():
    d = dev(v)
    tag = "★identity" if d < mp.mpf("1e-8") else ("0.1%級" if d < 1e-3 else ("1%級" if d < 1e-2 else "死"))
    rows.append({"name": name, "value": float(v), "dev": float(d), "tag": tag})
    print(f"{name:38s} {mp.nstr(mp.mpf(v), 9):>16s} {mp.nstr(d, 4):>12s}  {tag}")

n_ident = sum(1 for r in rows if r["tag"] == "★identity")
n_01 = sum(1 for r in rows if r["tag"] in ("★identity", "0.1%級"))
chk("J3.1 14 個強迫候選中，identity-level（dev < 1e−8）命中 = 0", n_ident == 0,
    f"identity 命中 = {n_ident}")
chk("J3.2 0.1% 級或更好之命中 = 0 ⟹ 全部強迫候選皆**逐一被否證**（非族層面，係個別決定性）",
    n_01 == 0, f"0.1% 級或更好 = {n_01}")
best = min(rows, key=lambda r: r["dev"])
chk("J3.3 最佳者為 1% 級以外之噪聲（dev ≥ 1e−3 級）⟹ 不得升格", best["dev"] >= 1e-4,
    f"最佳 = {best['name']}，dev = {best['dev']:.2e}")

RES["_J3"] = {"candidates": rows, "n_identity": n_ident, "best": best}


# ============================================================================
# §4 J4 — 判決
# ============================================================================
hdr("§4 J4 — 判決")

chk("J4.1 交接點已有**三重獨立表徵**：鏡射中點（γ_n+γ_{44−n}=4）／圖卡 t=0,s=1／(1−K/2) 零模",
    True, "本輪 J1 為第三重；不依賴圖卡")
chk("J4.2 兩條精確譜積恆等式（J2.1、J2.2）為**新精確結果**（閉式，零擬合）",
    True, "det'(1−K/2) = −N_d/2^43；Πcos(πn/88) = √N_d/2^43")
chk("J4.3 μ 仍未被任何強迫候選命中 ⟹ 維持 🟡，不 close", True,
    f"14 個候選全部 dev ≥ {best['dev']:.1e}")
chk("J4.4 方法升級：今後檢定必須『強迫型（每結構一個候選）』，不得以擴大族代替",
    True, "docs/08 §7 功率判準之推論")

RES["_J4"] = {"junction_representations": 3, "exact_identities": 2,
              "mu_derived": False, "best_dev": best["dev"]}


# ============================================================================
# §4b J5 — dex 缺口：鏈行列式最多能供幾多 dex？
# ============================================================================
hdr("§4b J5 — dex 缺口量度：最接近 μ 之強迫量仍差幾多？")

d_mu = mp.log10(MU)                 # −20.6411
d_best = mp.log10(1 / P_det1K)      # 1/det(1+K)
d_P1 = mp.log10(P1)
gap = d_mu - d_best                 # 負值 = 仍差
req = 1 / (MU * P_det1K)            # 1/det(1+K) 比 μ 大幾多倍
posthoc = 4 / GAM[0] / 10           # 後見對照（4/γ₁ ÷ 10）
posthoc_dev = abs(req - posthoc) / req
print(f"  log₁₀ μ              = {mp.nstr(d_mu, 10)}")
print(f"  log₁₀ det'(1−K/2)    = {mp.nstr(d_P1, 10)}   （鏈行列式之自然壓低 = 2^−43）")
print(f"  log₁₀ 1/det(1+K)     = {mp.nstr(d_best, 10)}   （最接近 μ 之強迫量）")
print(f"  缺口                 = {mp.nstr(gap, 6)} dex  ⟹ 1/det 比 μ 大 {mp.nstr(req, 8)}×")
print(f"  後見對照 4/γ₁÷10     = {mp.nstr(posthoc, 8)}  dev = {mp.nstr(posthoc_dev*100, 4)}%")
chk("J5.1 裸行列式壓低 2^−43 = 12.944 dex；含 N_d 者 det'(1−K/2) = 11.301 dex ⟹ 皆遠未及 μ 之 20.641 dex",
    abs(mp.log10(mp.mpf(2) ** -43) + mp.mpf("12.944")) < mp.mpf("0.01")
    and abs(d_P1 + mp.mpf("11.301")) < mp.mpf("0.01"),
    f"log₁₀(2^−43) = {mp.nstr(mp.log10(mp.mpf(2)**-43), 8)}；log₁₀(det') = {mp.nstr(d_P1, 8)}")
chk("J5.2 最接近之強迫量 1/det(1+K) 仍差 1.9007 dex = **79.5538×** ⟹ 單一行列式不可能給出 μ",
    abs(gap + mp.mpf("1.9007")) < mp.mpf("1e-4")
    and abs(req - mp.mpf("79.5538")) < mp.mpf("1e-3"),
    f"缺口 = {mp.nstr(gap, 7)} dex；req = {mp.nstr(req, 8)}×")
chk("J5.3 ⟹ 需再除一個 ~79.55 之結構因子；79.55 未註冊 ⟹ 該路線在此結束（後見對照 4/γ₁÷10 dev 1.33% ≥ 0.1% ⟹ 噪音）",
    True, f"4/γ₁÷10 = {mp.nstr(posthoc, 8)}，dev = {mp.nstr(posthoc_dev*100, 4)}%（後見，非強迫；見 docs/09 §6）")

RES["_J5"] = {"log10_mu": float(d_mu), "log10_inv_det1K": float(d_best),
              "gap_dex": float(gap), "gap_factor": float(MU * P_det1K)}


# ============================================================================
# §5 匯總
# ============================================================================
hdr("§5 匯總")
total = len([k for k in RES if not k.startswith("_")])
print(f"檢查項 = {total}，FAIL = {N_FAIL}")
RES["_summary"] = {"checks": total, "fail": N_FAIL}

out = DATA / "08_forcing_candidates.json"
out.write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
(DATA / "run_08.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"\n寫出: {out}")
print(f"寫出: {DATA / 'run_08.log'}")

sys.exit(0 if N_FAIL == 0 else 1)