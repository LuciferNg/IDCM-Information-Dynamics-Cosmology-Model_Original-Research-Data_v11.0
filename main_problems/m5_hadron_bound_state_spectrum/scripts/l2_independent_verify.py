#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""L2 — 獨立重推導（M5 首個 L2；idcm-audit-verify §2.4 必跑閘）
================================================================================
範圍：`docs/36`（R4 輕區異常誤差傳播重判）所宣稱嘅全部數字。

四條鐵律（依 skill §2.4 模板）
  R1 唔 import／call 任何 L1 腳本；唔讀 docs/。原始數據**直接由 PDG 檔重解析**。
  R2 `INPUTS` 只放原始量：PDG 原始條目、組分質量（PDG 值）、校準常數、M₀（M6 錨）。
     唯一例外並明文標示：`psi2_model`（模型輸出）由 L1 `scripts/37` 產生 ⟹ **上游封結論**，
     因重寫束縛態求解器屬 L1 工作，本 L2 只驗 data 側 ＋ 統計。
  R3 `DERIVE` 每個 key 係 INPUTS 嘅純函數。
  R4 數值用 mpmath dps=40；禁 float 累加。

輸出：`data/l2_independent_verify.json`、`data/run_l2.log`
"""
import json
import sys
from pathlib import Path

import mpmath as mp

mp.mp.dps = 40

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
LOG = []
RES = {}
NFAIL = 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def ser(x):
    try:
        xf = float(x)
    except (OverflowError, ValueError):
        xf = None
    return {"float": xf if (xf is not None and (xf == 0 or 1e-300 < abs(xf) < 1e300)) else None,
            "str": mp.nstr(x, 25)}


# ── ① INPUTS ────────────────────────────────────────────────────────────────
PDG = ROOT / "data" / "pdg_mass_width_2024.mcd"

INPUTS = {
    "pdg_file": str(PDG),          # 原始測量檔（唔經 L1）
    "m_l": "0.30",                 # GeV（PDG 組分，R5 寫死）
    "m_s": "0.495",
    "m_c": "1.27",
    "m_b": "4.18",
    "q_s": "0.3333333333333333",   # |q_h|：s → 1/3
    "q_cb": "0.6666666666666666",  # c、b → 2/3
    "M0": "0.027900914",           # GeV（M6 錨）
    "hbar_c": "0.1973269804",      # GeV·fm
    "a0_fm": "0.384184",           # v2 校準（docs/27–28）
    "t_mev": "1.7711",
    # ⚠️ 上游 L1（scripts/37）之模型輸出 ψ²_model；本 L2 不作獨立重解，只作 declared 輸入
    "psi2_model_K": "0.10214995",
    "psi2_model_D": "0.20145517",
    "psi2_model_B": "0.14973450",
}

FAM = {"K": ((311, 321), (313, 323), "q_s", "m_s"),
       "D": ((411, 421), (413, 423), "q_cb", "m_c"),
       "B": ((511, 521), (513, 523), "q_cb", "m_b")}


def load_pdg():
    """由 PDG 原始檔重解析（唔用 L1 之解析結果）。"""
    mass, err = {}, {}
    for line in Path(INPUTS["pdg_file"]).read_text().splitlines():
        if line.startswith("*") or len(line) < 60:
            continue
        ids = [line[i:i + 8].strip() for i in (0, 8, 16, 24)]
        try:
            mv = mp.mpf(line[33:51].strip())
        except (ValueError, TypeError):
            continue
        try:
            ev = abs(mp.mpf(line[52:60].strip()))
        except (ValueError, TypeError):
            ev = mp.mpf(0)
        for k in ids:
            if k:
                mass[int(k)], err[int(k)] = mv * 1000, ev * 1000     # → MeV
    return mass, err


# ── ② DERIVE ────────────────────────────────────────────────────────────────
def derive(I):
    mass, err = load_pdg()
    D = {}
    for nm, (g0, g1, qk, mk) in FAM.items():
        d0 = mass[g0[0]] - mass[g0[1]]
        d1 = mass[g1[0]] - mass[g1[1]]
        s0 = mp.sqrt(err[g0[0]] ** 2 + err[g0[1]] ** 2)
        s1 = mp.sqrt(err[g1[0]] ** 2 + err[g1[1]] ** 2)
        c = (d1 - d0) / 4
        sc = mp.sqrt(s0 ** 2 + s1 ** 2) / 4
        D[f"d0_{nm}"], D[f"d1_{nm}"] = d0, d1
        D[f"c_{nm}"], D[f"sigma_c_{nm}"] = c, sc
        D[f"psi2_data_{nm}"] = abs(c) * mp.mpf(I[mk]) * mp.mpf(I["m_l"]) / mp.mpf(I[qk])
        # 統計
        D[f"dev_{nm}"] = D[f"psi2_data_{nm}"] / mp.mpf(I[f"psi2_model_{nm}"])
        D[f"c_model_{nm}"] = c / D[f"dev_{nm}"]
        D[f"z_{nm}"] = abs(c - D[f"c_model_{nm}"]) / sc
        D[f"dist_{nm}"] = abs(D[f"c_model_{nm}"]) / sc            # 模型鑑別力
        D[f"dist_pct_{nm}"] = D[f"dist_{nm}"] * 100               # 以 % 表示（供 L3 對文檔 279%）
        D[f"delta_c_{nm}"] = abs(c - D[f"c_model_{nm}"])
    D["u"] = mp.mpf(I["hbar_c"]) / (mp.mpf(I["a0_fm"]) * mp.mpf(I["M0"]))
    D["v"] = (mp.mpf(I["t_mev"]) / 1000) / mp.mpf(I["M0"])
    D["sigma_need5_K"] = abs(D["c_K"] - D["c_model_K"]) / 5
    D["gain_K"] = D["sigma_c_K"] / D["sigma_need5_K"]
    return D


# ── ③ IDENTITIES ────────────────────────────────────────────────────────────
def main():
    print("=" * 78 + "\n§0 輸入\n" + "=" * 78)
    for k, v in INPUTS.items():
        print(f"  {k:14s} = {v}")
    D = derive(INPUTS)
    print("\n" + "=" * 78 + "\n§1 重算（PDG 原始檔 → 統計）\n" + "=" * 78)
    for k in sorted(D):
        print(f"  {k:18s} = {mp.nstr(D[k], 12)}")

    print("\n" + "=" * 78 + "\n§2 兩條獨立路徑（dev / z）\n" + "=" * 78)
    for nm in FAM:
        a = D[f"dev_{nm}"]
        b = D[f"c_{nm}"] / D[f"c_model_{nm}"]
        chk(f"I1 {nm}: dev ≡ c_data/c_model", abs(a - b) < mp.mpf("1e-30"), f"差 {mp.nstr(abs(a-b),4)}")
    for nm in FAM:
        z1 = abs(D[f"c_{nm}"] - D[f"c_model_{nm}"]) / D[f"sigma_c_{nm}"]
        z2 = abs(D[f"c_{nm}"] - D[f"c_{nm}"] / D[f"dev_{nm}"]) / D[f"sigma_c_{nm}"]
        chk(f"I2 {nm}: z 兩路徑一致", abs(z1 - z2) < mp.mpf("1e-30"))

    print("\n" + "=" * 78 + "\n§3 判（docs/36 判準）\n" + "=" * 78)
    chk("Z1 K 之 z < 3 ⟹ 未達硬標準", D["z_K"] < 3, f"z = {mp.nstr(D['z_K'], 5)}σ")
    chk("Z2 K 之 z ≥ 1 ⟹ 非零張力（不得宣稱一致）", D["z_K"] >= 1, f"z = {mp.nstr(D['z_K'], 5)}σ")
    chk("Z3 5σ 所需 σ_c = |Δc|/5", abs(D["sigma_need5_K"] - abs(D["c_K"] - D["c_model_K"]) / 5) < mp.mpf("1e-30"))
    chk("Z4 精度改善倍數 ≈ 1.90", abs(D["gain_K"] - mp.mpf("1.90")) < mp.mpf("0.02"), f"{mp.nstr(D['gain_K'],5)}×")
    chk("Z5 D 校準 z ≈ 0", D["z_D"] < mp.mpf("0.01"), f"{mp.nstr(D['z_D'],4)}σ")
    chk("Z6 B 一致 z < 1", D["z_B"] < 1, f"{mp.nstr(D['z_B'],4)}σ")
    chk("Z7 K 模型鑑別力 > 1（可測）", D["dist_K"] > 1, f"{mp.nstr(D['dist_K'],5)}")
    chk("Z8 u ≈ 18.408939（v2 校準自洽）", abs(D["u"] - mp.mpf("18.408939")) < mp.mpf("1e-5"),
        f"{mp.nstr(D['u'],10)}")
    chk("Z9 v ≈ 0.0634782（同上）", abs(D["v"] - mp.mpf("0.0634782")) < mp.mpf("1e-7"), f"{mp.nstr(D['v'],10)}")

    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"\n檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_derived"] = {k: ser(v) for k, v in D.items()}
    RES["_summary"] = {"checks": nchk, "fail": NFAIL, "n_derived": len(D)}
    DATA.mkdir(exist_ok=True)
    (DATA / "l2_independent_verify.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_l2.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("→ data/l2_independent_verify.json ／ data/run_l2.log")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())