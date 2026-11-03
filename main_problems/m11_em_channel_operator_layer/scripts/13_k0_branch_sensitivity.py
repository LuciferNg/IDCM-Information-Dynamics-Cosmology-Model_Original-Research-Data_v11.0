#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M11 scripts/13 — k₀ 兩支敏感度（M5 docs/55 裁定後之並報實算）
================================================================
問題：M5 docs/55 判 `w=|v|²`（v2）為正典，並指 M11 之「v2 支污染」判詞範圍過寬。
      本腳本量化：**M11 全部判詞在 k₀ 兩支下會否翻轉**（唔靠斷言，逐格實算）。

輸入（只讀）：
  data/01_r1_criteria.json  之 `C9.channels[ch].box_min/box_max`（箱界，Fraction 級精確）
  data/04_r4_necessity.json 之 `models[].verdicts[ch].em`（模型 em，⟨O⟩ 已由錨釘死 ⟹ 與 k₀ 無關）
判準（事前寫死）：
  P1 逐格判 `lo ≤ em/k₀ ≤ hi`；兩支判定不同 ⟹ 翻轉。
  P2 判詞層面：模型「盲測可行」⟺ 其**全部**有預言之盲通道皆在箱內。
  P3 `k₀^crit`、`a^crit`、C9 符號封鎖與 k₀ 無關 ⟹ 預期零移動（作為對照，唔作判詞）。
輸出：data/13_k0_branch_sensitivity.json、log 由 driver/shell redirect 產生。
"""
import json
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

K0 = {"v1_docs13": F(38152622441, 10000000000),          # 3.8152622441 MeV（M11 docs/06 正典）
      "v2_m5docs28": F(3748112, 1000000)}                # 3.748112 MeV（M5 docs/28，docs/55 定為正典）
BLIND = ["S0-S+", "S--S+", "D+-D0"]


def main():
    d1 = json.load(open(os.path.join(DATA, "01_r1_criteria.json"), encoding="utf-8"))
    d4 = json.load(open(os.path.join(DATA, "04_r4_necessity.json"), encoding="utf-8"))
    ch = d1["C9"]["channels"]

    print("=" * 96)
    print("M11 scripts/13 — k₀ 兩支敏感度（v1 = M11 正典；v2 = M5 docs/55 新正典）")
    print("=" * 96)
    print(f"k₀(v1) = {float(K0['v1_docs13']):.10f} MeV ｜ k₀(v2) = {float(K0['v2_m5docs28']):.10f} MeV"
          f" ｜ 位移 = {(float(K0['v2_m5docs28'])/float(K0['v1_docs13'])-1)*100:+.3f}%")

    rows, flips, permodel = [], 0, {}
    for m in d4["models"]:
        lab = m["label"].split()[0]
        ok = {k: [] for k in K0}
        for c in BLIND:
            v = m["verdicts"].get(c)
            if not v or v.get("em") is None:
                continue
            em = v["em"]
            lo, hi = F(str(ch[c]["box_min"])), F(str(ch[c]["box_max"]))
            t1, t2 = em / float(K0["v1_docs13"]), em / float(K0["v2_m5docs28"])
            f1, f2 = float(lo) <= t1 <= float(hi), float(lo) <= t2 <= float(hi)
            flips += (f1 != f2)
            ok["v1_docs13"].append(f1)
            ok["v2_m5docs28"].append(f2)
            rows.append(dict(model=lab, ch=c, em=em, lo=float(lo), hi=float(hi),
                             t_v1=t1, t_v2=t2, feasible_v1=f1, feasible_v2=f2,
                             flip=bool(f1 != f2)))
            print(f"  {lab} {c:7s} em={em:+9.4f} 箱[{float(lo):+.3f},{float(hi):+.3f}]"
                  f"  t(v1)={t1:+8.4f} {'✅' if f1 else '✗'}   t(v2)={t2:+8.4f} {'✅' if f2 else '✗'}"
                  + ("   ★翻轉" if f1 != f2 else ""))
        permodel[lab] = {k: (all(v) if v else None) for k, v in ok.items()}

    print("\n[B] 模型層判詞（盲測可行 ⟺ 全部有預言之盲通道在箱內）")
    for lab, d in permodel.items():
        print(f"  {lab}: v1={d['v1_docs13']} → v2={d['v2_m5docs28']}"
              + ("   ★翻轉" if d["v1_docs13"] != d["v2_m5docs28"] else "   不變"))

    print("\n[C] 與 k₀ 無關之量（對照，不作判詞）")
    k0c = json.load(open(os.path.join(DATA, "05_r5_k0_identity.json"), encoding="utf-8"))
    print(f"  k₀^crit: " + " ".join(f"{k}={v:.6f}" for k, v in k0c["k0_crit"].items())
          + f" ｜ a^crit(all)={k0c['a_crit_all']:.8f} fm ⟹ 定義為 |em|/|箱界|，零 k₀ 依賴")
    print(f"  C9 符號封鎖：{d1['C9'].get('sign_blocked')} ⟹ 與 k₀、κ 無關")

    nch = len(rows)
    print(f"\n[D] 匯總：盲測格 {nch} 個｜翻轉 {flips} 個｜模型層翻轉 "
          f"{sum(1 for l, d in permodel.items() if d['v1_docs13'] != d['v2_m5docs28'])} 個")
    print("⟹ M11 之 R4 判詞在 k₀ 兩支下" + ("**不變**" if flips == 0 else "**有翻轉，須改判**"))

    out = dict(k0_v1=float(K0["v1_docs13"]), k0_v2=float(K0["v2_m5docs28"]),
               shift_pct=(float(K0["v2_m5docs28"]) / float(K0["v1_docs13"]) - 1) * 100,
               rows=rows, model_level=permodel, n_cells=nch, n_flips=flips,
               k0_crit=k0c["k0_crit"], a_crit_all=k0c["a_crit_all"],
               sign_blocked=d1["C9"].get("sign_blocked"),
               verdict="判詞不變（0 翻轉）" if flips == 0 else "判詞翻轉")
    with open(os.path.join(DATA, "13_k0_branch_sensitivity.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())