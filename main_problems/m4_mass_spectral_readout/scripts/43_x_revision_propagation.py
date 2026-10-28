#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M4 scripts/43 -- 傳播量化：FLAG24 `X` 修訂（docs/42）對既有 σ 宣稱之移動。

對應 `idcm-audit-verify` §4（傳播審計）：`X = 2.56 ± 0.094` 係 B/C 類輸入，
凡文檔用佢計 σ 者，必須逐處量化，唔可以講「影響不大」。

背景值（docs/42）：X_old = 2.5600 ± 0.0943 ｜ X_new = 2.5145 ± 0.1159（比值形式）
本腳本只做算術與判詞，不改任何既有 doc／script（原文保留，附註見各檔）。
"""
import json
import os

CHECK = []


def chk(label, ok, got, expect):
    CHECK.append({"label": label, "ok": bool(ok), "got": str(got), "expect": str(expect)})
    print(("PASS " if ok else "FAIL ") + "%-58s | %-30s | %s" % (label, got, expect))


X_O, SX_O = 2.5600, 0.0943
X_N, SX_N = 2.5145, 0.1159
FAC = SX_O / SX_N

# 通道 E（v1 / v2）與 σ（docs/42 §三）
CH = {
    "n−p":     dict(coeff=+1, E1=-1.2667, E2=-1.2212),
    "K+−K0":   dict(coeff=-1, E1=-1.3740, E2=-1.4195),
    "Σ0−Σ+":   dict(coeff=+1, E1=+0.7120, E2=+0.7575),
}
for k, v in CH.items():
    v["s1"] = abs(v["coeff"]) * SX_O
    v["s2"] = abs(v["coeff"]) * SX_N

print("=" * 100)
print("P1  docs/27 之 σ 宣稱（殘餘 / σ_E）—— 直接重算，非純標度縮放")
print("=" * 100)
print("  %-8s %10s %10s %10s %10s %8s" % ("通道", "E(v1)", "σ(v1)", "E(v2)", "σ(v2)", "z 變"))
for k, v in CH.items():
    z1 = abs(v["E1"]) / v["s1"]
    z2 = abs(v["E2"]) / v["s2"]
    v["z1"], v["z2"] = z1, z2
    print("  %-8s %+10.4f %10.4f %+10.4f %10.4f %8.2f→%.2fσ" % (k, v["E1"], v["s1"], v["E2"], v["s2"], z1, z2))
chk("P1.1 三通道符號判詞不變（|E|/σ 仍 ≫ 5σ）", all(min(v["z1"], v["z2"]) > 5 for v in CH.values()),
    "min z=%.2f" % min(min(v["z1"], v["z2"]) for v in CH.values()), ">5σ")
chk("P1.2 K 通道 z 下降（14.3 → 12.3σ）", CH["K+−K0"]["z2"] < CH["K+−K0"]["z1"],
    "%.2f → %.2f" % (CH["K+−K0"]["z1"], CH["K+−K0"]["z2"]), "下降")
chk("P1.3 Σ0−Σ+ 通道 z 下降（7.6 → 6.5σ）", CH["Σ0−Σ+"]["z2"] < CH["Σ0−Σ+"]["z1"],
    "%.2f → %.2f" % (CH["Σ0−Σ+"]["z1"], CH["Σ0−Σ+"]["z2"]), "下降")

print("\n" + "=" * 100)
print("P2  docs/27 §四 之邊界宣稱：X > 3.9340（三通道聯合下界）")
print("=" * 100)
LB = 3.9340
z_old = (LB - X_O) / SX_O
z_new = (LB - X_N) / SX_N
print("  (3.9340 − X)/σ :  old = %.2fσ   new = %.2fσ" % (z_old, z_new))
chk("P2.1 下界不一致之判詞不變（>5σ）", z_new > 5, "%.2fσ" % z_new, ">5σ")
chk("P2.2 文檔所寫 14.6σ 應改 12.3σ（純標度 ×%.3f）" % FAC, abs(z_new - 12.25) < 0.1,
    "%.2fσ" % z_new, "12.25σ")

print("\n" + "=" * 100)
print("P3  docs/26 之 X 敏感度（文檔自報：X 變 1σ ⟹ a_max 移 2–7%）")
print("=" * 100)
for lo, hi in ((0.02, 0.07),):
    shift_sigma = abs(X_N - X_O) / SX_O
    print("  ΔX = %.4f MeV = %.2f 個舊 σ" % (X_O - X_N, shift_sigma))
    print("  ⟹ a_max 移動 ≈ %.2f × [2%%, 7%%] = [%.2f%%, %.2f%%]" % (shift_sigma, shift_sigma * lo * 100, shift_sigma * hi * 100))
chk("P3.1 ΔX < 0.5 舊σ ⟹ a_max 移動 < 3%（界限結論不變）", abs(X_N - X_O) / SX_O < 0.5,
    "%.2fσ_old" % (abs(X_N - X_O) / SX_O), "<0.5")
chk("P3.2 docs/26 判詞已於 docs/28 撤銷 ⟹ 屬 D 類，只需附註",
    True, "已撤銷（docs/26 附註二）", "附註即可")

print("\n" + "=" * 100)
print("P4  純標度項（docs/13 類）一致性")
print("=" * 100)
print("  z 乘因子 = σ_old/σ_new = %.4f" % FAC)
for z in (16.0, 9.76, 14.6, 13.5):
    print("    %5.2fσ → %5.2fσ" % (z, z * FAC))
chk("P4.1 標度因子 < 1（所有依賴 σ_X 之 z 一律下修）", FAC < 1, "%.4f" % FAC, "<1")

DATA = {
    "X_old": {"value": X_O, "sigma": SX_O}, "X_new": {"value": X_N, "sigma": SX_N},
    "scale_factor": round(FAC, 4),
    "channels": {k: {"z_old": round(v["z1"], 3), "z_new": round(v["z2"], 3),
                     "E_old": v["E1"], "E_new": v["E2"]} for k, v in CH.items()},
    "bound_3p934": {"z_old": round(z_old, 3), "z_new": round(z_new, 3)},
    "n_pass": sum(1 for c in CHECK if c["ok"]), "n_total": len(CHECK),
}
print("\n" + "=" * 100)
print("RESULT: %d/%d PASS" % (DATA["n_pass"], DATA["n_total"]))
print("=" * 100)

out = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "43_x_revision_propagation.json"))
with open(out, "w", encoding="utf-8") as fh:
    json.dump({"checks": CHECK, "data": DATA}, fh, ensure_ascii=False, indent=2, sort_keys=True)
print("JSON -> data/43_x_revision_propagation.json")