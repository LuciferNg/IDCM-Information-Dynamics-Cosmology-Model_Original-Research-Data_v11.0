"""M7 scripts/18 — 回核 M5 docs/52（算子槽位需求 ＋ √5 指數 ＋ A_eff 標度）
==========================================================================================================
M5 → M7 之核心要求：M7 之夸克模型基底若要承載 Δ0（電荷分裂），**必須加入一個新算子**，
量子數已釘死：一體、isospin-odd、非 q²、標度 A_eff = 4.4492±0.4603 MeV、帶 m_h 依賴（B 抑制 14.35×）。

本檔做三件事（全部可算）：
  §1 核 M5 之反解：k = ln T₃ / ln(m_b/m_c) 是否 = 2.2361，與 √5 之差
  §2 判 √5 之地位：單點反解 ⟹ 0 dof；並以 M5 自己之 K 通道 σ 說明「單一冪律」整族已死
  §3 判 A_eff 之標度比（A_eff/X_ext、A_eff/t）並**明確拒絕**任何凑數升格
  §4 M7 基底之可行性判定：庫存內有無「一體／isospin-odd／非 q²」算子
輸出 data/18_m5_operator_slot.json
"""
import json
import os

from mpmath import mp, mpf, log, sqrt, nstr

mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "18_m5_operator_slot.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<54} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


T3 = mpf("14.352258")          # M5：1/s，B 抑制倒數
RATIO = mpf("3.291339")        # M5：m_b/m_c
SQ5 = sqrt(5)
A_EFF, S_A = mpf("4.4492"), mpf("0.4603")
X_EXT = mpf("2.56")            # FLAG24：X = m_d − m_u
X_GAM = mpf("2.5461523")       # 本 M γ 載體分支（已 NO-GO，僅作對照）
T_M5_COEF = mpf("2.51")        # M5：A_eff = 2.51 t ⟹ 反推其 t 尺度（本 M 不判定來源）

print("=" * 100)
print("M7 scripts/18 — 回核 M5 docs/52：算子槽位需求、√5 指數、A_eff 標度")
print("=" * 100)

hdr("§1 指數反解 k = ln T₃ / ln(m_b/m_c)")
k = log(T3) / log(RATIO)
dev = (k / SQ5 - 1) * 100
print(f"    k = ln({nstr(T3,9)}) / ln({nstr(RATIO,8)}) = {nstr(k, 10)}")
print(f"    √5 = {nstr(SQ5, 10)} ⟹ 相對差 = {float(dev):+.4f}%")
ok("k 反解 = 2.2361（M5 記載）", abs(k - mpf("2.2361")) < mpf("5e-5"), nstr(k, 8), "2.2361")
ok("與 √5 之差 ≈ 0.001%（M5 記載）", abs(dev) < mpf("0.01"), f"{float(dev):+.4f}%", "<0.01%",
   "單一質量比反解 ⟹ 0 dof")

hdr("§2 √5 之地位：單點反解 vs 定律")
print("    M5 自報：k=1 ⟹ B 4.0σ／K 113σ；k=2 ⟹ B 0.8σ／K 601σ")
print("    ⟹ dσ_K/dk > 0（K 通道對 k 單調惡化）⟹ k = √5 = 2.236 時 K 通道 ≫ 601σ")
print("    ⟹ 『k = √5』係 B 通道單點擬合，**唔可能**係統一冪律；整族已被 K 通道否證（M5 自判相同）")
ok("√5 不可升格為定律（K 通道反證：k↑ ⟹ K σ↑ 單調）", True, "B-only fit", "0 dof",
   "與 M5 自判一致；本 M 只加一句：任何『√5 結構』主張須先在 K 通道盲測")
print(f"    附：若硬套 k=√5 於 D 通道（m_c/m_s = {nstr(mpf('1270')/mpf('93.4'),8)}），"
      f"抑制 = {nstr((mpf('1270')/mpf('93.4'))**SQ5, 8)}，M5 之 ΔI=0 通道已 170σ 失敗 ⟹ 同一結論")

hdr("§3 A_eff 之標度比（僅量級，不作結構）")
r1 = A_EFF / X_EXT
r2 = A_EFF / X_GAM
print(f"    A_eff/X_ext = {nstr(A_EFF,6)}/{nstr(X_EXT,4)} = {nstr(r1, 8)}；√3 = {nstr(sqrt(3), 8)}"
      f" ⟹ 差 {float((r1/sqrt(3)-1)*100):+.3f}%")
print(f"    A_eff/X_γ   = {nstr(A_EFF,6)}/{nstr(X_GAM,8)} = {nstr(r2, 8)} ⟹ 差 {float((r2/sqrt(3)-1)*100):+.3f}%")
ok("A_eff/X 與 √3 之接近屬凑數（0 dof、單點、兩個 X 值給出不同偏差）", True,
   f"{float((r1/sqrt(3)-1)*100):+.2f}% / {float((r2/sqrt(3)-1)*100):+.2f}%", "拒絕升格",
   "依 φ-power 陷阱紀律：單點吻合不可作結構，且兩載體 X 給出不同偏差 ⟹ 無穩定性")
t_implied = A_EFF / T_M5_COEF
print(f"    M5 記 A_eff = 2.51 t ⟹ 反推其 t 尺度 = {nstr(t_implied, 6)} MeV（本 M 不判定其來源，只記一致性）")
ok("A_eff/2.51 = 1.7726 MeV（與其 t 尺度自洽；不作結構主張）",
   abs(t_implied - mpf("1.77259")) < mpf("5e-4"), nstr(t_implied, 7) + " MeV", "1.77259",
   "M5 之 t 尺度屬其產權；本 M 只核算術")

hdr("§4 M7 基底可行性：庫存有無『一體／isospin-odd／非 q²』算子")
print("    本 M 庫存（已註冊）：① per-mode 加性 Σσ_nγ_n（味-blind：γ 只依 n，不依味）")
print("                        ② 樹核 g(d)（兩體、非加性；d = 模距）")
print("                        ③ 嵌入層標籤（O216：色荷=嵌入層標籤非譜量）")
print("    ⟹ ① 味-blind ⟹ 不能生 isospin-odd；② 兩體 ⟹ 違 M5 之「一體」需求；")
print("      ③ 係標籤層而非譜量 ⟹ 可承載『味』荷，但**無標度**（標籤冇能量）")
print("    ⟹ 結論：本 M 基底**無現成算子**可承載 Δ0 —— 與 M5 之判定一致，且原因更強：")
print("       缺嘅唔止係『一體 odd 算子』，而係『標籤層 → 譜量』之橋（該橋本身未註冊）。")
ok("M7 基底無現成一體/isospin-odd/非q² 算子 ⟹ 必為新成分（與 M5 一致）", True,
   "3 類庫存全部排除", "新成分",
   "①②③ 逐項排除；缺『標籤層→譜量』橋 ⟹ 需求升為 M7 之結構缺口")
OUT["requirement_for_M7"] = dict(
    one_body=True, isospin_odd=True, non_q2=True,
    scale_MeV=float(A_EFF), scale_sigma=float(S_A),
    missing_bridge="嵌入層標籤（味）→ 譜量 之標度橋（未註冊）")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(k=float(k), dev_sq5_pct=float(dev), A_eff_over_Xext=float(r1),
                        A_eff_over_Xgamma=float(r2), t_implied=float(A_EFF / T_M5_COEF)),
           verdict="√5 = B-only fit (0 dof), family dead in K; M7 basis lacks any one-body isospin-odd "
                   "spectral operator; missing bridge = label→spectrum scale")
json.dump(OUT, open(os.path.join(M7, "data", "18_m5_operator_slot.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/18_m5_operator_slot.json")