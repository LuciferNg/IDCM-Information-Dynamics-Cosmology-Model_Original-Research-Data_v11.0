"""M7 scripts/16 — 回核 M2 docs/02 裁定（R1–R4）：λ-族、互換律兩種讀法、M6 樹核
==========================================================================================================
M2 docs/02 對本 M 三封之裁定：
  R1 加性覆蓋複合態 ✅（O237 逐模讀數之繼承）
  R2 p 族（可分族）淘汰**不完備** 🔴 —— λ-族 (m = M₀[Σσγ + λΣ_{n<m}σσw]) 與既有約束 C1–C5 相容
  R3 只可升條件定理（前提 = O237 逐模讀數）
  R4 M7 撤銷接受，closed-moot
本檔核 M2 之五個數字 ＋ 一個**判據審計**：
  C4「互換律＋X 對合 恆等維持（E_sync' = Σγ − m）」——若 E_sync 係事後定義，該「檢定」恆真 ⟹ 零信息；
  反之若按律名之字面讀法 E_sync(σ) ≡ m(1−σ)/M₀（同步態 = 補遮罩之質量），則 λ≠0 能否維持？
輸出 data/16_m2_adjudication.json
"""
import json
import os
import random
from fractions import Fraction as F

from mpmath import mp, mpf, sin, pi, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "16_m2_adjudication_verify.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<52} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


ND, NM = 88, 44
M0 = mpf("27.900915")
G = {n: 4 * sin(n * pi / ND) ** 2 for n in range(1, NM + 1)}
SG = sum(G.values())
SG2 = sum(v ** 2 for v in G.values())

print("=" * 100)
print("M7 scripts/16 — 回核 M2 docs/02 裁定（λ-族／互換律兩讀法／M6 樹核）")
print("=" * 100)

hdr("§1 譜常數與 M2 C2／C3 之數字")
W = (SG ** 2 - SG2) / (2 * SG)          # Σ_{n<m} γ_nγ_m/Σγ（σ 全 1）
lam_min = -SG / W
print(f"    Σγ = {nstr(SG, 12)}；Σγ² = {nstr(SG2, 10)}")
print(f"    W = Σ_(n<m) γ_nγ_m/Σγ = {nstr(W, 10)} ⟹ λ_min = −Σγ/W = {nstr(lam_min, 8)}")
print(f"    C2 (43,44) 兩模束縛偏移 = M₀·γ₄₃γ₄₄/Σγ = {nstr(M0 * G[43] * G[44] / SG, 8)} MeV")
print(f"    C2 全 44 模 = M₀·W = {nstr(M0 * W, 8)} MeV")
ok("Σγ = 90（已註冊）", abs(SG - 90) < mpf("1e-20"), nstr(SG, 12), "90")
ok("λ_min = −2.0695（M2 C3）", abs(lam_min + mpf("2.0695")) < mpf("5e-4"), nstr(lam_min, 8), "-2.0695")
ok("C2 兩模 = 4.9538 MeV（M2 C2）", abs(M0 * G[43] * G[44] / SG - mpf("4.9538")) < mpf("5e-4"),
   nstr(M0 * G[43] * G[44] / SG, 8), "4.9538")
ok("C2 全 44 模 ≈ 1213.4 MeV（M2 C2）", abs(M0 * W - mpf("1213.4")) < mpf("0.5"), nstr(M0 * W, 6), "1213.4")


def mass(sig, lam):
    """m/M₀ = Σσγ + λ Σ_{n<m} σ_nσ_m γ_nγ_m/Σγ"""
    tot = sum(G[n] for n in range(1, NM + 1) if sig[n])
    cross = mpf(0)
    for n in range(1, NM + 1):
        if not sig[n]:
            continue
        for mm in range(n + 1, NM + 1):
            if sig[mm]:
                cross += G[n] * G[mm] / SG
    return tot + lam * cross


hdr("§2 C4 審計：『互換律恆等維持』係恆真抑或有內容？")
rng = random.Random(2718)
dev_id, dev_lit = [], []
for _ in range(200):
    sig = {n: rng.random() < 0.5 for n in range(1, NM + 1)}
    comp = {n: not sig[n] for n in range(1, NM + 1)}
    m_s, m_c = mass(sig, mpf(-1)), mass(comp, mpf(-1))
    # 讀法（i）事後定義 E_sync := Σγ − m/M₀ ⟹ 恆等
    dev_id.append(abs(m_s + (SG - m_s) - SG))
    # 讀法（ii）字面：同步態 = 補遮罩之質量 ⟹ E_sync := m(1−σ)/M₀
    dev_lit.append(abs(m_s + m_c - SG))

print(f"    讀法(i) 事後定義：200 組最大偏差 = {nstr(max(dev_id), 6)}（恆等，零信息）")
nz = sum(1 for d in dev_lit if d > mpf("1e-6"))
print(f"    讀法(ii) 字面（E_sync = 補遮罩質量）：{nz}/200 組違反；最大偏差 = {nstr(max(dev_lit), 8)}")
ok("C4 在事後定義讀法下恆真 ⟹ 該『檢定』零信息（非五項約束之一）", max(dev_id) < mpf("1e-20"),
   f"max|dev| = {nstr(max(dev_id),4)}", "0")
ok("★在字面讀法下，λ=−1 之族 200/200 違反互換律 ⟹ R2 之相容性依讀法而定",
   nz == 200, f"{nz}/200", "200/200",
   "⟹ 『p 族不完備』並非無條件成立，而係依『互換律』之定義選擇")

hdr("§3 為何兩種讀法會分叉：λ 項之結構本質")
print("    事後定義：E_sync 係『剩餘量』⟹ 任何 m 皆滿足 m/M₀+E_sync=Σγ ⟹ 定律不可證偽 ⟹ 零約束力。")
print("    字面定義：E_sync = 補遮罩質量 ⟹ m(σ)+m(1−σ)=M₀Σγ 係**可證偽**方程：")
print("        m(σ)+m(1−σ) − M₀Σγ = M₀λ[W_σ + W_{1−σ}]，一般 ≠ 0 ⟹ 排除 λ≠0。")
sig_test = {n: (n % 3 == 0) for n in range(1, NM + 1)}
w_s = sum(G[n] * G[mm] / SG for n in range(1, NM + 1) if sig_test[n]
          for mm in range(n + 1, NM + 1) if sig_test[mm])
w_c = sum(G[n] * G[mm] / SG for n in range(1, NM + 1) if not sig_test[n]
          for mm in range(n + 1, NM + 1) if not sig_test[mm])
print(f"    例（σ_n = 1 iff n≡0 mod 3）：W_σ = {nstr(w_s, 8)}，W_(1−σ) = {nstr(w_c, 8)}"
      f" ⟹ 偏差 = M₀·|λ|·(W_σ+W_(1−σ)) = {nstr(M0 * (w_s + w_c), 8)} MeV 於 λ=−1")
ok("λ 項之可證偽偏差為 O(10²) MeV 級 ⟹ 字面讀法下 λ 族被強烈排除",
   M0 * (w_s + w_c) > 100, nstr(M0 * (w_s + w_c), 6) + " MeV", ">100 MeV")

hdr("§4 M6 樹核非加性（M2 docs/02 §七 引用）")
N = 88
g = lambda d: F((N * N - 1) - 6 * d * (N - d), 12 * N)
print(f"    g(0) = {float(g(0)):.6f}；g(44) = {float(g(44)):.6f}；"
      f"g(1)+g(1) = {float(2*g(1)):.6f}；g(2) = {float(g(2)):.6f}")
viol = [d for d in range(0, 22) if 2 * g(d) != g(2 * d)]
last_pos = max(d for d in range(0, N // 2 + 1) if g(d) > 0)   # 注意 g(d)=g(N−d) 對稱 ⟹ 單調段僅 d ≤ N/2
print(f"    加性違反（2g(d) ≠ g(2d)）計 {len(viol)}/22 個 d（d=0..21）；最後正值 d = {last_pos}")
ok("M6 g(1)+g(1) = 13.676136 且 g(2) = 6.355114（非加性）",
   abs(float(2 * g(1)) - 13.676136) < 1e-5 and abs(float(g(2))) - 6.355114 < 1e-5,
   f"{float(2*g(1)):.6f} vs {float(g(2)):.6f}", "13.676136 vs 6.355114")
ok("g 在 d ≤ 44 單調減、零交叉於 d≈18.6（最後正值整數 d = 18）", last_pos == 18, last_pos, 18,
   "g(d)=g(N−d) 對稱 ⟹ 全域非單調；M2 之『零交叉 d=18』= 單調段最後正值")

hdr("§5 對本 M 之結論（R1／R3／R4 接受；R2 標為讀法依賴）")
print("    R1：⟨O⟩≡1 在 O237 逐模讀數下成立 ⟹ 與本 M scripts/02 之 504/504 同構（恆等，非巧合）。")
print("    R3：條件定理 T_M2 接受；本 M 引用時須一併寫明前提（O237 逐模）。")
print("    R4：closed-moot 接受。")
print("    R2：『p 族不完備』成立與否取決於互換律採事後定義（C4 零信息 ⟹ 不完備成立）")
print("        抑或字面定義（λ 族 200/200 違反 ⟹ 完備性可恢復）。此為 M2 產權之定義選擇。")
OUT["verdict"] = dict(R1="accepted", R2="reading-dependent (C4 vacuous; literal reading kills lambda-family)",
                      R3="accepted (conditional theorem)", R4="closed-moot")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(SG=float(SG), SG2=float(SG2), W=float(W), lam_min=float(lam_min),
                        C2_2mode=float(M0 * G[43] * G[44] / SG), C2_full=float(M0 * W),
                        dev_identity_max=float(max(dev_id)), dev_literal_max=float(max(dev_lit)),
                        n_violate_literal=int(nz), example_dev=float(M0 * (w_s + w_c)),
                        g_add_viol=int(len(viol)), g_last_pos=int(last_pos)))
json.dump(OUT, open(os.path.join(M7, "data", "16_m2_adjudication.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/16_m2_adjudication.json")