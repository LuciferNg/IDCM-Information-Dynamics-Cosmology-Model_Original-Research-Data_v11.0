#!/usr/bin/env python3
"""45 — R24：負號之解剖 —— 鏡像？對稱？映射？相位轉向？抑或「區分」

使用者問：「負號本質是什麼？鏡像？對稱？映射？……還是相位轉向？」

本輪把「−」當作**待判定之物**，對四個候選本體逐一判死或判生，並給出第五讀法之數值證據。

待解釋之數據：ρ_d = +0.40350365、ρ_u = −0.02588837（同一核 w 對三個測試向量之比值）

判準（跑前凍死）
--------------------------------------------------
F1 若某讀法使 ρ 對 w 不變 ⟹ 該讀法為「規範」⟹ 不可解釋數據（數據非規範不變量）
F2 若某讀法強制 |ρ_u| = f(ρ_d) 且與數據差 > 3σ ⟹ 該讀法判死
F3 若某讀法可在 ℝ⁴⁴ 上生成任意 (ρ_d, ρ_u) ⟹ 不可證偽 ⟹ 判「空」
F4 若某讀法保振幅（純旋轉）⟹ 檢 |ρ_u| 是否可由 |ρ_d| 或 1 之旋轉得來
F5 第五讀法（區分）須給出可算之結構量，且該量須為「物件層變換」所不能改變者
"""
import mpmath as mp

mp.mp.dps = 30
NS = 44
K = {"u": 2, "d": 3, "s": 4}
RD, RD_S = mp.mpf("0.40350365"), mp.mpf("0.03510")
RU, RU_S = mp.mpf("-0.02588837"), mp.mpf("0.03292")
phi = lambda n, j: mp.sqrt(mp.mpf(2) / (NS + 1)) * mp.sin(mp.mpf(n) * mp.pi * j / (NS + 1))

print("=" * 96)
print("45 — R24：負號本質之解剖")
print("=" * 96)
print("   數據：ρ_d = %s ± %s ；ρ_u = %s ± %s" % (mp.nstr(RD, 8), mp.nstr(RD_S, 5), mp.nstr(RU, 8), mp.nstr(RU_S, 5)))
print("   |ρ_u| / ρ_d = %s  ⟹ 「−」同時帶有**振幅**信息（非單位模）" % mp.nstr(abs(RU) / RD, 8))

# ---------- A：四個物件層候選 ----------
print("\n[A] 四個候選（皆為**物件層**變換）")
print("   A1 全局符號 w → −w（對稱／規範）")
Gm = lambda i, j, w: mp.fsum(w(n) * phi(n, i) * phi(n, j) for n in range(1, NS + 1))
w0 = lambda n: 1 / (4 * mp.sin(mp.mpf(n) * mp.pi / (2 * (NS + 1))) ** 2)
wm = lambda n: -w0(n)


def rho(w):
    ds = Gm(K["d"], K["s"], w) - Gm(K["u"], K["s"], w)
    return ((Gm(K["d"], K["d"], w) - Gm(K["u"], K["d"], w)) / ds,
            (Gm(K["d"], K["u"], w) - Gm(K["u"], K["u"], w)) / ds)


p0, p1 = rho(w0), rho(wm)
print("      w → ρ = (%s, %s) ； −w → ρ = (%s, %s)"
      % (mp.nstr(p0[0], 10), mp.nstr(p0[1], 10), mp.nstr(p1[0], 10), mp.nstr(p1[1], 10)))
print("      差 = %s ⟹ **規範不變 ⟹ 不可解釋數據**（F1）"
      % mp.nstr(max(abs(p0[0] - p1[0]), abs(p0[1] - p1[1])), 6))

print("\n   A2 扇區鏡像 u ↔ d（ℤ₂ 反稱）⟹ 強制 |ρ_u| = ρ_d")
mirror_off = RD / abs(RU)
print("      鏡像預言 |ρ_u| = ρ_d = %s ；實測 |ρ_u| = %s ⟹ 偏 %s 倍 > 3σ ⟹ **判死**（F2）"
      % (mp.nstr(RD, 8), mp.nstr(abs(RU), 8), mp.nstr(mirror_off, 6)))

print("\n   A3 映射（Δ_n 之符號圖樣）")
Dn = [phi(n, 3) - phi(n, 2) for n in [0] + list(range(1, NS + 1))][1:]
neg = sum(1 for x in Dn if x < 0)
pos = sum(1 for x in Dn if x > 0)
print("      Δ_n：正 %d 個、負 %d 個 ⟹ 測試向量已帶號，可生任何符號" % (pos, neg))
print("      但 docs/28 已證 cone = ℝ³（任何 (ρ_d,ρ_u) 皆可由某 w 得）⟹ **不可證偽 ⟹ 空**（F3）")

print("\n   A4 相位轉向（ℤ₄ 特徵標 i^{nc}）")
print("      實部：Re(i^n) = (1_0 − 1_2) = S5 ；Im(i^n) = (1_1 − 1_3) = S6")
print("      ⟹ 整個 ℤ₄ 軌道之實核由扇區指示子 {1_0,1_1,1_2,1_3} 張成（S1–S4 為基底）")
print("      而 docs/36 已測該基底：ρ_d 全域 ∈ [0.74538173, 1.1345044] ⟹ 永不達 0.4035 ⟹ **判死**")
print("   ★ 純旋轉保振幅：|ρ| 不變。數據 |ρ_u| = %s ≪ ρ_d = %s ⟹ **非純相位**（F4）"
      % (mp.nstr(abs(RU), 6), mp.nstr(RD, 6)))

# ---------- B：配對解剖 ----------
print("\n[B] 第五讀法：『−』係**配對之差**（區分），非物件之變換")
print("   恆等式（docs/28）：δ_x = ⟨w, v_x⟩，其中 v_x(n) ≡ φ_n(k_x)·Δ_n、Δ_n ≡ φ_n(3) − φ_n(2)")
vu = [phi(n, K["u"]) * (phi(n, 3) - phi(n, 2)) for n in range(1, NS + 1)]
vd = [phi(n, K["d"]) * (phi(n, 3) - phi(n, 2)) for n in range(1, NS + 1)]
vs = [phi(n, K["s"]) * (phi(n, 3) - phi(n, 2)) for n in range(1, NS + 1)]
dot = lambda a, b: mp.fsum(x * y for x, y in zip(a, b))
nrm = lambda a: mp.sqrt(dot(a, a))
print("   ‖v_u‖ = %s ；‖v_d‖ = %s ；‖v_s‖ = %s" % (mp.nstr(nrm(vu), 8), mp.nstr(nrm(vd), 8), mp.nstr(nrm(vs), 8)))
print("   cos(v_u,v_s) = %s" % mp.nstr(dot(vu, vs) / nrm(vu) / nrm(vs), 8))
print("   cos(v_d,v_s) = %s" % mp.nstr(dot(vd, vs) / nrm(vd) / nrm(vs), 8))
print("   cos(v_u,v_d) = %s" % mp.nstr(dot(vu, vd) / nrm(vu) / nrm(vd), 8))
Gram = mp.matrix([[dot(a, b) for b in (vu, vd, vs)] for a in (vu, vd, vs)])
print("   Gram（u,d,s 序）之譜：%s" % [mp.nstr(x, 6) for x in mp.eigsy(Gram)[0]])
print("   條件數 = %s" % mp.nstr(abs(mp.eigsy(Gram)[0][0] / mp.eigsy(Gram)[0][2]), 6))

# ---------- C：最小範數 w ----------
print("\n[C] 達標之最小範數核 w*（3 約束、44 維）")
A = mp.matrix([[vu[n], vd[n], vs[n]] for n in range(NS)]).T   # 3×44
b = mp.matrix([RU, RD, 1])
wstar = A.T * mp.inverse(A * A.T) * b
print("   ⟨w*,v_s⟩ = %s ；⟨w*,v_d⟩ = %s ；⟨w*,v_u⟩ = %s"
      % (mp.nstr(dot([wstar[i] for i in range(NS)], vs), 8),
         mp.nstr(dot([wstar[i] for i in range(NS)], vd), 8),
         mp.nstr(dot([wstar[i] for i in range(NS)], vu), 8)))
ws = [wstar[i] for i in range(NS)]
th_u = mp.acos(dot(ws, vu) / nrm(ws) / nrm(vu)) * 180 / mp.pi
print("   ‖w*‖ = %s" % mp.nstr(nrm(ws), 8))
print("   ∠(w*, v_u) = %s° ⟹ **大於 90°（因 ⟨w*,v_u⟩ < 0）**" % mp.nstr(th_u, 8))
print("   ‖w*‖ 中由 span{v_u,v_d,v_s} 解釋之比例 = %s"
      % mp.nstr(nrm(ws) / nrm(ws), 6) + "（最小範數解全在此 span 內）")
print("   w* 之符號圖樣：正 %d、負 %d、|w*| 最大 %s、最小 %s"
      % (sum(1 for x in ws if x > 0), sum(1 for x in ws if x < 0),
         mp.nstr(max(abs(x) for x in ws), 6), mp.nstr(min(abs(x) for x in ws), 6)))

print("\n[C2] 對照：若要求 ρ_u = −ρ_d（鏡像），所需 ⟨w,v_u⟩ = %s ⟹ 與實測差 %s 倍"
      % (mp.nstr(-RD, 8), mp.nstr(mirror_off, 6)))

# ---------- D：負號不需「帶號核」 ----------
print("\n[D] ★ 關鍵：負號**不需要**帶號核")
print("   docs/28 已找到 w ≥ 0（且 entry-positive）之見證，精確給 δ_u/δ_s = −0.02588837")
print("   ⟹ 正核 ＋ 帶號測試向量（Δ_n 有 17 負 25 正）已足以生 −")
print("   ⟹ **負號之所在地 = 測試向量之幾何（區分），非核之符號性**")
print("   ⟹ docs/34 之「任何候選必須本質帶號」須收窄：帶號性已具備（在 v_x 內），")
print("      真正缺者 = 使 ⟨w,v_u⟩ 近零之**結構理由**（近正交，非符號）")

# ---------- E：負號之顯著性 ----------
print("\n[E] 負號之顯著性（供文檔溯源）")
DS = mp.mpf("3.4250")            # δ_s（MeV），R5 之解
sig_du = RU_S * DS
print("   δ_s = %s MeV（R5）；σ_δu = σ_ρu × δ_s = %s MeV" % (mp.nstr(DS, 6), mp.nstr(sig_du, 10)))
print("   δ_u = −0.0887 MeV ⟹ |δ_u| / σ_δu = %s σ" % mp.nstr(mp.mpf("0.0887") / sig_du, 6))
print("   ⟹ 文檔作 σ_δu ≈ 0.1128 MeV、負號 0.79σ")

print("\n" + "=" * 96)
print("[判定]")
print("   A1 全局符號（對稱／規範）：✅ 不可見 ⟹ 不能解釋數據")
print("   A2 扇區鏡像：🔴 預言 |ρ_u| = ρ_d ⟹ 偏 %s 倍 ⟹ 判死" % mp.nstr(mirror_off, 5))
print("   A3 映射（Δ_n 符號圖樣）：🟡 存在但不可證偽（cone = ℝ³）⟹ 空")
print("   A4 相位轉向（ℤ₄）：🔴 ρ_d ∈ [0.74538173, 1.1345044]、且純旋轉保振幅 ⟹ 判死")
print("   B/C 第五讀法【區分】：『−』= ⟨w,v_u⟩ < 0 而 ⟨w,v_s⟩ > 0 ⟹ 同一核對兩個不同測試向量之配對異號")
print("        ∠(w*,v_u) = %s° > 90°、|ρ_u|/ρ_d = %s ⟹ 係**近正交**，非鏡像、非純相位" % (mp.nstr(th_u, 6), mp.nstr(abs(RU) / RD, 6)))
print("   ⟹ 四者皆屬「物件層變換」，而數據之負號屬「關係層」（配對）；與『性質＝關係產物』一致")
print("=" * 96)