#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M14 scripts/04 — Q2②：`3! = 6` 個對應可否**只用結構**砍到唯一？

問題（record §二 Q2②）：
判準：A1 出處／A2 非循環（禁以觀測 δ 反解）／A3 零參數／A5 範圍／A6 味盲。

四項檢定：
  §一 軌道   ：6 個對應在 𝒜 下是否同一 gauge 軌道
  §二 輸出之標籤結構：結構所給者係「無標籤之多集」抑或「帶標籤之三元組」
  §三 單射   ：6 個對應之 (ρ_u,ρ_d) 是否互異（全報 N_trials = 6）
  §四 對號步驟：6 個之中「命中」者，係由結構抑或由觀測選出（A2 之形式化）

α = (1,2,4)：docs/10 之 x=2 唯一雙通成員（**引用**，不重掃；禁區 2）。
"""
from fractions import Fraction as F
from itertools import permutations
import mpmath as mp

mp.mp.dps = 40
print("=" * 78)
print("M14 scripts/04 — Q2②：3! 對應之可唯一化性")
print("=" * 78)

RD = F("0.40350365")
RU = F("-0.02588837")
SD = mp.mpf("0.044785178")     # M12 docs/42 R27 正典
SU = mp.mpf("0.033583899")


def rho(a):
    """ρ_u = (2α_u−α_d)/(2α_s)、ρ_d = (2α_d−α_u)/(2α_s)（A₂ Cartan ⊕ 1 之閉式）。"""
    au, ad, as_ = a
    return (2 * au - ad) / (2 * as_), (2 * ad - au) / (2 * as_)


def z_of(a):
    ru, rd = rho(a)
    ru_m, rd_m = mp.mpf(ru.numerator) / ru.denominator, mp.mpf(rd.numerator) / rd.denominator
    zu = (ru_m - mp.mpf(RU.numerator) / RU.denominator) / SU
    zd = (rd_m - mp.mpf(RD.numerator) / RD.denominator) / SD
    return ru, rd, zu, zd, zu * zu + zd * zd


# =====================================================================
print("\n【§一】6 個對應之地位：同屬一個 gauge 軌道")
print("  依 scripts/02 §四（L3）：𝒜 對味標籤置換精確不變（10⁷ null，bit-identical）")
print("  依 M13 docs/05 §S2（本輪 scripts/02 §二 已獨立重現）：ρ_x 只依 α，與三元組無關")
print("  ⟹ 『邊個味住邊個站』屬 𝒜 不變群之軌道 ⟹ 結構支應 0 bit")
print("  ⟹ |S₃| = 6 需 log₂6 = 2.585 bit；結構支應 0 ⟹ 缺口須由外部輸入補齊")

# =====================================================================
print("\n【§二】結構輸出之標籤結構：帶標籤？抑或無標籤多集？")
a = (F(1), F(2), F(4))
Mscal = mp.mpf(1) / 90
delta = [Mscal * (2 * a[0] - a[1]), Mscal * (2 * a[1] - a[0]), Mscal * 2 * a[2]]
print(f"  α=(1,2,4)、M=(1/90)(A₂⊕1) ⟹ δ = Mα = "
      f"({mp.nstr(delta[0],4)}, {mp.nstr(delta[1],4)}, {mp.nstr(delta[2],4)})")
print(f"  其比值多集 = {{{mp.nstr(delta[0]/delta[2],4)}, {mp.nstr(delta[1]/delta[2],4)}, 1.0}}")
print("  ⚠ 但 δ 之三個分量係按 **槽位** 索引（第 1/2/3 槽），而『槽位 ↔ 味名』本身即待定者")
print("  ⟹ 結構所能說者 = 『三個讀出值之比為 0 : 3/8 : 1』—— 一個**無味名之多集**")
print("  ⟹ 把 0 認作 δ_u、3/8 認作 δ_d，係**對號步驟**，不在結構之內")

# =====================================================================
print("\n【§三】6 個對應之預言（N_trials = 6，全報，不篩選）")
print(f"  觀測：ρ_d = {float(RD):.8f} ± {mp.nstr(SD,6)} ； "
      f"ρ_u = {float(RU):.8f} ± {mp.nstr(SU,5)}")
print(f"  {'#':>2s} {'(α_u,α_d,α_s)':>16s} {'ρ_u':>11s} {'ρ_d':>11s} {'z_u':>9s} {'z_d':>9s} {'χ²':>10s}  判")
base = [F(1), F(2), F(4)]
rows = []
for i, perm in enumerate(permutations(base), 1):
    ru, rd, zu, zd, chi2 = z_of(perm)
    rows.append((i, perm, ru, rd, zu, zd, chi2))
    ok = "✅ |z|<1" if max(abs(zu), abs(zd)) < 1 else ("🟡 2σ 內" if chi2 < 4 else "🔴 排除")
    print(f"  {i:>2d} {str(tuple(int(x) for x in perm)):>16s} {float(ru):>11.5f} {float(rd):>11.5f} "
          f"{mp.nstr(zu,5):>9s} {mp.nstr(zd,5):>9s} {mp.nstr(chi2,5):>10s}  {ok}")
pairs = {(r[2], r[3]) for r in rows}
print(f"  ⟹ 相異預言對 = {len(pairs)} / 6 ⟹ 映射**單射**（6 個樣本不可合併）")
surv = [r for r in rows if max(abs(r[4]), abs(r[5])) < 1]
print(f"  ⟹ 1σ 存活 = {len(surv)} / 6；存活者之 α 次序 = {tuple(int(x) for x in surv[0][1])}")

# =====================================================================
print("\n【§四】對號步驟之定性：存活者係由結構抑或由觀測選出？")
print("  結構之全部內容（§二）＝『比值多集 {0, 3/8, 1}』；此多集**不含味名**")
print(f"  觀測之比值 = (ρ_u, ρ_d) = ({float(RU):.6f}, {float(RD):.6f})，帶味名")
print("  對號所需 = 決定『0 ↔ ρ_u』與『3/8 ↔ ρ_d』⟹ 用『誰最接近 0』")
print("  候選：0 與 |ρ_u| 之距 = "
      f"{mp.nstr(abs(mp.mpf(RU.numerator)/RU.denominator),4)}；3/8 與 ρ_d 之距 = "
      f"{mp.nstr(abs(mp.mpf(3)/8 - mp.mpf(RD.numerator)/RD.denominator),4)}")
print("  ⟹ 對號判據 = 觀測值本身 ⟹ **A2（以觀測反解）**，非結構推導")
print("  ⟹ 結構對『哪一個對應』之貢獻 = 0 bit；6 個對應之篩選全由數據完成")
print("     亦即：本節所報之 1/6『命中』係**後驗對號**，不得作結構結果（記錄容許標記後驗，")
print("     但 Q2② 之問法係『只用結構砍』⟹ 答案為否）")

# =====================================================================
print("\n【§五】獨立重跑 M13 docs/07 之 216 掃描（6 規則 × 36 指派；全報）")
KS = {1: 3, 4: 2, 5: 0}          # v9 z2_map（mode → d）；docs/03 §一之 d 值
km = [1, 4, 5]
B = 4
RULES = {
    "R1 α=k":       lambda k, d: F(k),
    "R2 α=1/k":     lambda k, d: F(1, k),
    "R3 α=2^d":     lambda k, d: F(2) ** d,
    "R4 α=2^(−d)":  lambda k, d: F(1, 2 ** d),
    "R5 α=B^d":     lambda k, d: F(B) ** d,
    "R6 α=B^(−d)":  lambda k, d: F(1, B ** d),
}
win = (mp.mpf("3.3048"), mp.mpf("4.1300"))     # docs/05 §S3 之 1σ 窗
n_tot = 0
n_win = 0
n_joint = 0
win_vals = set()
for rname, g in RULES.items():
    for pk in permutations(range(3)):          # 哪個 k 給哪個味
        for pd in permutations(range(3)):      # 哪個 d 給哪個味
            a = tuple(g(km[pk[i]], KS[km[pd[i]]]) for i in range(3))
            if a[2] == 0:
                continue
            n_tot += 1
            ru, rd = rho(a)
            rus = mp.mpf(ru.numerator) / ru.denominator
            rds = mp.mpf(rd.numerator) / rd.denominator
            r = mp.mpf((a[2] / a[0]).numerator) / (a[2] / a[0]).denominator
            inwin = win[0] <= r <= win[1]
            zu = (rus - mp.mpf(RU.numerator) / RU.denominator) / SU
            zd = (rds - mp.mpf(RD.numerator) / RD.denominator) / SD
            if inwin:
                n_win += 1
                win_vals.add(str(a[2] / a[0]))
            if abs(zu) < 1 and abs(zd) < 1:
                n_joint += 1
print(f"  有效組合 = {n_tot}（216 減去 α_s = 0 之退化者）")
print(f"  α_s/α_u 落 1σ 窗 {win} 者 = {n_win} / {n_tot}"
      f"（{100*n_win/max(n_tot,1):.1f}%）；落窗之值 = {sorted(win_vals)}")
print(f"  **聯合閘** |z_u|<1 且 |z_d|<1 者 = {n_joint} / {n_tot}")
print(f"  ⟹ 與 M13 docs/07 §一 之所報（36/216 落窗，全部恰為 4.0；聯合閘 0/216）"
      f"{'逐項吻合' if (n_win == 36 and n_joint == 0 and win_vals == {'4'}) else '⚠ 不吻合，須裁差異'}")

# =====================================================================
print("\n" + "=" * 78)
print("Q2② 判決")
print("=" * 78)
print(f"""
  §一 6 個對應同屬一 gauge 軌道 ⟹ 結構支應 0 bit（需 2.585 bit，全數外來）
  §二 結構之輸出係**無味名之比值多集** {{0, 3/8, 1}} ⟹ 對號係額外一步
  §三 6 個預言互異（單射，{len(pairs)}/6）；1σ 存活 {len(surv)}/6
  §四 存活者之選出係靠『誰最接近 0』⟹ 用觀測 ⟹ A2

  ⟹ **Q2② = 死（不可只用結構砍至唯一）**
    範圍標籤（A5）：
      ① {{有限集：6 條已註冊規則 × 36 指派 = 216 —— M13 docs/07 報 0/216}}
      ② {{族：任何單函數 g:(k,d)→α —— M13 docs/09 §3.2 已證必帶 dof ≥ 1 ⟹ 違 A3}}
      ③ {{族：𝒜 之全部算子 —— 味置換精確對稱 ⟹ 0 bit（scripts/02 L3）}}
    四條路（規則枚舉／一般函數／算子層／數據對號）依次封閉，末者違 A2。
  ⟹ 唯一誠實表述：3! 之選擇 = **外部輸入**（Q1 之局部化；與 M7 docs/20 §七 u/d 指派
     NO-GO 同源）。
  N_trials = 6（全報）。零上游寫入。
""")
raise SystemExit(0)