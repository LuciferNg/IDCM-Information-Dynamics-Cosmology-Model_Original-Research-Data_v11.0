#!/usr/bin/python3
"""19 — handoff 回覆之算術核（精確 Fraction；浮點只作顯示）
對象：(a) M12 任務書收窄之 existence 回執；(b) M14 H1 之 0 bit 回執；
      (c) M15 兩條來件（Q2 已閉？／3 與 3/16 分開）
只讀 PDG vendored 檔（經 01 之同一解析器邏輯重算關鍵量），零 import 其他腳本。
"""
from fractions import Fraction as F
import json
import re
import sys
from pathlib import Path

mp = __import__("mpmath")
mp.mp.dps = 40

ROOT = Path(__file__).resolve().parent.parent
PDG = ROOT.parent / "m12_cy3_flavour_operator/data/pdg2024/mass_width_2024.txt"
B = 4
N_d = 44
N_gen = 3

out = []


def say(s=""):
    print(s)
    out.append(s)


def hdr(s):
    say()
    say("=" * 78)
    say(s)
    say("=" * 78)


def _row(pid):
    for line in PDG.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("*") or len(line) < 70:
            continue
        cols = [line[0:8], line[8:16], line[16:24], line[24:32]]
        if any(c.strip() == pid for c in cols) and line[33:51].strip():
            return line
    raise SystemExit(f"PDG 缺 {pid}")


def mass(pid, name=""):
    return mp.mpf(_row(pid)[33:51].strip()) * 1000


# ---------------------------------------------------------------- M15 ① 幾何級數
hdr("M15 ① 幾何級數 Σ B^{−k} = B/(B−1)")
for k in (0, 1, 4, 11):
    part = sum(F(1, B ** j) for j in range(k + 1))
    say(f"  k=0..{k:2d} 部分和 = {part}  = {mp.nstr(mp.mpf(part.numerator)/part.denominator, 12)}")
p12 = sum(F(1, B ** j) for j in range(12))
tail = F(B, B - 1) - p12
say(f"  部分和 k=0..11 = {p12}（M15 報 5592405/4194304 {'✓' if p12 == F(5592405, 4194304) else '✗'}）")
say(f"  餘項           = {tail}（M15 報 1/12582912 {'✓' if tail == F(1, 12582912) else '✗'}）")
say(f"  閉式           = {F(B, B-1)} = {mp.nstr(mp.mpf(4)/3, 12)}")
say("  ★ 判：恆等式成立，但對**任何** B 皆成立（見下），單獨不構成 4/3 之來源")
say("     B:  " + "  ".join(f"{b}→{F(b, b-1)}" for b in (2, 3, 4, 5, 6, 8)))

# ---------------------------------------------------------------- M15 ② 44/33
hdr("M15 ② 鏈／樹模數比 N_d/K_tree = 44/33")
say(f"  44/33 = {F(44, 33)} = {mp.nstr(mp.mpf(44)/33, 12)}")
say(f"  而 N_d·(B−1)/B = 44·3/4 = {F(N_d*(B-1), B)} ⟹ 33 恰等於 N_d(B−1)/B")
say(f"  (44/33) ÷ (B/(B−1)) = {F(44,33)/F(B,B-1)}  ⟹ 代數上恆等（同一條 1/(1−1/B)）")
say("  ★ 判：若 K_tree=33 之出處 = 「N_d 去掉 1/B」，則 ② ≡ ①，唔係第二條獨立來源")
say("     反之若 33 有獨立登錄（樹層模數直數），則 ② 為真約束 ⟹ 需 M15 指認")

# ---------------------------------------------------------------- 迹通道對 B 之鑑別力
hdr("迹通道 2Δ_N/Δ_Ξ = (B−1)/B² 對 B 之鑑別力（M13 自有）")
m_n = mass("2112", "MASS")
m_p = mass("2212", "MASS")
m_Sm = mass("3112", "MASS")
m_S0 = mass("3212", "MASS")
m_Sp = mass("3222", "")
m_Xm = mass("3312", "MASS")
m_X0 = mass("3322", "MASS")
say(f"  質量 (MeV): n={mp.nstr(m_n,12)} p={mp.nstr(m_p,12)} Σ⁻={mp.nstr(m_Sm,12)} "
    f"Σ⁰={mp.nstr(m_S0,12)} Σ⁺={mp.nstr(m_Sp,12)} Ξ⁻={mp.nstr(m_Xm,12)} Ξ⁰={mp.nstr(m_X0,12)}")
dN = m_n - m_p
dX = m_Xm - m_X0
tr = 2 * dN / dX
say(f"  Δ_N = {mp.nstr(dN,12)}  Δ_Ξ = {mp.nstr(dX,12)}")
say(f"  trace = 2Δ_N/Δ_Ξ = {mp.nstr(tr,12)}")
# σ（沿用 l2 之傳播：σ(m) 由 PDG 誤差欄）
def sdev(pid):
    for line in PDG.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("*") or len(line) < 70:
            continue
        cols = [line[0:8], line[8:16], line[16:24], line[24:32]]
        if any(c.strip() == pid for c in cols):
            s = line[51:60].strip()          # 第一欄 = 對稱誤差（+X −X 之前者）
            return mp.mpf(s) * 1000 if s else mp.mpf(0)
    return mp.mpf(0)
s_tr = mp.sqrt((2 * sdev("2112") / dX) ** 2 + (2 * sdev("2212") / dX) ** 2
               + (2 * dN * sdev("3312") / dX ** 2) ** 2 + (2 * dN * sdev("3322") / dX ** 2) ** 2)
say(f"  σ(trace) = {mp.nstr(s_tr,8)}")
q = tr / 2                      # Δ_N/Δ_Ξ（= trace 之半）
s_q = s_tr / 2
say(f"  Δ_N/Δ_Ξ = {mp.nstr(q,12)}  σ = {mp.nstr(s_q,8)}   （trace = 2q）")
say(f"  對 3/16:  z = {mp.nstr((q - mp.mpf(3)/16)/s_q,6)}")
say("  B 之整數鑑別（以 (B−1)/B² 比 Δ_N/Δ_Ξ）：")
for b in (2, 3, 4, 5, 6):
    v = mp.mpf(b - 1) / b ** 2
    say(f"    B={b}: (B−1)/B² = {mp.nstr(v,8):>10s}   z = {mp.nstr((q - v)/s_q, 6):>9s}")
say("  ★ 判：迹通道**會**鑑別 B（B=5 已 >4σ）；M15 ① 之 B/(B−1) 形式則見下")

# ---------------------------------------------------------------- α 通道對 B 之鑑別力
hdr("α 組合 α_s/(α_u+α_d)：與迹之恆等 ＋ 對 B 之鑑別力")
s_A = mp.sqrt(sdev("3112") ** 2 + sdev("3212") ** 2)          # σ(Δ_Σ)
s_X = mp.sqrt(sdev("3312") ** 2 + sdev("3322") ** 2)          # σ(Δ_Ξ)
rho_d_prop = 2 * (m_Sm - m_S0) / dX - 1
sig_rho_d = mp.sqrt((2 * s_A / dX) ** 2 + (2 * (m_Sm - m_S0) * s_X / dX ** 2) ** 2)
say(f"  核對：ρ_d = 2Δ_Σ/Δ_Ξ − 1 = {mp.nstr(rho_d_prop,12)}（canonical 0.403503649635）")
say(f"        σ_ρd = {mp.nstr(sig_rho_d,12)}（canonical 0.04478517819）")
m_Sm_ = m_Sm
dS = m_Sm_ - m_S0                      # Δ_Σ = m(Σ⁻) − m(Σ⁰)
d_s = dX / 2
d_d = dS - d_s
d_u = dN - d_d
rho_u = d_u / d_s
rho_d = d_d / d_s
say(f"  δ_s = Δ_Ξ/2 = {mp.nstr(d_s,10)}  δ_d = {mp.nstr(d_d,10)}  δ_u = {mp.nstr(d_u,10)}")
say(f"  ρ_u = {mp.nstr(rho_u,12)}   ρ_d = {mp.nstr(rho_d,12)}   trace = {mp.nstr(rho_u+rho_d,12)}")
# x = α_u/α_s, y = α_d/α_s；ρ_u = x − y/2, ρ_d = y − x/2 之精確反解
xx = (2 * rho_d + 4 * rho_u) / 3
yy = (4 * rho_d + 2 * rho_u) / 3
say(f"  反解 x=α_u/α_s = {mp.nstr(xx,12)}   y=α_d/α_s = {mp.nstr(yy,12)}")
say(f"  回代：ρ_u = {mp.nstr(xx - yy/2,12)}（差 {mp.nstr(xx-yy/2-rho_u,3)}）、"
    f"ρ_d = {mp.nstr(yy - xx/2,12)}（差 {mp.nstr(yy-xx/2-rho_d,3)}）")
alpha_sum = 1 / (xx + yy)
trace_ = rho_u + rho_d
say(f"  ★ 恆等：α_s/(α_u+α_d) = {mp.nstr(alpha_sum,12)}   1/(2·trace) = {mp.nstr(1/(2*trace_),12)}"
    f"   差 = {mp.nstr(alpha_sum - 1/(2*trace_),3)}")
say("    ⟹ 此「α 組合」與迹係**同一條**量，非獨立通道")
r_rho = 3 / (2 * rho_d_prop)           # ρ_u=0 分支下之 r（= docs/05 之 r）
sig_r = 3 * sig_rho_d / (2 * rho_d_prop ** 2)
say(f"  （參考）ρ_u=0 分支：r = 3/(2ρ_d) = {mp.nstr(r_rho,10)} ± {mp.nstr(sig_r,6)}；"
    f"精確反解 r = 1/x = {mp.nstr(1/xx,10)} —— 兩者不同，因實測 ρ_u≠0")
sig_as = s_tr / (2 * trace_ ** 2)
say(f"  σ(α_s/(α_u+α_d)) = {mp.nstr(sig_as,8)}（由 σ(trace) 傳播）")
say("  family A：α_s/(α_u+α_d) vs B/(B−1)（≡ trace vs (B−1)/(2B)）")
for b in (2, 3, 4, 5, 6):
    v = mp.mpf(b) / (b - 1)
    say(f"    B={b}: {mp.nstr(v,8):>9s}   z = {mp.nstr((alpha_sum - v)/sig_as, 6):>9s}")
say("  ★ 判：B=4 最接近（|z|≈0.22）；B=5 只 1.8σ，未達 3σ ⟹ 此**弱族**可偏好 B=4 但不能排除 B=5")
say("      （前次回覆用 r/3 比較並報『B=5 更近』—— 該算法混用 ρ_u=0 分支與 ρ_d 之 σ，已廢）")

# ---------------------------------------------------------------- γ_max−1 vs B−1
hdr("M15 促請：3 = γ_max−1 與迹 law 之 3/16 分開")
gamma_max = 4                      # γ44 = 4 = B
say(f"  γ_max = γ44 = {gamma_max}；γ_max−1 = {gamma_max-1}（O215 docs/03 T7 之獨立譜源，k-free）")
say(f"  N_gen = B−1 = {B-1}")
say(f"  同值：{'是' if gamma_max-1 == N_gen else '否'}；同源：否（一為譜端點，一為 B 之計數）")
say(f"  迹 law 之常數：3/16 = (B−1)/B² = {F(B-1, B*B)}")
say("  ★ 判：接受 M15 之促請 —— 迹 law 內之 3 一律標 N_gen=B−1，不得引 γ_max−1")

# ---------------------------------------------------------------- M12 任務書 existence test
hdr("M12 任務書收窄之 existence 回執（兩選一）")
say("  (a) K≥5 離散結構：K=5 類大小 34,752,256；雙 1σ 空容量預期 10,773；0.25σ 695")
say("      ⟹ 任何命中皆屬空容量 ⟹ **判死**（docs/02、docs/03）")
say("  (b) 3 模式 + 非 ±1 權重：三模耗 1 連續 dof（T2 不滿足）⟹ 收斂為單一標量")
say(f"      r = α_s/α_u = {mp.nstr(r_rho,10)} ± {mp.nstr(sig_r,6)}；窗 [3.34606,4.18155] 內唯一整數 4")
say("      ⟹ 權重比由 ρ_u=0 固定為 α=(1,2,4)（A₂ Cartan 零方向），但**次序 provenance 未閉**")
say("  ⟹ existence test 結果：**未能產出結構算子**；證據上限 1 bit（通道相依）")

# ---------------------------------------------------------------- M14 H1 回執
hdr("M14 H1 回執（α 之來歷 = 0 bit）")
say("  M13 標籤層定理（docs/11）：α 只依 k ⟹ ρ_d=3/4 強制 +7.74σ；α 只依 d ⟹ ρ_d=0 強制 −9.01σ；")
say("  任意 h(k,d) ⟹ dof=2。216 個 (k,d) 函數族 **0/216** 落窗 ⟹ 與 M14 之 0 bit 一致。")
say("  ⟹ M13 不再掃 (k,d)→α；殘餘（次序即指派）已發還 M14 Q2。")

# ---------------------------------------------------------------- 匯總
hdr("§匯總")
say("  ① 幾何級數：恆等成立、B-blind（∀B 同一式）⟹ 單獨 0 bit")
say("  ② 44/33：與 ① 代數恆等，除非 K_tree=33 有獨立登錄 ⟹ 需 M15 指認")
say("  ③ 兩族鑑別力不同（同一數據）：family A（B/(B−1)，≡ trace vs (B−1)/(2B)）弱 —— B=4 −0.225σ、B=5 +1.809σ；")
say("     family B（(B−1)/B² vs Δ_N/Δ_Ξ）強 —— B=5 +4.93σ、B=3 −5.72σ。α 組合 = 1/(2·trace) 恆等 ⟹ **非**獨立通道")
say("  ④ 3 vs 3/16：接受分開（N_gen=B−1 用於迹 law；γ_max−1 屬 O215 T7，不互引）")
say("  ⑤ M12 任務書：兩選一皆未能產出結構算子 ⟹ existence test 誠實報負")
say("  ⑥ M14 H1：M13 0/216 與 M14 0 bit 同向 ⟹ 殘餘已發還")

RESULT = {"r": str(r_rho), "sigma_r": str(sig_r), "trace": str(tr), "sigma_trace": str(s_tr),
          "geo_partial_12": str(p12), "geo_tail": str(tail), "ratio_44_33": str(F(44, 33)),
          "z_trace_B4": float((q - mp.mpf(3) / 16) / s_q),
          "z_trace_B5": float((q - mp.mpf(4) / 25) / s_q),
          "alpha_sum_ratio": str(alpha_sum),
          "sigma_alpha_sum_ratio": str(sig_as),
          "r_rho_branch": str(r_rho),
          "r_exact": str(1 / xx),
          "z_alpha_B4": float((alpha_sum - mp.mpf(4) / 3) / sig_as),
          "z_alpha_B5": float((alpha_sum - mp.mpf(5) / 4) / sig_as),
          "z_trace_B3": float((q - mp.mpf(2) / 9) / s_q),
          "z_familyB_B5": float((q - mp.mpf(4) / 25) / s_q)}
say()
say("JSON → data/run_19_handoff_replies.json")
(ROOT / "data/run_19_handoff_replies.json").write_text(
    json.dumps(RESULT, ensure_ascii=False, indent=2), encoding="utf-8")
(ROOT / "data/run_19_handoff_replies.log").write_text("\n".join(out) + "\n", encoding="utf-8")
sys.exit(0)