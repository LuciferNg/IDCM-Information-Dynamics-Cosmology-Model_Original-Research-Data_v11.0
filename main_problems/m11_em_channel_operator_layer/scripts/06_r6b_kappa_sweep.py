"""M11 scripts/06 — R6b：κ 敏感度掃描（Q11 之決策面）
================================================================================================
問題（Q11）：M5 之 κ 有兩支流通（κ=0 主校準 a₀=0.377422 vs κ=0.215 天花板點 a=0.243274），
            何者正典未裁。M11 之**全部判詞**都經 k₀ = α·ħc/a₀ 依賴 a₀ ⟹ 必須先量出「判詞之 κ 穩健性」。

方法（零新參數、零重推、零改 M5）：
  1. **逐字移植** M5 `scripts/21_R_ceiling.py` 之模型（amps：41 層截斷鏈 ＋ 接觸項；w ∝ N_nodes v²）
     與校準方程：`A_D = A(a,t) + κ·g(7)·q_C`、`c(a,t) = C_D`（A_D=1.2025、C_D=0.3525）。
  2. 解 (a,t)：以 **κ=0 為起點之連續延拓**（前一 κ 之解作初值，`scipy.optimize.root` on log-log），
     殘差 |ΔA|+|Δc| < 1e−9（M5 docs/19 R1 同判準）；失敗則多起點重試。
  3. **移植驗證**：與 M5 已發表之 `data/21_R_ceiling.json`（κ ∈ [0.18,0.32]）逐點比對 (a,t)。
  4. 對每個 κ 重算 M11 之判詞：
     R1（C9，8 通道）：t'(κ) = em/k₀(κ)，em ≡ t_ref·k₀_ref（em 不含 k₀ ⟹ 純縮放）；
     R4（四模型盲測）：t'(κ) = em/k₀(κ)，feasible ⟺ t' ∈ [box_min, box_max]。
  5. 物理閘（用 M5 自己之判準）：接觸非負（κ < A_D/(g(0)q_D) = 0.24599，`docs/34`）；
     束縛守門（D、K 兩族 P₀ > 1e−3 ∧ argmax w ≤ 20，`docs/24` R1）。

判準（事前寫死）：
  K1  t 只經 k₀ 依賴 κ ⟹ 判詞係 a(κ) 之階梯函數。
  K2  κ* ＝使任一判詞改變之最小 κ（網格上）。
  K3  **穩健** ⟺ 在「物理閘全過」之 κ 區間內無 κ*。
  K4  禁改 M5 之校準／M11 之箱／禁為翻轉而選 κ。
輸出：`data/06_r6b_kappa_sweep.json`
"""
import json
import os

import numpy as np
from scipy.optimize import root

ALPHA = 1 / 137.035999084
HBARC = 0.1973269804                    # GeV·fm（與 M5 scripts/21 一致；amps ×1000 → MeV）
B, L = 4, 40
NODES = [1.0] + [B * (B - 1) ** (d - 1) for d in range(1, L + 1)]
POW = [B ** (-d) for d in range(L + 1)]
OFF = np.array([-np.sqrt(B)] + [-np.sqrt(B - 1)] * (L - 1))
M_C, M_S, M_L, Q_C, Q_S = 1.27, 0.495, 0.30, 2 / 3, 1 / 3
G0, G7 = 2581 / 352, 1447 / 352
A_D, C_D = 1.2025, 0.3525
KAPPA_CONTACT_MAX = A_D / (G0 * Q_C)    # M5 docs/34：0.24599…
# K 之兩觀測量（M5 scripts/21：E1K／E2K／σ 逐字）
E1K, E2K = 1.3740, 1.3200
S1K, S2K = 0.0960728889958036, 0.34122719704032967
_O_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # O 根（cwd 無關）
_DATA = os.path.join(_O_DIR, "data")
M5_JSON = os.path.join(os.path.dirname(_O_DIR), "m5_hadron_bound_state_spectrum", "data", "21_R_ceiling.json")

LOG = []


def emit(s=""):
    print(s, flush=True)
    LOG.append(s)


def amps(a, t, q, mh, ml):
    """M5 scripts/21 amps() 逐字移植（A、c 單位 MeV；bd = 束縛守門）"""
    gcoef = ALPHA * q * HBARC / a / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / a) ** 3 / t
    diag = np.array([-gcoef * POW[d] - (h if d == 0 else 0.0) for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    w = w / w.sum()
    rA, P0 = float((w * np.array(POW)).sum()), float(w[0])
    bd = bool(P0 > 1e-3 and int(np.argmax(w)) <= 20)
    return ALPHA * q * (HBARC / a) * rA * 1000, Cq * (HBARC / a) ** 3 * P0 * 1000, bd


def _F(x, kp):
    a, t = float(np.exp(x[0])), float(np.exp(x[1]))
    if not (1e-4 < a < 30 and 1e-13 < t < 1e4):
        return [1e3, 1e3]
    A, c, _ = amps(a, t, Q_C, M_C, M_L)
    return [np.log(max(A + kp * G7 * Q_C, 1e-300) / A_D), np.log(max(c, 1e-300) / C_D)]


def solve_kappa(kp, guess=None, multistart=True):
    """解 D 校準；回傳 (a_fm, t_MeV, resid) 或 None"""
    starts = []
    if guess is not None:
        starts.append([np.log(guess[0]), np.log(guess[1] / 1000)])
    if multistart:
        starts += [[-1.0, 0.0], [-0.5, -7.0], [-2.0, 1.0], [-1.5, -2.0], [0.5, -5.0], [-3.0, -3.0]]
    starts = starts[:1] if (guess is not None and not multistart) else starts
    best = None
    for g in starts:
        try:
            s = root(_F, g, args=(kp,), method="hybr", tol=1e-11)
        except Exception:
            continue
        if not s.success:
            continue
        r = float(np.max(np.abs(_F(s.x, kp))))
        if r > 1e-9:
            continue
        a, t = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
        if not (1e-3 < a < 20 and 1e-12 < t < 1e3):
            continue
        A, c, _ = amps(a, t, Q_C, M_C, M_L)
        if abs(A + kp * G7 * Q_C - A_D) + abs(c - C_D) > 1e-8:
            continue
        if best is None or r < best[2]:
            best = (a, t, r)
    return best


def verdicts_at(k0, C9, K0_REF, D4):
    """回傳 (r1_nbad, r1_detail, r4 dict)"""
    r1, nbad = {}, 0
    for ch, v in C9.items():
        em = v["t"] * K0_REF
        tp = em / k0
        lo, hi = min(v["box_min"], v["box_max"]), max(v["box_min"], v["box_max"])
        feas = bool(lo - 1e-9 <= tp <= hi + 1e-9)
        nbad += int(not feas)
        r1[ch] = dict(t=tp, feasible=feas)
    r4 = {}
    for m in D4["models"]:
        nfail, det = 0, {}
        for ch, v in m["verdicts"].items():
            if v.get("O") is None or v.get("em") is None:
                continue
            tp = v["em"] / k0
            lo, hi = min(v["box"]), max(v["box"])
            feas = bool(lo - 1e-9 <= tp <= hi + 1e-9)
            nfail += int(not feas)
            det[ch] = dict(t=tp, feasible=feas)
        r4[m["label"][:1]] = dict(blind_fails=nfail, detail=det,
                                  grade=("OUT（盲測否證）" if nfail else m["grade"]))
    return nbad, r1, r4


def main():
    D1 = json.load(open(os.path.join(_DATA, "01_r1_criteria.json")))
    D4 = json.load(open(os.path.join(_DATA, "04_r4_necessity.json")))
    K0_REF = D1["C9"]["k0_MeV"]
    A0_REF = HBARC * ALPHA / (K0_REF / 1000)
    C9 = D1["C9"]["channels"]

    emit("=" * 108)
    emit("M11 scripts/06 — R6b：κ 敏感度掃描（Q11 之決策面）")
    emit("=" * 108)
    emit("模型：M5 scripts/21（amps／A_D=1.2025／C_D=0.3525／g(7)=1447/352／g(0)=2581/352）逐字移植")
    emit(f"參考點：data/01 k₀ = {K0_REF:.6f} MeV ⟹ a₀ = {A0_REF:.8f} fm（κ=0 應重現此值 ⟹ 移植自測）")
    emit(f"物理閘：接觸非負 κ < A_D/(g(0)q_D) = {KAPPA_CONTACT_MAX:.5f}（M5 docs/34）；"
         "束縛守門 P₀>1e−3 ∧ argmax≤20（docs/24 R1）")

    s0 = solve_kappa(0.0)
    emit(f"\n[A] 移植自測 κ=0：a = {s0[0]:.8f} fm、t = {s0[1]*1000:.5f} MeV（殘差 {s0[2]:.1e}）"
         f"｜M5 docs/13 正典 a₀ = 0.377422 fm、t = 0.7884 MeV ⟹ Δa = {abs(s0[0]/0.377422-1)*100:.5f}%")

    if os.path.exists(M5_JSON):
        ref = json.load(open(M5_JSON))["rows"]
        m5 = {round(float(r["kappa"]), 4): (float(r["a0_fm"]), float(r["t_MeV"])) for r in ref}
        dev = []
        for kp in sorted(m5):
            s = solve_kappa(kp, guess=(0.3774, 0.7884), multistart=False) or solve_kappa(kp)
            a5, t5 = m5[kp]
            dev.append((kp, None if s is None else abs(s[0] / a5 - 1),
                        None if s is None else abs(s[1] * 1000 / t5 - 1)))
        bad = [d for d in dev if d[1] is None or d[1] > 1e-6 or d[2] > 1e-6]
        mx = max((d[1] for d in dev if d[1] is not None), default=0.0)
        mx2 = max((d[2] for d in dev if d[2] is not None), default=0.0)
        emit(f"\n[B] 移植驗證 vs M5 data/21_R_ceiling.json：{len(dev)} 點，"
             f"max|Δa/a| = {mx:.2e}、max|Δt/t| = {mx2:.2e} ⟹ "
             + ("✅ 逐點一致（同一支、同一模型）" if not bad else f"⚠ {len(bad)} 點不符：{bad[:3]}"))

    grid = [0.0]
    kp, guess, failed = 0.0, (s0[0], s0[1] * 1000), []
    KSCAN_MAX = 0.3300        # 超出接觸界，僅為重現 M5 docs/19 T5（χ² 最小 @0.3）作 z 通道驗證
    while kp < KSCAN_MAX - 1e-12:
        kp = round(min(kp + 0.0025, KSCAN_MAX), 6)
        s = solve_kappa(kp, guess=guess, multistart=False) or solve_kappa(kp)
        if s is None:
            failed.append(kp)
            continue
        guess = (s[0], s[1] * 1000)
        grid.append(kp)
    for extra in (0.1825, 0.215, 0.26986):
        if all(abs(extra - g) > 1e-9 for g in grid):
            grid.append(extra)
    grid = sorted(grid)
    emit(f"\n[C] 掃描 {len(grid)} 個 κ 點（κ ∈ [0, {max(grid):.5f}]；步 0.0025 ＋ κ=0.1825／0.215／0.26986）"
         + (f"；{len(failed)} 點無解" if failed else ""))

    rows, guess = [], (s0[0], s0[1] * 1000)
    for kp in grid:
        s = solve_kappa(kp, guess=guess, multistart=False) or solve_kappa(kp)
        if s is None:
            continue
        a, t, res = s
        guess = (a, t * 1000)
        k0 = ALPHA * HBARC / a * 1000
        _, _, bd_D = amps(a, t, Q_C, M_C, M_L)
        A_K, c_K, bd_K = amps(a, t, Q_S, M_S, M_L)
        contact_ok = (A_D - kp * G7 * Q_C) > 0
        nbad, r1, r4 = verdicts_at(k0, C9, K0_REF, D4)
        AK_tot = A_K + kp * G0 * Q_S
        z1 = (AK_tot + 3 * c_K - E1K) / S1K      # M5 docs/19 T5 同式
        z2 = (AK_tot - c_K - E2K) / S2K
        rows.append(dict(kappa=kp, a0_fm=a, t_MeV=t * 1000, k0_MeV=k0, resid=res,
                         bound_D=bd_D, bound_K=bd_K, contact_ok=contact_ok,
                         A_K=AK_tot, c_K=c_K, z1=z1, z2=z2, chi2=z1 * z1 + z2 * z2,
                         r1_nbad=nbad, r1=r1, r4=r4,
                         gate=bool(contact_ok and bd_D and bd_K)))

    for r in rows:
        fails = "/".join(f"{k}:{v['blind_fails']}" for k, v in r["r4"].items())
        emit(f"  κ={r['kappa']:.5f}  a={r['a0_fm']:.5f}  t={r['t_MeV']:.3f}  k₀={r['k0_MeV']:.4f}  "
             f"D束={int(r['bound_D'])} K束={int(r['bound_K'])} 接觸≥0={int(r['contact_ok'])} "
             f"閘={'✅' if r['gate'] else '✗'}  R1壞={r['r1_nbad']}/8  "
             f"K之z=({r['z1']:+.3f},{r['z2']:+.3f}) χ²={r['chi2']:6.2f}  盲測失敗 {fails}")

    base = next(r for r in rows if abs(r["kappa"]) < 1e-12)
    emit("\n[D] 文獻點與翻轉點")
    for tgt in (0.1825, 0.215, 0.26986):
        r = next((x for x in rows if abs(x["kappa"] - tgt) < 1e-9), None)
        if r:
            emit(f"  κ={tgt}: a={r['a0_fm']:.5f} fm、k₀={r['k0_MeV']:.4f} MeV、R1 壞 {r['r1_nbad']}/8、"
                 "盲測失敗 " + "/".join(f"{k}:{v['blind_fails']}" for k, v in r["r4"].items())
                 + f"、閘={'✅' if r['gate'] else '✗'}（D束={int(r['bound_D'])} K束={int(r['bound_K'])} "
                 f"接觸≥0={int(r['contact_ok'])}）")
    flips = []
    for r in rows:
        ch = []
        if r["r1_nbad"] != base["r1_nbad"]:
            ch.append(f"R1:{base['r1_nbad']}→{r['r1_nbad']}")
        for k, v in r["r4"].items():
            if v["blind_fails"] != base["r4"][k]["blind_fails"]:
                ch.append(f"R4-{k}:{base['r4'][k]['blind_fails']}→{v['blind_fails']}")
        if ch:
            flips.append(dict(kappa=r["kappa"], change=ch, gate=r["gate"]))
    if flips:
        for f in flips[:16]:
            emit(f"  ⚠ κ={f['kappa']:.5f}：{'；'.join(f['change'])}｜閘={'✅' if f['gate'] else '✗'}")
        if len(flips) > 16:
            emit(f"  …另有 {len(flips)-16} 個變動點")
    else:
        emit("  ✅ 全掃描範圍內判詞不變（R1 可行集合、四模型盲測失敗數皆同）")

    adm = [r for r in rows if r["gate"]]
    segs = []
    for r in adm:
        if segs and abs(r["kappa"] - segs[-1][1]) < 0.003:
            segs[-1][1] = r["kappa"]
        else:
            segs.append([r["kappa"], r["kappa"]])
    emit(f"\n[E] 物理閘全過之 κ：{len(adm)}/{len(rows)} 點"
         + ("，區間聯集 " + " ∪ ".join(f"[{a:.5f}, {b:.5f}]" for a, b in segs) if adm else "（空集）"))
    if adm:
        emit(f"    ⟹ a₀ ∈ [{min(r['a0_fm'] for r in adm):.5f}, {max(r['a0_fm'] for r in adm):.5f}] fm；"
             f"k₀ ∈ [{min(r['k0_MeV'] for r in adm):.4f}, {max(r['k0_MeV'] for r in adm):.4f}] MeV")
        emit(f"    ⟹ 區內 R1 壞 {sorted({r['r1_nbad'] for r in adm})}/8；"
             f"盲測失敗 C={sorted({r['r4']['C']['blind_fails'] for r in adm})}、"
             f"D={sorted({r['r4']['D']['blind_fails'] for r in adm})}")
        lo, hi = min(r["kappa"] for r in adm), max(r["kappa"] for r in adm)
        fin = [f for f in flips if lo <= f["kappa"] <= hi]
        emit(f"    ⟹ 閘內翻轉點：{'**無** ⟹ 判詞對 κ 穩健（Q11 對 M11 失效）' if not fin else str([f['kappa'] for f in fin]) + ' ⟹ 判詞不穩健，Q11 必須裁'}")
    else:
        emit("    🔴 閘交集為空 ⟹ M5 之 κ 族在物理上無可行點（另案）")

    # ── [F] κ 可否由 M5 之第三族（K）獨立釘死 ──
    emit("\n[F] κ 之資料釘死（用 M5 已有之 K 兩數；M11 不參與 ⟹ 對 M11 為盲）")
    for kp in (0.0, 0.0375, 0.09, 0.1825, 0.215, 0.26986):
        r = next((x for x in rows if abs(x["kappa"] - kp) < 1e-9), None)
        if r:
            emit(f"  κ={kp:<7} K: z₁={r['z1']:+.3f} z₂={r['z2']:+.3f} χ²={r['chi2']:6.2f} "
                 f"K束={int(r['bound_K'])} ⟹ {'✅ 合格' if r['gate'] else '✗ 閘外'}")
    gk = [r for r in adm]
    if gk:
        b = min(gk, key=lambda r: r["chi2"])
        emit(f"  ⟹ 閘內（D+K 皆束縛）之 K-χ² 最小：κ = {b['kappa']:.5f}（χ² = {b['chi2']:.3f}、"
             f"a = {b['a0_fm']:.5f} fm、k₀ = {b['k0_MeV']:.4f} MeV、√χ² = {b['chi2']**0.5:.2f}σ）")
        emit(f"     該點 M11 判詞：R1 壞 {b['r1_nbad']}/8、盲測失敗 "
             + "/".join(f"{k}:{v['blind_fails']}" for k, v in b["r4"].items()))
        b2 = min(rows, key=lambda r: r["chi2"])
        emit(f"  ⟹ 全掃描之 K-χ² 最小：κ = {b2['kappa']:.5f}（χ² = {b2['chi2']:.3f}、"
             f"✓χ² = {b2['chi2']**0.5:.2f}σ）、K束={int(b2['bound_K'])} ⟹ "
             f"{'**與束縛守門相容**（即 M5 docs/19 T5 之最优点，M5 docs/24 亦記該點 d_peak=2 束縛）' if b2['bound_K'] else '⚠ 與束縛守門不相容'}")
        band = [r for r in rows if r["chi2"] <= b2["chi2"] + 1.0]
        if band:
            emit(f"      Δχ² ≤ 1 帶：κ ∈ [{min(r['kappa'] for r in band):.5f}, {max(r['kappa'] for r in band):.5f}]"
                 f" ⟹ κ 由 K 兩數釘死至 ±{(max(r['kappa'] for r in band)-min(r['kappa'] for r in band))/2:.4f}"
                 f"（該帶內 a₀ ∈ [{min(r['a0_fm'] for r in band):.5f}, {max(r['a0_fm'] for r in band):.5f}] fm、"
                 f"k₀ ∈ [{min(r['k0_MeV'] for r in band):.3f}, {max(r['k0_MeV'] for r in band):.3f}] MeV）")
            emit(f"      該帶 M11 判詞：R1 壞 {sorted({r['r1_nbad'] for r in band})}/8；"
                 f"盲測失敗 A={sorted({r['r4']['A']['blind_fails'] for r in band})}"
                 f" B={sorted({r['r4']['B']['blind_fails'] for r in band})}"
                 f" C={sorted({r['r4']['C']['blind_fails'] for r in band})}"
                 f" D={sorted({r['r4']['D']['blind_fails'] for r in band})}")

    # ── [G] 各通道之精確翻轉點 κ*（可行 ⟺ t ∈ 箱；t(κ) 單調 ⟹ 唯一根）──
    emit("\n[G] 各通道之精確翻轉點 κ*（t(κ) = em/k₀(κ) 單調遞減 ⟹ 唯一）")
    targets = []
    for ch, v in C9.items():
        em = v["t"] * K0_REF
        lo, hi = min(v["box_min"], v["box_max"]), max(v["box_min"], v["box_max"])
        t_ref = em / K0_REF
        dist = max(lo - t_ref, t_ref - hi, 0.0)          # 出箱距離（箱內 = 0）
        targets.append((f"R1 {ch}", em, lo, hi, dist))
    for m in D4["models"]:
        for ch, v in m["verdicts"].items():
            if v.get("O") is None or v.get("em") is None:
                continue
            lo, hi = min(v["box"]), max(v["box"])
            t_ref = v["em"] / K0_REF
            dist = max(lo - t_ref, t_ref - hi, 0.0)
            targets.append((f"R4-{m['label'][:1]} {ch}", v["em"], lo, hi, dist))

    _cache = {}

    def a_of(kp):
        if kp in _cache:
            return _cache[kp]
        s = solve_kappa(kp, guess=(0.3774, 0.7884), multistart=False) or solve_kappa(kp)
        _cache[kp] = s[0] if s else None
        return _cache[kp]

    def kstar_for_scale(s_tgt):
        """解 a(κ) = αħc/(k₀*) ，其中 k₀* = K0_REF/s_tgt（t = t_ref·K0_REF/k₀ ⟹ t = bnd ⟺ s = s_tgt）"""
        a_star = HBARC * ALPHA / ((K0_REF / s_tgt) / 1000)
        lo_k, hi_k = 0.0, 0.33
        f_lo = a_of(lo_k) - a_star
        f_hi = a_of(hi_k) - a_star
        if f_lo * f_hi > 0:
            return None
        for _ in range(44):
            mid = 0.5 * (lo_k + hi_k)
            f_m = a_of(mid) - a_star
            if f_m * f_lo <= 0:
                hi_k = mid
            else:
                lo_k, f_lo = mid, f_m
        return 0.5 * (lo_k + hi_k)

    def t_at(kp, t_ref):
        s = solve_kappa(kp, guess=(0.3774, 0.7884), multistart=False)
        return None if s is None else t_ref * K0_REF / (ALPHA * HBARC / s[0] * 1000)

    for lab, em, lo, hi, dist in sorted(targets, key=lambda z: z[4]):
        t_ref = em / K0_REF
        feas0 = bool(lo - 1e-9 <= t_ref <= hi + 1e-9)
        cand = []
        for s_tgt in (lo / t_ref, hi / t_ref):
            if not (0 < s_tgt <= 1 - 1e-12):
                continue
            ks = kstar_for_scale(s_tgt)
            if ks is not None and ks > 1e-6:
                cand.append((ks, s_tgt))
        flip = None
        for ks, _ in sorted(cand):
            tv = t_at(ks * 1.002 + 1e-6, t_ref)
            if tv is None:
                continue
            if bool(lo - 1e-9 <= tv <= hi + 1e-9) != feas0:
                flip = ks
                break
        pct = (dist / max(abs(lo), abs(hi), 1e-9)) * 100 if dist > 0 else 0.0
        emit(f"  {lab:22s} t_ref={t_ref:+8.4f} 箱[{lo:+.4f},{hi:+.4f}] "
             f"{'箱外 ' + format(pct, '6.3f') + '%' if dist > 0 else '箱內'} "
             f"⟹ 翻轉 κ* = " + ("（掃描內不翻）" if flip is None else f"{flip:.6f}"))

    out = dict(source="M5 scripts/21_R_ceiling.py（逐字移植 amps/校準方程/判準）",
               kappa_contact_max=KAPPA_CONTACT_MAX, k0_ref=K0_REF, a0_ref=A0_REF,
               selfcheck=dict(a_fm=s0[0], t_MeV=s0[1] * 1000, resid=s0[2]),
               rows=rows, flips=flips, admissible=[r["kappa"] for r in adm])
    with open(os.path.join(_DATA, "06_r6b_kappa_sweep.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2, default=str)
    emit("\n-> data/06_r6b_kappa_sweep.json")


main()