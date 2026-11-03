#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M13 scripts/02_criteria_T1_T5.py — 判準本體（先立判準，後算）

本腳本**不產生任何新候選**。它只做三件事：
  A. 把 SESSION_START／01_problem_record §三 之 T1–T5 ＋ §4.2 兩定理，寫成可執行判準（門檻全部硬編碼、跑前凍結）。
  B. 以 M12 已判死之三個核作**自檢**：pipeline 必須逐位重現 M12 之精確有理值
     （`1/γ → 42/41, −2/41`；`P+ → −63/41, 43/41`；`P− → −28/41, +28/41`）。
     重現成功 = harness 與 M12 之 `δ_x = Σ_n w_n φ_n(k_x)Δ_n` 係同一對象 ⟹ 後續候選判定可通。
  C. 印出 T5 可裁決性門檻（σ 需降至何值）—— 呢個係 Q2 之定量答案。

判準凍結（跑前寫死，不得於結果後調整）：
  G0 出處閘      : builder 源碼之識別符 ⊆ 註冊白名單（AST 檢查），且須列 provenance。
  G1 符號閘      : 需 ρ_u < 0；若核恆正（min w > 0）且 ρ_u ≥ 0 ⟹ FAIL。
  G2 味盲閘      : |1 − ρ_d| > 3σ_ρd（常數項壓倒 ⟹ ρ→1 ⟹ 不能帶味）。
  T1 非可分解    : ① 源碼須引用 K_tree 側符號；② 對 K 參數之響應族秩 ≥ 2（秩 1 且 ∝ w ⟹ 因子相消）。
  T2 離散性      : 對 builder 內每個連續常數 ±1e-6 相對擾動，|Δρ_d| ≤ 1e-9（否則 dof > 0）。
                   並須聲明 N_trials（離散選擇乘積），全報不得只報命中。
  T3 非 ℤ₂-odd   : max|w_n + w_{45−n}| / max|w_n| > 1e-30（=0 ⟹ 奇宇稱 ⟹ ρ_u = −ρ_d 定理鎖死）。
  T4 零參數      : builder 源碼不得引用目標常數（AST 名檢 'TARGET'、'target'、'rho_d_target'）；
                   且候選須單一預先指定；若 N_trials > 1，須報 look-elsewhere 校正。
  T5 可裁決      : 報 σ_ρd/格距、半格 vs σ、以及 5σ 判別所需 σ。σ ≤ 半格 = 可裁決（現況不可能）。

用法：本檔可被 03+ 之候選腳本 import（`from importlib import ...` 或 sys.path），
      或直接執行以重印自檢與門檻。
"""
import ast
import sys
import mpmath as mp

mp.mp.dps = 40

# ---------------- 註冊常數（白名單；來源見 01_problem_record §4.2／§六） ----------------
N = 44                      # 模式數（= N_d/2 = 88/2）
ND = 88                     # N_d
SUM_GAMMA = 90              # Σγ_even（/90 慣例）
B = 4                       # B = |A₅|/15，鏈分支
NGEN = 3                    # N_gen = B − 1
KTREE = 33                  # K_tree 深度
MU2 = 2                     # μ²
H11, H21 = 36, 98           # CY₃ Hodge 數
SEC88 = 1 / mp.cos(mp.pi / 88)   # sec(π/88)

WHITELIST = {
    "gamma", "g", "n", "N", "ND", "SUM_GAMMA", "B", "NGEN", "KTREE", "MU2",
    "H11", "H21", "SEC88", "M0", "phi", "Delta", "V", "A5", "sin", "cos", "tan",
    "pi", "sqrt", "mp", "sum", "range", "len", "abs", "min", "max", "enumerate",
    "float", "int", "list", "zip", "sorted", "round", "fsum", "power", "exp", "log",
}
FORBIDDEN_TARGET = {"TARGET", "target", "rho_d_target", "RHO_D", "rho_target", "obs"}

# ---------------- 已測量目標（只可事後比較；T4 靠 AST 封鎖 builder 引用） ----------------
RHO_D_OBS = mp.mpf("0.40350365")
RHO_U_OBS = mp.mpf("-0.02588837")
SIG_RD = mp.mpf("0.044785178")      # 正典（R27／docs/42；docs/38 之 0.03510 已 SUPERSEDED）
SIG_RU = mp.mpf("0.033583899")
GRID = mp.mpf(1) / 41               # 載體自然格距
HALF = GRID / 2
SIG_GATE = {"rho_d": SIG_RD, "rho_u": SIG_RU}
N_GRID_NEIGH = 2                    # 最近兩個格點

# ---------------- 幾何（φ_n(j) = √(2/45) sin(π n j / 45)） ----------------
def phi(n, j):
    return mp.sqrt(mp.mpf(2) / 45) * mp.sin(mp.pi * n * j / 45)

DELTA = [phi(n, 3) - phi(n, 2) for n in range(1, N + 1)]     # Δ_n = φ_n(3) − φ_n(2)
K_STATION = {"u": 2, "d": 3, "s": 4}                          # 已裁站指派（docs/14／docs/28）
V = {x: [phi(n, k) * DELTA[n - 1] for n in range(1, N + 1)] for x, k in K_STATION.items()}


def gamma(n, conv="/90"):
    """γ_n：/90 配 45−n（M12 慣例）；/88 配 44−n（M6 慣例）。"""
    den = 90 if conv == "/90" else 88
    return 4 * mp.sin(mp.pi * n / den) ** 2


def rho_of(w):
    """ρ_x = ⟨w,v_x⟩/⟨w,v_s⟩（v_x = φ(k_x)Δ）。"""
    f = lambda x: mp.fsum(w[n - 1] * V[x][n - 1] for n in range(1, N + 1))
    ds = f("s")
    return {"rho_u": f("u") / ds, "rho_d": f("d") / ds,
            "delta_u": f("u"), "delta_d": f("d"), "delta_s": ds}


def z_of(got, exp, sig):
    return (mp.mpf(got) - mp.mpf(exp)) / mp.mpf(sig)


# ---------------- 判準（每個 gate 回 (PASS/FAIL, 說明)） ----------------
def gate_G0(cand):
    src = cand.get("builder_src", "")
    names = {node.id for node in ast.walk(ast.parse(src))
             if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)}
    bad = {x for x in names if x not in WHITELIST and not x.startswith("_")}
    missing = [p for p in cand.get("provenance", []) if not p]
    if bad:
        return False, f"越界識別符 {sorted(bad)}（須全部可追溯至註冊常數）"
    if not cand.get("provenance"):
        return False, "無 provenance 清單"
    if missing:
        return False, "provenance 有空項"
    return True, f"識別符 {len(names)} 全部在冊；provenance {len(cand['provenance'])} 條"


def gate_G1(cand, r):
    if r["rho_u"] < 0:
        return True, f"ρ_u = {mp.nstr(r['rho_u'],6)} < 0 ✅"
    if min(cand["w"]) > 0:
        return False, f"核恆正 ⟹ 範疇性 ρ_u ≥ 0（實得 {mp.nstr(r['rho_u'],6)}）"
    return False, f"ρ_u = {mp.nstr(r['rho_u'],6)} ≥ 0（與目標異號）"


def gate_G2(cand, r):
    dev = abs(1 - r["rho_d"])
    if dev <= 3 * SIG_RD:
        return False, f"|1−ρ_d| = {mp.nstr(dev,4)} ≤ 3σ ⟹ 味盲極限（常數項壓倒）"
    return True, f"|1−ρ_d| = {mp.nstr(dev,4)} > 3σ ✅"


def gate_T1(cand):
    src = cand.get("builder_src", "")
    names = {n.id for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Name)}
    tree_syms = {"KTREE", "H11", "H21", "M0", "A5", "B", "NGEN"}
    if not (names & tree_syms):
        return False, "源碼無 K_tree 側符號 ⟹ 純鏈因子 ⟹ 於比值相消（docs/31）"
    fam = cand.get("k_family")
    if fam is None:
        return False, "未提供 K 參數響應族（無法驗可分性）"
    M = mp.matrix(len(fam), N)
    for i, w in enumerate(fam):
        for j in range(N):
            M[i, j] = w[j]
    # 秩檢定（相對門檻）
    s = mp.svd_r(M, compute_uv=False)
    tol = s[0] * mp.mpf("1e-12") if s else 0
    rank = sum(1 for x in s if x > tol)
    if rank < 2:
        return False, f"響應族秩 = {rank} ⟹ 可寫成 a(γ)·b(K) 單乘積 ⟹ 因子相消"
    return True, f"響應族秩 = {rank} ≥ 2 ⟹ 非 Kronecker-可分解 ✅"


def gate_T2(cand):
    if cand.get("continuous_params"):
        return False, f"含連續參數 {cand['continuous_params']} ⟹ dof > 0"
    p = cand.get("dof_probe")
    if p is None:
        return False, "未提供 dof 擾動探針（無法證明 dof = 0）"
    worst = max(abs(d) for d in p)
    if worst > mp.mpf("1e-9"):
        return False, f"±1e-6 相對擾動引起 |Δρ_d| = {mp.nstr(worst,4)} > 1e-9 ⟹ 有隱含連續參數"
    ntrials = cand.get("n_discrete_choices")
    if not isinstance(ntrials, int) or ntrials < 1:
        return False, f"N_trials 未寫死（得 {ntrials}）"
    return True, f"dof = 0（worst |Δρ| = {mp.nstr(worst,3)}）；N_trials = {ntrials}"


def gate_T3(cand):
    w = cand["w"]
    odd = max(abs(w[n - 1] + w[45 - n - 1]) for n in range(1, 45))
    scale = max(abs(x) for x in w)
    rat = odd / scale
    if rat < mp.mpf("1e-30"):
        return False, "w_n + w_(45−n) ≡ 0 ⟹ ℤ₂-odd ⟹ 奇宇稱定理 ⟹ ρ_u = −ρ_d 精確（死）"
    return True, f"oddness = {mp.nstr(rat,6)} > 1e-30 ⟹ 非 ℤ₂-odd ✅"


def gate_T4(cand):
    src = cand.get("builder_src", "")
    names = {n.id for n in ast.walk(ast.parse(src))
             if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    hit = names & FORBIDDEN_TARGET
    if hit:
        return False, f"builder 引用目標 {sorted(hit)} ⟹ 反解 ⟹ 擬合"
    if cand.get("post_hoc_variants", 1) > 1 and not cand.get("look_elsewhere_corrected"):
        return False, "N > 1 變體而未做 look-elsewhere 校正"
    return True, "builder 未見目標；單一預先指定 ✅"


def gate_T5(cand, r):
    """可裁決性 = 候選值與最近結構競爭者之距 ÷ σ > 5。
    結構競爭者集 = ℚ/41 格點（載體自然分母，docs/38 §D）。
    若落格點 ⟹ 距最近格點 ≪ σ ⟹ 必然 1σ 內 ⟹ 零資訊 ⟹ 不可裁決。"""
    ratio = SIG_RD / GRID
    need_sig = GRID / 5
    k = mp.nint(r["rho_d"] * 41)
    nearest = mp.mpf(k) / 41
    dist_sig = abs(r["rho_d"] - nearest) / SIG_RD
    if dist_sig < 5 and mp.nint(r["rho_d"] * 41) == mp.nint(RHO_D_OBS * 41):
        return False, (f"落格點 {k}/41（距 {mp.nstr(dist_sig,3)}σ < 5σ）⟹ 不可裁決；"
                       f"σ_ρd/格距 = {mp.nstr(ratio,6)}，需 σ ≤ {mp.nstr(need_sig,4)}（1/{mp.nstr(SIG_RD/need_sig,3)}）")
    if dist_sig < 5:
        return False, (f"與最近格點 {k}/41 相距 {mp.nstr(dist_sig,3)}σ < 5σ ⟹ 不可裁決；"
                       f"σ_ρd/格距 = {mp.nstr(ratio,6)}")
    return True, (f"與最近格點 {k}/41 相距 {mp.nstr(dist_sig,3)}σ ≥ 5σ ⟹ 可裁決（當前數據即已判）；"
                  f"σ_ρd/格距 = {mp.nstr(ratio,6)}")


GATES = [("G0", gate_G0), ("G1", gate_G1), ("G2", gate_G2), ("T1", gate_T1),
         ("T2", gate_T2), ("T3", gate_T3), ("T4", gate_T4)]


def judge(cand):
    r = rho_of(cand["w"])
    out = {}
    for name, fn in GATES:
        if name in ("G1", "G2"):
            ok, msg = fn(cand, r)
        else:
            ok, msg = fn(cand)
        out[name] = (ok, msg)
    ok5, msg5 = gate_T5(cand, r)
    out["T5"] = (ok5, msg5)
    zd = z_of(r["rho_d"], RHO_D_OBS, SIG_RD)
    zu = z_of(r["rho_u"], RHO_U_OBS, SIG_RU)
    alive = all(v[0] for k, v in out.items() if k != "T5")
    return {"rho": r, "gates": out, "z_d": zd, "z_u": zu, "alive": alive}


# ================= B. harness 自檢（M12 已判死三核） =================
def selftest():
    print("=== B. harness 自檢：重現 M12 R25 之精確有理值 ===")
    cases = {
        "1/γ（控制）":      dict(w=[1 / gamma(n) for n in range(1, N + 1)]),
        "P+ 偶":           dict(w=[1 / gamma(n) + 1 / gamma(45 - n) for n in range(1, N + 1)]),
        "P− 奇":           dict(w=[1 / gamma(n) - 1 / gamma(45 - n) for n in range(1, N + 1)]),
    }
    expect = {"1/γ（控制）": (mp.mpf(-2) / 41, mp.mpf(42) / 41),
              "P+ 偶": (mp.mpf(43) / 41, mp.mpf(-63) / 41),
              "P− 奇": (mp.mpf(28) / 41, mp.mpf(-28) / 41)}
    allok = True
    for k, v in cases.items():
        r = rho_of(v["w"])
        eu, ed = expect[k]
        du = abs(r["rho_u"] - eu)
        dd = abs(r["rho_d"] - ed)
        ok = du < mp.mpf("1e-25") and dd < mp.mpf("1e-25")
        allok &= ok
        print(f"  {k:12s} ρ_u={mp.nstr(r['rho_u'],10):>13s} ρ_d={mp.nstr(r['rho_d'],10):>13s} "
              f"|Δ|={mp.nstr(max(du,dd),3):>9s} {'PASS' if ok else 'FAIL'}")
        print(f"      ×41: ρ_u={mp.nstr(r['rho_u']*41,6)}  ρ_d={mp.nstr(r['rho_d']*41,6)}")
    print(f"  ⟹ pipeline ≡ M12 之 δ_x = Σ_n w_n φ_n(k_x)Δ_n：{'PASS' if allok else 'FAIL'}")

    print("\n=== C. gate 行為自檢（已知判死核）===")
    demo = [
        dict(name="1/γ", w=cases["1/γ（控制）"]["w"],
             builder_src="w = [1/gamma(n) for n in range(1, N+1)] * KTREE",
             provenance=["K_chain γ_n (/90)", "K_tree(標度, 因相消)"],
             k_family=[[1 / gamma(n) * c for n in range(1, N + 1)] for c in (1, 2)],
             dof_probe=[mp.mpf(0)], n_discrete_choices=1, continuous_params=[],
             post_hoc_variants=1),
        dict(name="P− 奇", w=cases["P− 奇"]["w"],
             builder_src="w = [1/gamma(n) - 1/gamma(45-n) for n in range(1, N+1)] * KTREE",
             provenance=["K_chain γ_n", "half-integer parity"],
             k_family=[[(1 / gamma(n) - 1 / gamma(45 - n)) * c for n in range(1, N + 1)] for c in (1, 2)],
             dof_probe=[mp.mpf(0)], n_discrete_choices=1, continuous_params=[],
             post_hoc_variants=1),
    ]
    for c in demo:
        res = judge(c)
        print(f"  [{c['name']}] ρ_d={mp.nstr(res['rho']['rho_d'],8)}  z_d={mp.nstr(res['z_d'],4)}  "
              f"z_u={mp.nstr(res['z_u'],4)}  alive={res['alive']}")
        for g, (ok, msg) in res["gates"].items():
            print(f"      {g:3s} {'PASS' if ok else 'FAIL'}  {msg}")
    print("\n  ★ 預期行為：1/γ 之 T1 應 FAIL（族秩 1 ⟹ 可分離 ⟹ 相消）；P− 之 T3 應 FAIL（奇宇稱）。")

    print("\n=== D. T5 可裁決性門檻（Q2 定量答案）===")
    print(f"  現行：σ_ρd = {mp.nstr(SIG_RD,8)}，σ_ρu = {mp.nstr(SIG_RU,8)}")
    print(f"  格距 1/41 = {mp.nstr(GRID,8)}；半格 = {mp.nstr(HALF,8)}")
    print(f"  σ_ρd ÷ 格距 = {mp.nstr(SIG_RD/GRID,6)}；σ_ρu ÷ 格距 = {mp.nstr(SIG_RU/GRID,6)}")
    print(f"  單一格點判別（5σ）需 σ_ρd ≤ 格距/5 = {mp.nstr(GRID/5,8)} ⟹ 需改善 {mp.nstr(SIG_RD/(GRID/5),4)}×")
    print(f"  半格門檻需改善 {mp.nstr(SIG_RD/HALF,4)}×；close 觸發器 (σ ≲ 0.005) 需改善 {mp.nstr(SIG_RD/mp.mpf('0.005'),4)}×")
    print("  ⟹ 結論：現精度下任何有理/格點命中皆不可裁決（100% 落在某格點 1σ 內）；")
    print("     唯一可裁決路徑 = 候選值被結構強制為單點，且未來測量達上述門檻。")
    return allok


if __name__ == "__main__":
    ok = selftest()
    print("\n=== E. 判準狀態 ===")
    print(f"  判準凍結版本: 2026-10-29（M13 解閘首輪）；自檢 {'PASS' if ok else 'FAIL'}")
    print("  本腳本不作任何候選判定、不改上游檔。N_trials = 0（判準自檢）。")
    sys.exit(0 if ok else 1)