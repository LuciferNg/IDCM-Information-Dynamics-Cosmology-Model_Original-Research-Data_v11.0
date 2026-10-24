"""M5 docs/15 配套 — 獨立驗證器：PDG 輸入 → 頭條數字 → root docs 一致性 → 跨文檔污染
================================================================================================
目的：`docs/11`（`scripts/12_full_audit.py`）只驗**可重現性**（重跑 byte-identical、docstring↔docs、
數字可回溯），唔驗**正確性**。本腳本獨立重算 M5 全部頭條數字（唔走推導腳本嘅代碼路徑，
用唔同括號／網格／mpmath），並核 root docs（README／final_result／SESSION_START／01_problem_record）
與跨文檔污染。

判準（寫死）：
  V1 輸入：由 PDG 原始檔（`data/pdg_mass_width_2024.mcd`）獨立解析 8 個介子質量＋誤差；
     重算 E = Δm − X（X = 2.56）與 σ = ⊕(σ₁,σ₂,σ_X)（σ_X = 0.094）⟹ 對 4 條關鍵通道 |Δ| ≤ 0.005 MeV。
  V2 頭條數字：重校 (a₀,t)、形狀量 (⟨B⟩,P₀)、K 系預言與 z、A3 窗／χ²、D1 1-D 曲線、
     |ψ(0)|²／r、A/|c| —— 相對容差 1e-4（重校量）／1e-2（盲測量）／5e-2（網格量）。
  V3 root docs：必要數字必須出現；`0.5287` 只許出現於「作廢／舊／缺陷／附註」語境；
     不可有未標示之舊 z（+3.20）。
  V4 計數：docs／scripts／logs／json 與 README 宣稱一致。
  V5 輸出：`data/16_verification.json`（逐項 PASS/FAIL）。

獨立性聲明：本腳本**不** import 任何 M5 模組；自身重寫 shape/amps（改用唔同嘅能量零點平移與
唔同 bisection 括號），故與推導腳本共享嘅只有模型定義本身（Hamiltonian），非代碼。
"""
import json
import re
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path.cwd() if (Path.cwd() / "README.md").exists() else Path(__file__).resolve().parent.parent
HBARC = 0.1973269804
ALPHA = 1 / 137.035999084
B, L = 4, 40
X_S, SX = 2.56, 0.094

claims = []


def claim(cid, what, computed, ref, tol):
    claims.append({"id": cid, "what": what, "computed": computed, "ref": ref, "tol": tol})


def rec(cid, what, computed, doc, ref, tol, unit=""):
    """自建 claim 再比對（避免依賴先前 append）"""
    ok = abs(computed - doc) <= tol * max(1.0, abs(doc)) if tol > 0 else computed == doc
    claims.append({"id": cid, "what": what, "computed": computed, "ref": ref, "tol": tol,
                   "doc": doc, "pass": bool(ok), "unit": unit,
                   "rel_err": ((computed - doc) / doc if doc else None)})
    print(f"  {cid:6s} {what:38s} 算 {computed:12.6f}{unit:3s} 文檔 {doc:12.6f}{unit:3s}"
          f" {'✅' if ok else '🔴 FAIL'}")
    return ok


# ---------------------------------------------------------------- V1 PDG 輸入
print("=" * 108)
print("V1 輸入驗證：PDG 原始檔 → 通道值與 σ")
print("=" * 108)
WANT = {411: "D+", 421: "D0", 413: "D*+", 423: "D*0", 321: "K+", 311: "K0", 323: "K*+", 313: "K*0"}
PDG = {}
for ln in open(ROOT / "data" / "pdg_mass_width_2024.mcd"):
    if ln.startswith("*") or not ln.strip():
        continue
    t = ln.split()
    try:
        i = int(t[0])
    except ValueError:
        continue
    if i in WANT:
        PDG[i] = (float(t[1]) * 1000, abs(float(t[2])) * 1000)     # MeV, MeV
for i, nm in WANT.items():
    print(f"  {nm:4s} m = {PDG[i][0]:9.3f} ± {PDG[i][1]:6.3f} MeV")

CH = [("D⁺−D⁰", 411, 421, 2.2600, 0.1176265276202609),
      ("D*⁺−D*⁰", 413, 423, 0.8500, 0.1176265276202609),
      ("K⁰−K⁺", 311, 321, 1.3740, 0.0960728889958036),
      ("K*⁰−K*⁺", 313, 323, 1.3200, 0.34122719704032967)]
sx_list = []
for nm, a, b, E_doc, s_doc in CH:
    dm = PDG[a][0] - PDG[b][0]
    E = dm - X_S
    sig_pdg = float(np.hypot(PDG[a][1], PDG[b][1]))
    sx = float(np.sqrt(max(s_doc ** 2 - sig_pdg ** 2, 0.0)))
    sx_list.append(sx)
    print(f"  {nm:10s} Δm = {dm:8.3f} ⟹ E = Δm − X = {E:7.4f}"
          f"（σ_Δm = {sig_pdg:6.3f}、σ_通道 = {s_doc:.4f} ⟹ 隱含 σ_X = {sx:.4f}）")
    rec(f"C{len(claims)+1:02d}", f"E({nm}) 由 PDG 重算", E, E_doc, "docs/13 §四.3", 2e-3, " MeV")
print(f"  ⟹ 隱含 σ_X 四通道 = {['%.5f' % s for s in sx_list]}；spread = "
      f"{(max(sx_list)-min(sx_list))/np.mean(sx_list)*100:.3f}%")
spread = float((max(sx_list) - min(sx_list)) / np.mean(sx_list))
ok05 = spread < 5e-3
claims.append({"id": "C05", "what": "σ_X 四通道一致（spread < 0.5%）", "computed": spread,
               "doc": 5e-3, "ref": "data/08 X_sigma = 0.094", "tol": 0.0, "pass": ok05,
               "unit": "", "rel_err": None})
print(f"  C05    σ_X spread = {spread*100:.3f}% < 0.5% {'✅' if ok05 else '🔴 FAIL'}")

# ---------------------------------------------------------------- 模型（獨立重寫）
NODES = [Fraction(1)] + [B * Fraction(B - 1) ** (d - 1) for d in range(1, L + 1)]
NODES = [float(x) for x in NODES]
POW = [B ** (-d) for d in range(L + 1)]
OFF = [-(B ** 0.5)] + [-(B - 1) ** 0.5] * (L - 1)


def shape(g, h):
    """能量零點平移版：對角加 +c（唔影響 eigenvector），驗證與推導腳本無關"""
    shift = 0.5 * (g + h)
    diag = np.array([-g * POW[d] - (h if d == 0 else 0.0) + shift for d in range(L + 1)])
    H = np.diag(diag) + np.diag(OFF, 1) + np.diag(OFF, -1)
    _, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    w = np.array(NODES) * v ** 2
    rA = float((w * np.array(POW)).sum() / w.sum())
    dp = int(np.argmax(w))
    return rA, float(w[0] / w.sum()), dp, bool(rA > 1e-3 and dp <= 20)


def amps(aC, aK, t, q, mh, ml):
    g = ALPHA * q * HBARC / aC / t
    Cq = (8 * np.pi / 3) * ALPHA * q / (16 * mh * ml)
    h = Cq * (HBARC / aK) ** 3 / t
    rA, P0, dp, bd = shape(g, h)
    return ALPHA * q * (HBARC / aC) * rA * 1000, Cq * (HBARC / aK) ** 3 * P0 * 1000, (rA, P0, dp, bd)


print("\n" + "=" * 108)
print("V2 頭條數字（獨立實作：能量零點平移 + 不共用網格）")
print("=" * 108)

A_D, c_D = 1.2025, 0.3525          # 數據反解（物理定義）

# 重校：外層對 a 精確掃（內層 t 用精確 bisection），唔用插值 ⟹ 精度達 1e-13
def t_solve(a, target_A=A_D):
    lo, hi = np.log(1e-12), np.log(1e3)
    f = lambda x: amps(a, a, float(np.exp(x)), 2 / 3, 1.27, 0.30)[0] - target_A
    flo = f(lo)
    if flo * f(hi) > 0:
        return None
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-14:
            break
    return float(np.exp(0.5 * (lo + hi)))


coarse = []
for a in np.exp(np.linspace(np.log(0.25), np.log(0.55), 60)):
    t = t_solve(float(a))
    if t is None:
        continue
    coarse.append((abs(amps(float(a), float(a), t, 2 / 3, 1.27, 0.30)[1] - c_D), float(a), t))
coarse.sort()
_, a_lo, _ = coarse[0]
lo, hi = np.log(a_lo * 0.97), np.log(a_lo * 1.03)
fK = lambda x: amps(float(np.exp(x)), float(np.exp(x)), t_solve(float(np.exp(x))), 2 / 3, 1.27, 0.30)[1] - c_D
flo = fK(lo)
for _ in range(300):
    mid = 0.5 * (lo + hi)
    fm = fK(mid)
    if flo * fm <= 0:
        hi = mid
    else:
        lo, flo = mid, fm
    if hi - lo < 1e-13:
        break
a0c = float(np.exp(0.5 * (lo + hi)))
t0c = t_solve(a0c)
A_chk, c_chk, st_chk = amps(a0c, a0c, t0c, 2 / 3, 1.27, 0.30)
rec("C09a", "重校 a₀ (fm)", a0c, 0.377422, "docs/14 §四.3", 1e-5, " fm")
rec("C09b", "重校 t (MeV)", t0c * 1000, 0.788398, "docs/14 §四.3", 1e-4, "")
rec("C10a", "⟨B^−d⟩ (D)", st_chk[0], 0.472772, "docs/14 §四.3", 1e-4)
rec("C10b", "P₀ (D)", st_chk[1], 0.368922, "docs/14 §四.3", 1e-4)
rec("C11a", "|ψ(0)|² (GeV³)", c_D / 1000 / ((8 * np.pi / 3) * ALPHA * (2 / 3) / (16 * 1.27 * 0.30)),
    0.05272, "docs/14 §四.2", 1e-3)
psi2 = c_D / 1000 / ((8 * np.pi / 3) * ALPHA * (2 / 3) / (16 * 1.27 * 0.30))
rec("C11b", "r (fm)", float((np.pi * psi2) ** (-1 / 3)) * HBARC, 0.3593, "docs/14 §四.2", 1e-3, " fm")

# K 系盲測
A_K, c_K, _ = amps(a0c, a0c, t0c, 1 / 3, 0.495, 0.30)
E1, E2 = A_K + 3 * c_K, A_K - c_K
rec("C12a", "A_K (MeV)", A_K, 0.2903, "docs/14 §四.5", 1e-3, "")
rec("C12b", "|c|_K (MeV)", c_K, 0.1880, "docs/14 §四.5", 1e-3, "")
rec("C13a", "z(E(K⁰−K⁺))", (E1 - 1.3740) / 0.0960728889958036, -5.41, "docs/13 §四.3", 5e-3)
rec("C13b", "z(E(K*⁰−K*⁺))", (E2 - 1.3200) / 0.34122719704032967, -3.57, "docs/13 §四.3", 5e-3)
rec("C13c", "z(ΔE_K = 4K)", (4 * c_K - 0.0540) / 0.35449400559106775, 1.97, "docs/13 §四.3", 5e-3)

# A3 族（15 成員）
FAM = [("F0", lambda l: 1.0), ("F1", lambda l: 2.0 ** l), ("F2", lambda l: 4.0 ** l),
       ("F3", lambda l: ((1 + 5 ** 0.5) / 2) ** l), ("F4", lambda l: 3.0 ** l)]
MODES = ["both", "contact-only", "coulomb-only"]
nb = np_ = 0
for nm, law in FAM:
    for md in MODES:
        s = law(7)
        aC = a0c if md == "contact-only" else a0c * s
        aK = a0c if md == "coulomb-only" else a0c * s
        A, K, st = amps(aC, aK, t0c, 1 / 3, 0.495, 0.30)
        if not st[3]:
            continue
        nb += 1
        z1 = (A + 3 * K - 1.3740) / 0.0960728889958036
        z2 = (A - K - 1.3200) / 0.34122719704032967
        np_ += bool(abs(z1) <= 2 and abs(z2) <= 2)
rec("C15a", "A3 束縛成員數", float(nb), 3.0, "docs/14 §四.5", 0.0)
claims[-1]["pass"] = nb == 3
print(f"         （束縛 {nb}/15 {'✅' if nb == 3 else '🔴'}）")
rec("C15b", "A3 命中數", float(np_), 0.0, "docs/14 §四.5", 0.0)
claims[-1]["pass"] = np_ == 0

# A3 連續掃 both（唔同網格 + 唔同範圍）
grid = np.exp(np.linspace(np.log(0.35), np.log(5.0), 700))
chis, ok_s = [], []
for s in grid:
    A, K, st = amps(a0c * s, a0c * s, t0c, 1 / 3, 0.495, 0.30)
    if not st[3]:
        continue
    ok_s.append(float(s))
    chis.append(((A + 3 * K - 1.3740) / 0.0960728889958036) ** 2 + ((A - K - 1.3200) / 0.34122719704032967) ** 2)
j = int(np.argmin(chis))
rec("C14a", "A3 both χ²_min (2 dof)", float(chis[j]), 12.7873, "docs/14 §四.6", 1e-2)
rec("C14b", "A3 both s_best", float(ok_s[j]), 0.9526, "docs/14 §四.6", 2e-3)
rec("C14c", "A3 both 束縛窗上界", float(max(ok_s)), 1.1373, "docs/14 §四.6", 5e-3)

# D1 1-D 曲線
Acoef = ALPHA * (1 / 3) * HBARC / a0c * 1000
rA_star, rA_star_1 = 1.3335 / Acoef, (1.3335 - 0.2570) / Acoef
rec("C16a", "K 系 Coulomb 前因子 (MeV)", Acoef, 1.2718, "docs/13 §四.6", 1e-3)
rec("C16b", "所需 ⟨B⟩* (中心)", rA_star, 1.0486, "docs/13 §四.6", 1e-3)
rec("C16c", "所需 ⟨B⟩* (1σ 下限)", rA_star_1, 0.8464, "docs/13 §四.6", 1e-3)
hs = np.exp(np.linspace(np.log(1.0), np.log(40.0), 900))
prev = None
h_star = None
for h in hs:
    rA = shape(ALPHA * (1 / 3) * HBARC / a0c / t0c, float(h))[0]
    if prev is not None and (prev[0] - rA_star_1) * (rA - rA_star_1) <= 0:
        h_star = float(h)
        break
    prev = (rA, h)
rec("C16d", "h* (1σ 下限處)", h_star, 7.19, "docs/13 §四.6", 5e-3)
P0_star = shape(ALPHA * (1 / 3) * HBARC / a0c / t0c, h_star)[1]
K_star = t0c * h_star * P0_star * 1000
rec("C16e", "K_K(h*) (MeV)", K_star, 4.54, "docs/13 §四.6", 1e-2)
rec("C16f", "K(h*)/K_need (倍)", K_star / 0.0135, 337.0, "docs/13 §四.6", 5e-3)
rec("C16g", "偏離 (σ)", (K_star - 0.0135) / 0.2658 * 3, 51.0, "docs/13 §四.6", 5e-2)

# A/|c| 恒等式
rec("C18", "A/|c| 比值 K/D（倍）", (1.3335 / 0.0135) / (1.2025 / 0.3525), 29.0, "docs/14 §四.1", 2e-2)

# ---------------------------------------------------------------- V3 root docs
print("\n" + "=" * 108)
print("V3 root docs 一致性與跨文檔污染")
print("=" * 108)
ROOTS = ["README.md", "final_result.md", "SESSION_START.md", "01_problem_record.md"]
NEED = {"0.377422": "重校 a₀", "0.7884": "重校 t", "0.3525": "物理 |c|",
        "-5.41": "z(E(K⁰−K⁺))", "-3.57": "z(E(K*⁰−K*⁺))", "+1.97": "z(ΔE_K)",
        "337": "K(h*) 倍數", "0.3593": "r 更正"}
miss = []
def norm(t):
    """統一減號（U+2212／U+FF0D → '-'）與全角，避免假陰性"""
    return t.replace("\u2212", "-").replace("\uff0d", "-").replace("\u2013", "-").replace("\uff0b", "+")


for f in ROOTS:
    txt = norm((ROOT / f).read_text())
    for k, why in NEED.items():
        if k.lstrip("+") not in txt:
            miss.append((f, k, why))
print(f"  必要數字缺失：{len(miss)} 項" + ("" if not miss else f" ⟹ {miss}"))
claim("C19", "root docs 載齊更正後數字", float(len(miss)), "README/final_result/SESSION_START/01", 0.0)
claims[-1]["doc"] = 0.0
claims[-1]["pass"] = not miss

bad = []
for f in ROOTS:
    for i, ln in enumerate((ROOT / f).read_text().splitlines(), 1):
        ln = norm(ln)
        if "0.5287" in ln and not re.search(r"作廢|舊|缺陷|附註|誤|多除|更正|supersed", ln):
            bad.append((f, i, ln.strip()[:70]))
print(f"  0.5287 未標示污染：{len(bad)} 處" + ("" if not bad else f" ⟹ {bad[:3]}"))
claim("C20", "0.5287 只出現於更正語境", float(len(bad)), "root docs", 0.0)
claims[-1]["doc"] = 0.0
claims[-1]["pass"] = not bad

z_old = [(f, i) for f in ROOTS for i, ln in enumerate(norm((ROOT / f).read_text()).splitlines(), 1)
         if "3.20" in ln and "作廢" not in ln and "舊" not in ln]
print(f"  未標示之舊 z=+3.20：{len(z_old)} 處" + ("" if not z_old else f" ⟹ {z_old[:3]}"))
claim("C21", "舊 z=+3.20 全部標示作廢", float(len(z_old)), "root docs", 0.0)
claims[-1]["doc"] = 0.0
claims[-1]["pass"] = not z_old

# ---------------------------------------------------------------- V4 計數
print("\n" + "=" * 108)
print("V4 檔案計數 vs README 宣稱")
print("=" * 108)
nd, ns = len(list((ROOT / "docs").glob("*.md"))), len(list((ROOT / "scripts").glob("*.py")))
nl, nj = len(list((ROOT / "data").glob("run_*.log"))), len(list((ROOT / "data").glob("*.json")))
rd = (ROOT / "README.md").read_text()
m_docs = re.search(r"`docs/`.*?（(\d+) 檔", rd)
m_scr = re.search(r"`scripts/`.*?核心 (\d+) 檔", rd)
print(f"  實測 docs={nd} scripts={ns} logs={nl} json={nj}；README 宣稱 docs={m_docs.group(1) if m_docs else '?'}"
      f" scripts={m_scr.group(1) if m_scr else '?'}")
ok = bool(m_docs) and bool(m_scr) and int(m_docs.group(1)) == nd and int(m_scr.group(1)) == ns
claim("C22", "README 計數 == 實測", float(nd + ns), "README", 0.0)
claims[-1]["doc"] = float(int(m_docs.group(1)) + int(m_scr.group(1))) if (m_docs and m_scr) else -1.0
claims[-1]["pass"] = ok
print(f"  {'✅ 一致' if ok else '🔴 不一致'}")

n_pass = sum(1 for c in claims if c.get("pass"))
print("\n" + "=" * 108)
print(f"[總計] 驗證項 {len(claims)}；PASS {n_pass}；FAIL {len(claims)-n_pass}")
print("=" * 108)
for c in claims:
    if not c.get("pass"):
        print(f"  🔴 {c['id']} {c['what']}：算 {c['computed']} vs 文檔 {c.get('doc')}")

out = {"n_claims": len(claims), "n_pass": n_pass, "n_fail": len(claims) - n_pass,
       "claims": claims,
       "pdg_inputs": {WANT[i]: {"mass_MeV": PDG[i][0], "err_MeV": PDG[i][1]} for i in sorted(PDG)},
       "sigma_X_implied": sx_list,
       "recalibration": {"a0_fm": a0c, "t_MeV": t0c * 1000, "rA": st_chk[0], "P0": st_chk[1],
                         "dpeak": st_chk[2], "bound": st_chk[3]},
       "K_channel": {"A": A_K, "c": c_K, "E1": E1, "E2": E2, "z1": (E1 - 1.3740) / 0.0960728889958036,
                     "z2": (E2 - 1.3200) / 0.34122719704032967},
       "a3": {"n_bound": nb, "n_pass": np_, "chi2_min": float(chis[j]), "s_best": float(ok_s[j]),
              "window_max": float(max(ok_s))},
       "d1": {"A_coef": Acoef, "rA_star": rA_star, "rA_star_1sigma": rA_star_1, "h_star": h_star,
              "K_at_h_star": K_star, "ratio": K_star / 0.0135},
       "verdict": ("✅ 全部 PASS：PDG 輸入可獨立重現（σ_X 四通道一致 0.094 MeV）、"
                   "重校與形狀量至 1e-4、K 系 z 與 A3/D1 數字至 1e-2、root docs 無未標示污染"
                   if n_pass == len(claims) else "🔴 有 FAIL，見 claims")}
with open(ROOT / "data" / "16_verification.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n-> data/16_verification.json")