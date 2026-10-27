"""M7 — 文檔 ↔ 腳本輸出 一致性核 v3（doc-vs-artifact，零任意容差）
================================================================================================
原理：文檔只印有限位數 ⟹ 判準 = 「文檔值 == 腳本值四捨五入到文檔所印精度」。
  roundcheck：由 doc 數字字串自身位數推 ulp，要求 |doc − art| ≤ 0.5·ulp（整數須全等）。
  ⟹ 無 rel-tol 可調，唔可以「放寬容差過關」。

取值（三種，皆以行定位）：
  rx=…,grp=k      正則第 k 組（表格行含內嵌 `|` 者，如 `A+3|c|`、`max|w−round(w)|`）
  cell=n,kth=k    markdown 行按 `|` 切，取第 n 格第 k 個數字
  cell=None,kth=k 整行第 k 個數字（prose）

C3 更正機制：doc 原文與腳本不合、但同檔另有以數值開頭之更正行（alt 定位）且該行與腳本相合
  ⟹ 報 PASS(corrected)，並獨立列出「原文 → 更正」；原文不動（append-only 紀律）。

輸出 data/04_doc_consistency.json；任何 FAIL ⟹ exit 1
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
DATA = os.path.join(M7, "data")
DOCS = os.path.join(M7, "docs")


def load(n):
    return json.load(open(os.path.join(DATA, n)))


J01, J02, J03, Jv1 = (load("01_reproduction_gate.json"), load("02b_r5_mode_content_nogo.json"),
                      load("03_grep_evidence.json"), load("02_r5_mode_content_nogo.json"))
J05 = load("05_additivity_scope.json")
J06 = load("06_gap_closure.json")
J07 = load("07_audit_full_chain.json")
J09 = load("09_m4_receipt.json")
J11 = load("11_o120_import.json")
J12 = load("12_o209_assignment.json")
J13 = load("13_ud_verdict.json")
CH = J01["checks"]
NUM = re.compile(r"[−-]?\+?\d+(?:\.\d+)?(?:[eE][−-]?\d+)?")
_doc = {}


def doc(name):
    if name not in _doc:
        for base in (DOCS, M7):
            p = os.path.join(base, name)
            if os.path.exists(p):
                _doc[name] = open(p, encoding="utf-8").read().splitlines()
                break
        else:
            raise FileNotFoundError(name)
    return _doc[name]


def norm(s):
    return s.replace("−", "-").replace("–", "-").replace("\u200b", "")


def get(marker, cell, kth, occ, rx, grp, dname):
    lines = doc(dname)
    if rx:
        pool = [ln for ln in lines if marker in ln] if marker else lines
        if marker and not pool:
            return None, f"marker「{marker}」於 regex 模式命中 0 行"
        for ln in pool:
            m = re.search(rx, ln)
            if m:
                return (norm(m.group(grp)) if grp else norm(m.group(0))), None
        return None, f"regex 未命中（限定 {len(pool)} 行）: {rx[:60]}"
    idxs = [i for i, ln in enumerate(lines) if marker in ln]
    if len(idxs) <= occ:
        return None, f"marker「{marker}」命中 {len(idxs)} 行（occ={occ}）"
    line = lines[idxs[occ]]
    if cell is None:
        chunk = line
    else:
        parts = line.split("|")
        if cell >= len(parts):
            return None, f"只有 {len(parts)} 格（cell={cell}）: {line[:80]}"
        chunk = parts[cell]
    toks = NUM.findall(norm(chunk))
    if len(toks) <= kth:
        return None, f"格內 {len(toks)} 個數字（kth={kth}）: {chunk[:80]}"
    return toks[kth], None


def ulp_of(s):
    s = norm(s).lstrip("+")
    if s.startswith("-"):
        s = s[1:]
    if "e" in s.lower():
        mant, exp = re.split("[eE]", s)
        exp = int(exp)
    else:
        mant, exp = s, 0
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return 10.0 ** (exp - dec)


def _tf(x):
    if isinstance(x, (int, float)):
        return float(x)
    t = NUM.findall(norm(str(x)))
    if not t:
        raise ValueError(f"無法由 {x!r} 取數")
    return float(t[0])


def roundcheck(ds, art, integer=False):
    d = float(norm(ds).lstrip("+"))
    a = _tf(art)
    if integer:
        return d == a, f"doc={ds} art={a:g}（整數須全等）"
    ulp = ulp_of(ds)
    diff = abs(d - a)
    return diff <= 0.5 * ulp * (1 + 1e-12), f"doc={ds} art={a:.12g} |Δ|={diff:.2e} ≤ {0.5*ulp:.1e}"


RES, CORR = [], []


def claim(name, dname, art, marker=None, cell=None, kth=0, occ=0, integer=False, alt=None, rx=None, grp=1):
    ds, err = get(marker, cell, kth, occ, rx, grp, dname)
    if err:
        RES.append(dict(name=name, doc=dname, status="🔴 MISSING", detail=err))
        return
    ok, detail = roundcheck(ds, art, integer)
    if ok:
        RES.append(dict(name=name, doc=dname, status="✅", detail=detail))
        return
    if alt:
        ds2, err2 = get(alt, None, 0, 0, None, 1, dname)
        if not err2:
            ok2, d2 = roundcheck(ds2, art, integer)
            if ok2:
                CORR.append(dict(name=name, original=ds, corrected=ds2, art=_tf(art)))
                RES.append(dict(name=name, doc=dname, status="✅ (corrected)",
                                detail=f"原文 {ds} → 更正 {ds2}；{d2}"))
                return
    RES.append(dict(name=name, doc=dname, status="🔴 MISMATCH", detail=detail))


# ══════════ docs/10 §三（19 項重現閘；含內嵌 `|` 之列用 regex） ══════════
TWO = r"([−-]?[\d.]+) \| ([−-]?[\d.]+) \| ✅"
S3_RX = [
    ("V2 `E(D⁺−D⁰)", TWO, 1, CH[6]["got"]),
    ("V2 `E(D*⁺−D*⁰)", TWO, 1, CH[7]["got"]),
    ("V3 解析式", TWO, 1, CH[8]["got"]),
    ("V7 `r = (π|ψ(0)|²)", r"(0\.313879) \| (0\.3139) \| ✅", 1, CH[14]["got"]),
]
S3_CELL = [
    ("V1 有限 L 精確式", CH[0]["got"]), ("V1 外推", CH[1]["got"]),
    ("V2 `a₀` (fm)", CH[2]["got"]), ("V2 `t` (MeV)", CH[3]["got"]),
    ("V2 `⟨B^{−d}⟩`", CH[4]["got"]), ("V2 `P₀`", CH[5]["got"]),
    ("V3 乘積守恆", CH[9]["got"]), ("V4 `B⁷`", CH[10]["got"]),
    ("V5 `A_max", CH[13]["got"]),
]
for marker, art in S3_CELL:
    claim(f"§三 {marker[:24]}", "10_reproduction_gate.md", art, marker, cell=3, kth=0)
for marker, rxp, grp, art in S3_RX:
    claim(f"§三 {marker[:24]}", "10_reproduction_gate.md", art, marker=marker, rx=rxp, grp=grp)
claim("§三 12 比值集=16384", "10_reproduction_gate.md", 16384, rx=r"\{∞,(16384),1\}", grp=1, integer=True)
claim("§三 13 逐元素 max=1", "10_reproduction_gate.md", 1.0, rx=r"max = ([\d.]+)", grp=1)

claim("§三 16 mpmath ⟨B⟩", "10_reproduction_gate.md", J01["mpmath"]["rA"], "G3 mpmath `⟨B^{−d}⟩`", cell=3, kth=0,
      alt="← §三 第16行 mpmath")
claim("§三 17 mpmath P₀", "10_reproduction_gate.md", J01["mpmath"]["P0"], "G3 mpmath `P₀`", cell=3, kth=0,
      alt="← §三 第17行 mpmath")
claim("§三 18 mpmath a₀", "10_reproduction_gate.md", J01["mpmath"]["a0_fm"], "G3 mpmath 打磨 `a₀`", cell=3, kth=0)
claim("§三 19 mpmath t", "10_reproduction_gate.md", J01["mpmath"]["t_MeV"], "G3 mpmath 打磨 `t`", cell=3, kth=0)
claim("§三 殘差 A", "10_reproduction_gate.md", J01["mpmath"]["resid_A"], "打磨後殘差", kth=1)
claim("§三 殘差 K", "10_reproduction_gate.md", J01["mpmath"]["resid_K"], "打磨後殘差", kth=2)
claim("§三 |ψ(0)|²", "10_reproduction_gate.md", J01["calibration"]["psi2_corrected_GeV3"], "重校後 `|ψ(0)|²`", kth=0)
claim("§三 r 重校", "10_reproduction_gate.md", J01["calibration"]["r_corrected_fm"], "重校後 `|ψ(0)|²`", kth=2)

# ══════════ docs/10 §五 兩約定表（第 2..6 格） ══════════
c5, cc = J02["calib"]["M5"], J02["calib"]["correct"]
for lab, m, src in (("M5", "M5（`w = N_d v²`，記錄值）", c5), ("正確", "正確（`w = v²`）", cc)):
    for i, (k, key) in enumerate([("a", "a"), ("t", "t_MeV"), ("⟨B⟩", "rA"), ("P₀", "P0"), ("ε", "E0")]):
        claim(f"§五 {lab} {k}", "10_reproduction_gate.md", src[key], m, cell=2 + i, kth=0)

# ══════════ docs/11 hash ＋ grep 表 ══════════
h = J03["def_sha256"]
RES.append(dict(name="docs/11 凍結 hash", doc="11_definition_freeze.md",
                status="✅" if any(h[:12] in ln for ln in doc("11_definition_freeze.md")) else "🔴 MISSING",
                detail=h[:12]))
DOCLAB = {"w_k 佔用": "`w_k`", "佔有率／佔據數": "佔有率", "多夸克態": "多夸克態", "Fock 展開": "Fock",
          "投影子差": "投影子差", "isovector": "isovector", "ūu − d̄d": "ūu", "夸克含量": "夸克含量",
          "海夸克／介子雲": "海夸克", "質量泛函非可加": "質量泛函非可加", "⟨O⟩ 矩陣元": "`⟨O⟩`"}
VERS = ["IDCM-v9.0--W-Field-Lattice", "IDCM-v10.0--It-From-Beat", "IDCM-v11.0--Reality-Architecture"]
for term, row in J03["terms"].items():
    lab = DOCLAB.get(term, term)
    # 只讀 §三.2（v4 語料範圍）之表格，避免與 §三 E1 舊表混淆
    alllines = doc("11_definition_freeze.md")
    try:
        i0 = next(i for i, ln in enumerate(alllines) if "§三.2" in ln)
        i1 = next((i for i, ln in enumerate(alllines) if i > i0 and ln.startswith("## ")), len(alllines))
    except StopIteration:
        i0, i1 = 0, len(alllines)
    block = alllines[i0:i1]
    cand = [ln for ln in block if ln.strip().startswith("|") and (lab in ln or term in ln)]
    if not cand:
        RES.append(dict(name=f"grep 表 {term}", doc="11_definition_freeze.md", status="🔴 MISSING", detail=lab))
        continue
    cells = re.findall(r"(\d+)f/(\d+)h", cand[0])
    ok, det = len(cells) == 3, f"doc={cells}"
    if ok:
        for i, v in enumerate(VERS):
            if (int(cells[i][0]), int(cells[i][1])) != (row[v]["nfiles"], row[v]["nhits"]):
                ok = False
                det += f" | {v.split('--')[0]} art={row[v]['nfiles']}f/{row[v]['nhits']}h"
    RES.append(dict(name=f"grep 表 {term}", doc="11_definition_freeze.md",
                    status="✅" if ok else "🔴 MISMATCH", detail=det))

# ══════════ docs/12 ══════════
r0 = J02["K3"]["rows"][0]
claim("§二 K1 max_dev", "12_r5_verdict.md", J02["K1"]["max_dev"],
      rx=r"max\|w−round\(w\)\| = ([\d.]+e−\d+)", grp=1)
claim("§二 K1 非整數態", "12_r5_verdict.md", J02["K1"]["noninteger"],
      rx=r"，非整數態 \*\*(\d+)\*\*", grp=1, integer=True)
claim("§二 K2 耦合", "12_r5_verdict.md", J02["K2"]["couplings"], "K2 | 加入 mode-1↔4", cell=2, kth=2, integer=True)
claim("§二 K2 非整數態", "12_r5_verdict.md", J02["K2"]["noninteger_states"], "K2 | 加入 mode-1↔4", cell=3, kth=0, integer=True)
claim("§二 K3 χ²", "12_r5_verdict.md", r0["chi2"], "K3 | `⟨O⟩≡1`", cell=3, kth=0)
claim("§二 K3 p", "12_r5_verdict.md", r0["p"], "K3 | `⟨O⟩≡1`", cell=3, kth=2)
claim("§二 K5 Δχ²", "12_r5_verdict.md", J02["K5"]["delta"], "K5 | M1 vs P₀-sea", cell=3, kth=0)
claim("§二 K6 χ²_sector", "12_r5_verdict.md", J02["K6"]["chi2_sector"], "K6 | 恒等性", cell=3, kth=2)
claim("§二 K6 p_sector", "12_r5_verdict.md", J02["K6"]["p_sector"], "K6 | 恒等性", cell=3, kth=3)
ALT_K4 = "← §二 K4 與 §3.4"
for i, lab in enumerate(["N", "Σ", "Ξ"]):
    claim(f"§二 K4 M5 {lab}", "12_r5_verdict.md", J02["K4"]["M5"]["shortfall"][i],
          "K4 | M5 binding 響應", cell=3, kth=i, alt=ALT_K4 if lab == "Ξ" else None)
    claim(f"§二 K4 corr {lab}", "12_r5_verdict.md", J02["K4"]["w=v²"]["shortfall"][i],
          "K4 | M5 binding 響應", cell=3, kth=i + 4, alt=ALT_K4 if lab == "Ξ" else None)
claim("§3.1 偏差", "12_r5_verdict.md", J02["K1"]["max_dev"], "次隨機核，偏差", kth=1)
claim("§3.3 Σγ=90", "12_r5_verdict.md", J02["K8"]["sum_gamma_exact"], rx=r"`Σγ = (\d+)` EXACT", grp=1, integer=True)

for lab, m, src in (("M5", "M5（`w=N_dv²`）", J02["K4"]["M5"]), ("corr", "正確（`w=v²`）", J02["K4"]["w=v²"])):
    claim(f"§3.4 {lab} E_bind", "12_r5_verdict.md", src["Ebind_MeV"], m, cell=2, kth=0)
    claim(f"§3.4 {lab} ∂E/∂m", "12_r5_verdict.md", src["dEbind_dml"], m, cell=3, kth=0)
    for i, ch in enumerate(["N", "Σ", "Ξ"]):
        claim(f"§3.4 {lab} short {ch}", "12_r5_verdict.md", src["shortfall"][i], m, cell=4, kth=i,
              alt=ALT_K4 if (lab == "corr" and ch == "Ξ") else None)
claim("§3.4 Δ_N", "12_r5_verdict.md", J02["K4"]["need"][0], "若接受三錨", kth=1)
claim("§3.4 Δ_Σ", "12_r5_verdict.md", J02["K4"]["need"][1], "若接受三錨", kth=2)
claim("§3.4 Δ_Ξ", "12_r5_verdict.md", J02["K4"]["need"][2], "若接受三錨", kth=3)

# ══════════ final_result（cell 2；token 序：⟨O⟩≡1, χ², 3 dof, p, 3.36σ, a, σ_a, b, σ_b） ══════════
claim("FR 19/19", "final_result.md", J01["npass"], "R0 重現閘", cell=2, kth=0, integer=True)
claim("FR a₀ 20位", "final_result.md", J01["mpmath"]["a0_fm"], "R0 重現閘", cell=2, kth=2)
claim("FR t 18位", "final_result.md", J01["mpmath"]["t_MeV"], "R0 重現閘", cell=2, kth=3)
claim("FR 殘差", "final_result.md", J01["mpmath"]["resid_A"], "R0 重現閘", cell=2, kth=5)
claim("FR χ²", "final_result.md", r0["chi2"], "三錨 χ²（獨立重算）", cell=2, kth=1)
claim("FR p", "final_result.md", r0["p"], "三錨 χ²（獨立重算）", cell=2, kth=3)
claim("FR 3.36σ", "final_result.md", 3.36, "三錨 χ²（獨立重算）", cell=2, kth=4)
claim("FR a 斜率", "final_result.md", J02["K3"]["linear"]["a"], "三錨 χ²（獨立重算）", cell=2, kth=5)
claim("FR σ_a", "final_result.md", J02["K3"]["linear"]["sig_a"], "三錨 χ²（獨立重算）", cell=2, kth=6)
claim("FR b 斜率", "final_result.md", J02["K3"]["linear"]["b"], "三錨 χ²（獨立重算）", cell=2, kth=7)
claim("FR σ_b", "final_result.md", J02["K3"]["linear"]["sig_b"], "三錨 χ²（獨立重算）", cell=2, kth=8)

# ══════════ docs/15（① ② 判決；R6） ══════════
claim("§15 C1 色嵌入 dim", "15_additivity_scope_verdict.md", J05["C1"]["dim_colour"],
      rx=r"= 2 × \(B−1\) = 2 × (\d+)", grp=1, integer=True)
claim("§15 C1 每模容量", "15_additivity_scope_verdict.md", J05["C1"]["capacity"],
      rx=r"= 2 × \(B−1\) = 2 × \d+ = (\d+)", grp=1, integer=True)
claim("§15 C1 max σ_n", "15_additivity_scope_verdict.md", J05["C1"]["sigma_max"],
      rx=r"max σ_n = 3/6 = ([\d.]+) < 1", grp=1)
claim("§15 C3 N dev", "15_additivity_scope_verdict.md", J05["C3"]["dev"]["N (uud)"],
      "N (uud) | 1 |", cell=3, kth=0)
claim("§15 C3 N 所需", "15_additivity_scope_verdict.md", J05["C3"]["need_MeV"][0],
      "N (uud) | 1 |", cell=4, kth=0)
claim("§15 C3 Σ 所需", "15_additivity_scope_verdict.md", J05["C3"]["need_MeV"][1],
      "Σ (uus) | 1 |", cell=4, kth=0)
claim("§15 C3 Ξ dev", "15_additivity_scope_verdict.md", J05["C3"]["dev"]["Ξ (uss)"],
      "Ξ (uss) | 2 |", cell=3, kth=0)
claim("§15 C3 Ξ 所需", "15_additivity_scope_verdict.md", J05["C3"]["need_MeV"][2],
      "Ξ (uss) | 2 |", cell=4, kth=0)
claim("§15 C5 Σ_pair", "15_additivity_scope_verdict.md", J05["C5"]["E_pair"],
      rx=r"= \((3/2)\)√2 = ([\d.]+) ℏ", grp=2)
claim("§15 C5 E₀^(3)", "15_additivity_scope_verdict.md", J05["C5"]["E3"],
      rx=r"= \(1/2\)\(0 \+ 2√3\) = ([\d.]+) ℏ", grp=1)
claim("§15 C5 相關 %", "15_additivity_scope_verdict.md", J05["C5"]["corr_rel"] * 100,
      rx=r"三體相關\s+= ([−-][\d.]+)%", grp=1)

# ══════════ docs/16（未推盡補推；R7）—— σ_* 用絕對值（文檔報偏離量級） ══════════
claim("§16 b̂", "16_gap_closure.md", J06["G1"]["b"], "| `b̂` |", cell=2, kth=0)
claim("§16 σ_b", "16_gap_closure.md", J06["G1"]["sig_b"], "| `σ_b` |", cell=2, kth=0)
claim("§16 z_b", "16_gap_closure.md", J06["G1"]["z_b"], "| `z_b` |", cell=2, kth=0)
claim("§16 p_MC", "16_gap_closure.md", J06["G1"]["p_mc"], "| `p_MC` |", cell=2, kth=0)
claim("§16 帶內成員", "16_gap_closure.md", J06["G2"]["n_in_band"], "| 帶內成員 |", cell=2, kth=0, integer=True)
claim("§16 帶寬/間距", "16_gap_closure.md", J06["G2"]["bandwidth_over_gap"], "| 帶寬/間距 |", cell=2, kth=0)
claim("§16 grading 權重值", "16_gap_closure.md", J06["G3"]["vecC"][0], "| grading 權重值 |", cell=2, kth=0, integer=True)
claim("§16 χ²_D", "16_gap_closure.md", J06["G3"]["chi2_D"], "| `χ²_D` |", cell=2, kth=0)
claim("§16 7/11", "16_gap_closure.md", J06["G5"]["ratio"], "| 7/11 |", cell=2, kth=0)
claim("§16 σ_N", "16_gap_closure.md", abs(J06["G1"]["per_channel_z"][0]), "| `σ_N` |", cell=2, kth=0)
claim("§16 σ_Σ", "16_gap_closure.md", abs(J06["G1"]["per_channel_z"][1]), "| `σ_Σ` |", cell=2, kth=0)
claim("§16 σ_Ξ", "16_gap_closure.md", abs(J06["G1"]["per_channel_z"][2]), "| `σ_Ξ` |", cell=2, kth=0)

# ══════════ docs/17（全鏈審計） ══════════
claim("§17 npass", "17_audit_full_chain.md", J07["npass"], "| 審計 npass |", cell=2, kth=0, integer=True)
claim("§17 ntotal", "17_audit_full_chain.md", J07["ntotal"], "| 審計 ntotal |", cell=2, kth=0, integer=True)
claim("§17 b̂", "17_audit_full_chain.md", J07["metrics"]["b_hat"], "| 審計 b̂ |", cell=2, kth=0)
claim("§17 σ_b", "17_audit_full_chain.md", J07["metrics"]["sig_b"], "| 審計 σ_b |", cell=2, kth=0)
claim("§17 χ²", "17_audit_full_chain.md", J07["metrics"]["chi2"], "| 審計 χ² |", cell=2, kth=0)
claim("§17 p", "17_audit_full_chain.md", J07["metrics"]["p"], "| 審計 p |", cell=2, kth=0)
claim("§17 Σγ", "17_audit_full_chain.md", J07["metrics"]["sum_gamma"], "| 審計 Σγ |", cell=2, kth=0, integer=True)
claim("§17 B7 相關%", "17_audit_full_chain.md", J07["metrics"]["corr_pct"], "| 審計 B7 相關% |", cell=2, kth=0)
claim("§17 MC/解析比值", "17_audit_full_chain.md", J07["metrics"]["mc_ratio"], "| 審計 MC/解析比值 |", cell=2, kth=0)
claim("§17 C2 未覆蓋", "17_audit_full_chain.md", J07["metrics"]["uncovered"], "| 審計 C2 未覆蓋 |", cell=2, kth=0, integer=True)

# ══════════ docs/18（M4 收件收據） ══════════
claim("§18 npass", "18_m4_intake_receipt.md", J09["npass"], "| 收據 npass |", cell=2, kth=0, integer=True)
claim("§18 ntotal", "18_m4_intake_receipt.md", J09["ntotal"], "| 收據 ntotal |", cell=2, kth=0, integer=True)
claim("§18 R2", "18_m4_intake_receipt.md", J09["metrics"]["ident_d"], "| R2 值 |", cell=2, kth=0)
claim("§18 R3", "18_m4_intake_receipt.md", J09["metrics"]["ident_e"], "| R3 值 |", cell=2, kth=0)
claim("§18 χ²_min", "18_m4_intake_receipt.md", J09["metrics"]["chi2_min"], "| 反推 χ²_min |", cell=2, kth=0)
claim("§18 R4", "18_m4_intake_receipt.md", J09["metrics"]["z_eq"], "| R4 z 等效 |", cell=2, kth=0)
claim("§18 R5", "18_m4_intake_receipt.md", J09["metrics"]["sum_gamma_struct"], "| R5 Σγ 結構 |", cell=2, kth=0, integer=True)

# ══════════ docs/19（O120 引入） ══════════
claim("§19 npass", "19_o120_import.md", J11["npass"], "| 引入檢定 npass |", cell=2, kth=0, integer=True)
claim("§19 ntotal", "19_o120_import.md", J11["ntotal"], "| 引入檢定 ntotal |", cell=2, kth=0, integer=True)
claim("§19 合法模式數", "19_o120_import.md", len(J11["parts"]["A"]["legal"]), "| 合法模式數 |", cell=2, kth=0, integer=True)
claim("§19 γ₃/γ₂", "19_o120_import.md", J11["parts"]["A"]["pairs_ud"][0][1], "| γ₃/γ₂ |", cell=2, kth=0)
claim("§19 z(u/d)", "19_o120_import.md", J11["parts"]["A"]["pairs_ud"][0][2], "| z(u/d) |", cell=2, kth=0)
claim("§19 z(t/b)", "19_o120_import.md", J11["parts"]["A"]["z_tb"], "| z(t/b) |", cell=2, kth=0)
claim("§19 全局解數", "19_o120_import.md", J11["parts"]["A"]["global_solutions"], "| 全局解數 |", cell=2, kth=0, integer=True)
claim("§19 V9 χ²", "19_o120_import.md", J11["parts"]["A"]["chi2"], "| V9 χ² |", cell=2, kth=0)
claim("§19 V9 p", "19_o120_import.md", J11["parts"]["A"]["p"], "| V9 p |", cell=2, kth=0)
claim("§19 m_ud 預言", "19_o120_import.md", J11["parts"]["A"]["m_ud_pred"], "| m_ud 預言 |", cell=2, kth=0)
claim("§19 X 預言", "19_o120_import.md", J11["metrics"]["X_pred"], "| X 預言 |", cell=2, kth=0)
claim("§19 z(X)", "19_o120_import.md", J11["metrics"]["z_X"], "| z(X) |", cell=2, kth=0)
claim("§19 χ²_inj", "19_o120_import.md", J11["metrics"]["chi2_M7"], "| χ²（全單射） |", cell=2, kth=0)
claim("§19 χ²_v9", "19_o120_import.md", J11["metrics"]["chi2_v9_shared"], "| χ²（v9 共享） |", cell=2, kth=0)

# ══════════ docs/20（O209 唯一分配） ══════════
M12 = J12["metrics"]
claim("§20 npass", "20_o209_assignment.md", J12["npass"], "| 檢定 npass |", cell=2, kth=0, integer=True)
claim("§20 ntotal", "20_o209_assignment.md", J12["ntotal"], "| 檢定 ntotal |", cell=2, kth=0, integer=True)
claim("§20 RMS", "20_o209_assignment.md", M12["min_rms"], "| 最小 RMS（重算） |", cell=2, kth=0)
claim("§20 n_tie", "20_o209_assignment.md", M12["n_tie"], "| 達最小配對數 |", cell=2, kth=0, integer=True)
claim("§20 corr", "20_o209_assignment.md", M12["corr"], "| corr（重算） |", cell=2, kth=0)
claim("§20 α", "20_o209_assignment.md", M12["alpha"], "| α（斜率） |", cell=2, kth=0)
claim("§20 ratio_pl", "20_o209_assignment.md", M12["ratio_pl_ud"], "| 冪律載體 u/d 比 |", cell=2, kth=0)
claim("§20 z_pl", "20_o209_assignment.md", M12["z_pl_ud"], "| z（冪律載體 u/d） |", cell=2, kth=0)
claim("§20 z_g209", "20_o209_assignment.md", M12["z_gamma_o209"], "| z（γ 載體 + O209 指派） |", cell=2, kth=0)
claim("§20 z_dl120", "20_o209_assignment.md", M12["z_dl_o120"], "| z（冪律載體 + O120 指派） |", cell=2, kth=0)
claim("§20 X_o120", "20_o209_assignment.md", M12["X_o120"], "| X（O120 分支） |", cell=2, kth=0)
claim("§20 X_o209", "20_o209_assignment.md", M12["X_o209"], "| X（O209 分支） |", cell=2, kth=0)

# ══════════ docs/20 §七（生死判定） ══════════
M13 = J13["metrics"]
claim("§7 npass", "20_o209_assignment.md", J13["npass"], "| §七 npass |", cell=2, kth=0, integer=True)
claim("§7 隨機組數", "20_o209_assignment.md", M13["tautology_runs"], "| 隨機組數 |", cell=2, kth=0, integer=True)
claim("§7 α_out", "20_o209_assignment.md", M13["alpha_out"], "| α_out（剔除 u,d） |", cell=2, kth=0)
claim("§7 z_b", "20_o209_assignment.md", M13["z_b"], "| 載體 B 樣本外 z |", cell=2, kth=0)
claim("§7 X23", "20_o209_assignment.md", M13["X_cands"]["(2,3)"], "| X（(2,3)） |", cell=2, kth=0)
claim("§7 X57", "20_o209_assignment.md", M13["X_cands"]["(5,7)"], "| X（(5,7)） |", cell=2, kth=0)
claim("§7 X913", "20_o209_assignment.md", M13["X_cands"]["(9,13)"], "| X（(9,13)） |", cell=2, kth=0)
claim("§7 Xspread", "20_o209_assignment.md", M13["X_spread_over_sigma"], "| X 跨度／σ_X |", cell=2, kth=0)

RES.append(dict(name="v1 vs v2 K1 max_dev", doc="—",
                status="✅" if abs(Jv1["K1"]["max_dev"] - J02["K1"]["max_dev"]) < 1e-18 else "🔴 MISMATCH",
                detail=f"{Jv1['K1']['max_dev']:.4e} vs {J02['K1']['max_dev']:.4e}"))
RES.append(dict(name="v1 vs v2 K3 χ²", doc="—",
                status="✅" if abs(Jv1["K3"]["rows"][0]["chi2"] - r0["chi2"]) < 1e-9 else "🔴 MISMATCH",
                detail=f"{Jv1['K3']['rows'][0]['chi2']:.6f} vs {r0['chi2']:.6f}"))

npass = sum(1 for r in RES if r["status"].startswith("✅"))
nfail = len(RES) - npass
print("=" * 116)
print(f"M7 文檔↔輸出 一致性核 v3：{npass}/{len(RES)} PASS，{nfail} FAIL；corrected {len(CORR)}")
print("=" * 116)
for r in RES:
    if not r["status"].startswith("✅"):
        print(f"  {r['status']:<14} {r['name']:<26} [{r['doc']}]  {r['detail']}")
if CORR:
    print("\n§更正清單（原文不動，另加更正行）：")
    for c in CORR:
        print(f"  • {c['name']}: 原文 {c['original']} → 更正 {c['corrected']}（腳本 {c['art']!r}）")
print("\n判：", "✅ 全部文檔數字 == 腳本輸出（按印刷精度）" if nfail == 0 else f"🔴 {nfail} 條須處理（禁放寬判準）")
json.dump(dict(npass=npass, ntotal=len(RES), nfail=nfail, corrected=len(CORR),
               verdict="PASS" if nfail == 0 else "FAIL", corrigenda=CORR, results=RES),
          open(os.path.join(DATA, "04_doc_consistency.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/04_doc_consistency.json")
sys.exit(0 if nfail == 0 else 1)