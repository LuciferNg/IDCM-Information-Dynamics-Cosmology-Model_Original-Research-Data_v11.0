"""M7 scripts/07 — 全鏈審計（idcm-audit-verify §2/§3 機械化）
================================================================================================
觸發：調律者 2026-10-26「先進行 audit-verify M7 session 全鏈」。
本腳本**只讀**（除自身 JSON/log），分三層：

A 結構層（§3.1 六項）
  A1 檔案對稱：README 所列 ↔ 磁碟實檔（兩側差集須空）
  A2 產物綁定：每個 data/*.json 有 producer；每個 script 宣告之輸檔存在（無孤兒）
  A3 路徑衛生：scripts/ 內 /tmp、scratch 寫入 = 0
  A4 driver：run_chain.sh 存在與否
  A5 版本檔：data/*.v1.* 是否被 doc 引用（孤兒版本檔偵測）

B 獨立重推導（L2：**不同算法**重算 headline 數字）
  B1 γ_n 兩法：|1−ζ^n|² vs 4sin²（50 dps）
  B2 Σγ = 90（n=0..44；22 對×4 + γ_22=2）
  B3 χ² = 18.221 與 p（mpmath 上不完全 Γ，獨立於 02b 之 p 實作）
  B4 b̂/σ_b：Fraction 精確正規方程（獨立於 numpy lstsq）
  B5 MC p（2.79e−5）vs 解析常態尾（跨方法互核）
  B6 G3 導數：sympy 對 M=Σn_q m_q **符號微分**（非手推公式）
  B7 諧振 −18.35%：sympy eigenvals（非手寫 √3）
  B8 容量 6 > N_q=3；B9 γ₄M₀ 與 docs/15 C3 兩值
C 記錄層
  C1 四件套狀態行一致 ＋ register state
  C2 水位：close 之後被改嘅檔清單
  C3 append-only：git numstat 掃歷史檔（01_problem_record／docs/*／SESSION_START）刪除行
判準全部寫於執行前。輸出：data/07_audit_full_chain.json
"""
import json
import os
import re
import subprocess
from fractions import Fraction

from mpmath import mp, mpf, gammainc, gamma, pi, sin, sqrt, cos, exp, inf
import sympy as sp

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
V11 = os.path.abspath(os.path.join(M7, "..", ".."))
REPO = os.path.abspath(os.path.join(V11, ".."))
REG = os.path.join(REPO, "problem_register", "register.jsonl")

OUT = {"module": "M7", "script": "07_audit_full_chain.py", "checks": [], "findings": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<50} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))
    if not cond:
        OUT["findings"].append(dict(item=name, got=str(got), want=str(want), note=note))


def load(n):
    return json.load(open(os.path.join(M7, "data", n)))


J01, J02b, J03 = load("01_reproduction_gate.json"), load("02b_r5_mode_content_nogo.json"), load("03_grep_evidence.json")
J04, J05, J06 = load("04_doc_consistency.json"), load("05_additivity_scope.json"), load("06_gap_closure.json")
AN = J01["anchors"]
Y = [Fraction(str(AN[k][0])) for k in ("N", "Sigma", "Xi")]
S = [Fraction(str(AN[k][1])) for k in ("N", "Sigma", "Xi")]
M0 = Fraction("27.9009145112")
MDMU = Fraction(str(AN["md_mu"][0]))

print("=" * 116)
print("M7 全鏈審計（scripts/07）—— A 結構 / B 獨立重推導 / C 記錄層")
print("=" * 116)

# ══════════ A1 檔案對稱 ══════════
print("\nA1 檔案對稱（README 所列 ↔ 磁碟）")
rd = open(os.path.join(M7, "README.md"), encoding="utf-8").read()
toks = set(re.findall(r"`((?:scripts|docs|data)/[\w.\-]+)`", rd))
disk = set()
for sub in ("scripts", "docs", "data"):
    for f in os.listdir(os.path.join(M7, sub)):
        if os.path.isfile(os.path.join(M7, sub, f)):
            disk.add(f"{sub}/{f}")
siblings = [d for d in os.listdir(V11) if os.path.isdir(os.path.join(V11, d))] + \
           [os.path.join("main_problems", x) for x in os.listdir(os.path.join(V11, "main_problems"))]
missing_on_disk, cross_m = [], []
docs_disk = os.listdir(os.path.join(M7, "docs"))
for t in sorted(toks):
    if t in disk:
        continue
    if t.startswith("docs/") and "." not in os.path.basename(t):
        base = os.path.basename(t)
        if any(f.startswith(base) for f in docs_disk):
            continue   # 本 M docs/NN 簡寫
        if any(os.path.isdir(os.path.join(V11, sib, "docs")) and
               any(f.startswith(base) for f in os.listdir(os.path.join(V11, sib, "docs")))
               for sib in siblings):
            cross_m.append(t)
            continue   # 跨 M docs/NN 簡寫（如 M4 之 `docs/28`）
    if any(os.path.exists(os.path.join(V11, sib, t)) for sib in siblings):
        cross_m.append(t)
    else:
        missing_on_disk.append(t)
unlisted = sorted(d for d in disk if d not in toks and not d.startswith("data/")
                  and (d.startswith("scripts/") or d.startswith("docs/")))
print(f"    README 提 {len(toks)} 項；磁碟 {len(disk)} 檔（含 log/版本檔）")
print(f"    跨 M／O 引用（存在於他目錄，合法）：{cross_m or '無'}")
print(f"    README 有而全樹各處都冇：{missing_on_disk or '無'}")
print(f"    磁碟有而 README 冇（除 run log／.v1）：{unlisted or '無'}")
ok("A1 README 引用全部存在", not missing_on_disk, missing_on_disk or "0", "0")
ok("A1 所有 scripts/*.py 與 docs/*.md 已列於 README", not unlisted, unlisted or "0", "0")

# ══════════ A2 產物綁定 ══════════
print("\nA2 產物綁定（data/*.json ↔ producer script）")
scripts = {f: open(os.path.join(M7, "scripts", f), encoding="utf-8").read()
           for f in os.listdir(os.path.join(M7, "scripts")) if f.endswith(".py")}
jsons = [f for f in os.listdir(os.path.join(M7, "data")) if f.endswith(".json")]
HIST = {"02_r5_mode_content_nogo.v1.json": "02 之 v1 歷史快照（原 run 改名保存；docs/17 已記）"}
from collections import defaultdict
CROSS = defaultdict(list)
for root, dirs, fs in os.walk(REPO):
    if root.count(os.sep) - REPO.count(os.sep) > 4:
        dirs[:] = []
        continue
    if os.path.basename(root) == "data":
        for f in fs:
            if f.endswith(".json"):
                CROSS[f].append(root)
orphan_json, orphan_out, hist_unref, cross_out = [], [], [], []
for j in jsons:
    if j in HIST:
        hist_unref.append(j)
        continue
    if not any(j in src for src in scripts.values()):
        orphan_json.append(j)
for f, src in scripts.items():
    # 每個 script 宣告嘅輸檔（docstring 內 "data/<name>.json"）必須存在
    for m in re.findall(r"data/([\w.\-]+\.json)", src):
        if m.startswith("07_") and f.startswith("07_"):
            continue   # 自身輸出（執行尾才寫）
        if not os.path.exists(os.path.join(M7, "data", m)):
            if m in CROSS:            # 跨 O／M 之 data 輸入（合法，需標明來源）
                cross_out.append(f"{f} -> {CROSS[m][0]}/{m}")
            else:
                orphan_out.append(f"{f} -> data/{m}")
print(f"    JSON {len(jsons)} 個；無 producer 者：{orphan_json or '無'}；歷史快照（白名單）：{hist_unref or '無'}")
print(f"    跨 O data 輸入（合法）：{len(cross_out)} 項：{cross_out[:3]}")
print(f"    script 宣告但缺檔：{orphan_out or '無'}")
ok("A2 無孤兒 JSON", not orphan_json, orphan_json or "0", "0")
ok("A2 無宣告缺檔", not orphan_out, orphan_out or "0", "0")

# ══════════ A3 路徑衛生 ══════════
print("\nA3 路徑衛生（scripts/ 內 /tmp、scratch 寫入）")
bad = {}
for f, src in scripts.items():
    if f.startswith("07_"):
        continue   # 自身含 pattern 字面，排除
    hits = [ln.strip() for ln in src.splitlines()
            if re.search(r"['\"]/tmp/|scratch", ln) and not ln.strip().startswith("#")]
    if hits:
        bad[f] = hits[:3]
print(f"    命中檔：{list(bad) or '無'}")
ok("A3 scripts 無 /tmp／scratch 寫入", not bad, len(bad), 0)

# ══════════ A4 driver ══════════
print("\nA4 driver")
drv = os.path.exists(os.path.join(M7, "scripts", "run_chain.sh"))
print(f"    scripts/run_chain.sh：{'存在' if drv else '**不存在**（§3.1 ⑥ gap）'}")
OUT["A4"] = dict(run_chain_sh=drv)
if not drv:
    OUT["findings"].append(dict(item="A4 driver", got="missing", want="run_chain.sh", note="§3.1 ⑥"))

# ══════════ A5 版本檔 ══════════
print("\nA5 版本檔（.v1）是否被引用")
v1 = [f for f in jsons if ".v1" in f]
_all = [os.path.join(M7, "docs", x) for x in os.listdir(os.path.join(M7, "docs"))] + [os.path.join(M7, "README.md")]
refd = {f: any(f in open(d, encoding="utf-8").read() for d in _all) for f in v1}
print(f"    {v1} ⟹ 被 doc 引用: {refd}")
ok("A5 版本檔有記錄歸屬（被 doc 提及）", all(refd.values()) if v1 else True, refd, "全部 True",
   "v1＝02 之舊版；docs/12 命名為 v1/v2 對照")

# ══════════ B1 γ_n 兩法 ══════════
print("\nB1 γ_n 兩法（|1−ζ^n|² vs 4sin²；50 dps）")
devs = []
for n in range(0, 45):
    z = exp(2j * pi * n / 88) if False else mp.mpc(cos(2 * pi * n / 88), sin(2 * pi * n / 88))
    a = abs(1 - z) ** 2
    b = 4 * sin(pi * n / 88) ** 2
    devs.append(abs(a - b))
print(f"    n=0..44；max dev = {mp.nstr(max(devs), 5)}")
ok("B1 兩法一致（<1e−45）", max(devs) < mpf("1e-45"), mp.nstr(max(devs), 3), "<1e−45")
# 配對恆等 γ_n + γ_{44−n} = 4
pd = [abs(4 * sin(pi * n / 88) ** 2 + 4 * sin(pi * (44 - n) / 88) ** 2 - 4) for n in range(0, 45)]
ok("B1 配對恆等 γ_n+γ_{44−n}=4（n=0..44）", max(pd) < mpf("1e-45"), mp.nstr(max(pd), 3), "<1e−45")

# ══════════ B2 Σγ = 90 ══════════
print("\nB2 Σγ（n=0..44）")
s45 = sum(4 * sin(pi * n / 88) ** 2 for n in range(0, 45))
s22 = 22 * 4 + 2  # 22 對×4 + γ_22=2（結構計數，另一法）
print(f"    數值 Σ_{{n=0..44}} γ_n = {mp.nstr(s45, 12)}；結構計數 22×4+2 = {s22}")
ok("B2 Σγ = 90（兩法）", abs(s45 - 90) < mpf("1e-40") and s22 == 90, mp.nstr(s45, 8), 90)

# ══════════ B3 χ² 與 p（獨立於 02b 實作） ══════════
print("\nB3 χ²(⟨O⟩≡1) 與 p（mpmath 上不完全 Γ）")
chi2 = sum(float((Fraction(1) - y) ** 2 / s ** 2) for y, s in zip(Y, S))
p3 = gammainc(mpf(3) / 2, mpf(chi2) / 2, inf) / gamma(mpf(3) / 2)
print(f"    χ² = {chi2:.6f}；p(3 dof) = {mp.nstr(p3, 8)}")
print(f"    02b K3 報：χ² = {J02b['K3']['rows'][0]['chi2']:.6f}，p = {J02b['K3']['rows'][0]['p']:.4e}")
ok("B3 χ² 重算一致（<1e−9）", abs(chi2 - J02b["K3"]["rows"][0]["chi2"]) < 1e-9, f"{chi2:.9f}",
   f"{J02b['K3']['rows'][0]['chi2']:.9f}")
ok("B3 p 重算一致（相對 <2%）", abs(float(p3) - J02b["K3"]["rows"][0]["p"]) / float(p3) < 0.02,
   mp.nstr(p3, 6), f"{J02b['K3']['rows'][0]['p']:.4e}")

# ══════════ B4 b̂ / σ_b（Fraction 精確正規方程） ══════════
print("\nB4 b̂/σ_b（Fraction 精確；獨立於 numpy lstsq）")
w = [1 / (s * s) for s in S]
xs = [Fraction(0), Fraction(1), Fraction(2)]
nbar = sum(wi * xi for wi, xi in zip(w, xs)) / sum(w)
Sxx = sum(wi * (xi - nbar) ** 2 for wi, xi in zip(w, xs))
bhat = sum(wi * (xi - nbar) * yi for wi, xi, yi in zip(w, xs, Y)) / Sxx
ahat = sum(wi * yi for wi, yi in zip(w, Y)) / sum(w) - bhat * nbar
sig_b = 1 / sqrt(mpf(Sxx.numerator) / mpf(Sxx.denominator))
sig_a_num = mpf(1) / sqrt(mpf(sum(w).numerator) / mpf(sum(w).denominator) + (mpf(nbar.numerator) / nbar.denominator) ** 2 / (mpf(Sxx.numerator) / Sxx.denominator))
print(f"    b̂ = {float(bhat):.9f}  σ_b = {mp.nstr(sig_b, 8)}  â = {float(ahat):.9f}  σ_a = {mp.nstr(sig_a_num, 6)}")
print(f"    06 G1 報：b̂ = {J06['G1']['b']:.9f}  σ_b = {J06['G1']['sig_b']:.9f}  â = {J06['G1']['a']:.9f}")
ok("B4 b̂ 一致（<1e−9）", abs(float(bhat) - J06["G1"]["b"]) < 1e-9, f"{float(bhat):.9f}", f"{J06['G1']['b']:.9f}")
ok("B4 σ_b 一致（<1e−8）", abs(float(sig_b) - J06["G1"]["sig_b"]) < 1e-8, mp.nstr(sig_b, 8), f"{J06['G1']['sig_b']:.9f}")

# ══════════ B5 MC p vs 解析常態尾 ══════════
print("\nB5 MC p（2.79e−5）vs 解析常態尾（跨方法）")
z = abs(float(bhat)) / float(sig_b)
p_an = mp.gammainc(mpf(1) / 2, mpf(z) ** 2 / 2, inf) / gamma(mpf(1) / 2)  # = erfc(z/√2) = 2Φ(−z)
print(f"    z = {z:.6f}；解析雙尾 p = Q(1/2,z²/2) = {mp.nstr(p_an, 6)}；MC(06) = {J06['G1']['p_mc']:.4e}")
ratio = float(p_an) / J06["G1"]["p_mc"]
ok("B5 MC 與解析一致（比值 0.9–1.1）", 0.9 < ratio < 1.1, f"比值 {ratio:.4f}", "≈1.0",
   "MC 誤差 ~1/√(10⁷p)≈6% ⟹ 容差 10%")

# ══════════ B6 G3 導數：sympy 符號微分 ══════════
print("\nB6 G3 導數（sympy 對 M=Σn_q m_q 符號微分）")
mu, md_, ms = sp.symbols("m_u m_d m_s", positive=True)
PAIRS = {"N": ((2, 1, 0), (1, 2, 0), 1), "Σ": ((2, 0, 1), (0, 2, 1), 2), "Ξ": ((1, 0, 2), (0, 1, 2), 1)}
vals = {}
for h, (hi, lo, coeff) in PAIRS.items():
    Mhi = sum(n * m for n, m in zip(hi, (mu, md_, ms)))
    Mlo = sum(n * m for n, m in zip(lo, (mu, md_, ms)))
    d = sp.diff(Mhi - Mlo, mu)  # ∂/∂m_u（M4 符號取正）
    v = sp.nsimplify(d / coeff)
    vals[h] = int(v)
print(f"    ∂(ΔM)/∂m_u / coeff = {vals}")
ok("B6 三通道皆 1（與 06 讀法 A/B 一致）", list(vals.values()) == [1, 1, 1], vals, [1, 1, 1])

# ══════════ B7 諧振 −18.35%（sympy eigenvals） ══════════
print("\nB7 pairwise 三體相關（sympy eigenvals，非手寫 √3）")
K = sp.Matrix([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]])  # 等質量 pairwise 彈簧 k=m=1
ev = sorted(K.eigenvals().items())
omegas = sorted([sp.sqrt(sp.Rational(e)) for e, m in K.eigenvals().items() for _ in range(m)])
E3 = sum(o for o in omegas) / 2
Epair = sp.Rational(3, 2) * sp.sqrt(2)
corr = (E3 - Epair) / Epair
print(f"    K 本徵值 {[e for e, _ in ev]} ⟹ E₀^(3) = {sp.nsimplify(E3)} = {float(E3):.6f}")
print(f"    Σ_pair = {float(Epair):.6f}（3×½√2）；相關 = {float(corr)*100:+.4f}%")
ok("B7 與 05 C5 記錄一致（−18.35%）", abs(float(corr) - J05["C5"]["corr_rel"]) < 1e-9,
   f"{float(corr)*100:+.6f}%", f"{J05['C5']['corr_rel']*100:+.2f}%")

# ══════════ B8/B9 容量與 γ₄M₀ ══════════
print("\nB8/B9 容量與 docs/15 C3 值")
cap = 2 * (4 - 1)
g4 = 4 * sin(pi * 4 / 88) ** 2
dev_xi = -g4 * mpf(M0.numerator) / M0.denominator
need_xi = mpf(1293) / 1000 * mpf(MDMU.numerator) / MDMU.denominator
print(f"    cap = 2×(B−1) = {cap} > N_q = 3；γ₄M₀ = {mp.nstr(dev_xi, 8)} MeV；需 +{mp.nstr(need_xi, 6)} MeV")
ok("B8 容量 6 > 3", cap == 6, cap, 6)
ok("B9 γ₄M₀ 與 05 C3 一致（−2.2605）", abs(float(dev_xi) - J05["C3"]["dev"]["Ξ (uss)"]) < 1e-9,
   mp.nstr(dev_xi, 8), f"{J05['C3']['dev']['Ξ (uss)']:.6f}")
ok("B9 所需 = Δ×2.56 一致（+3.310）", abs(float(need_xi) - J05["C3"]["need_MeV"][2]) < 1e-3,
   mp.nstr(need_xi, 6), f"{J05['C3']['need_MeV'][2]:.3f}")

# ══════════ C1 四件套 ＋ register state ══════════
print("\nC1 四件套狀態行 ＋ register state")
def first_status(p, pat=r"狀態"):
    for ln in open(p, encoding="utf-8"):
        if pat in ln:
            return ln.strip()[:120]
    return "(無)"
four = {f: first_status(os.path.join(M7, f)) for f in
        ("01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md")}
for f, s in four.items():
    print(f"    {f:<26} {s}")
regstate = None
for ln in open(REG, encoding="utf-8"):
    try:
        d = json.loads(ln)
    except Exception:
        continue
    if d.get("pid") == "M7" and d.get("state"):
        regstate = (d["state"], d.get("verdict"), d.get("ts"))
regstate = regstate or (None, None, None)
print(f"    register: state={regstate[0]} verdict={regstate[1]} ts={regstate[2]}")
consist = regstate[0] == "closed" and regstate[1] == "negative"
ok("C1 register = closed/negative", consist, f"{regstate[0]}/{regstate[1]}", "closed/negative")
blocks = {}
for f in ("01_problem_record.md", "SESSION_START.md"):
    txt = open(os.path.join(M7, f), encoding="utf-8").read()
    m = re.search(r"state:\s*(\w+)", txt)
    v = re.search(r"verdict:\s*(\w+)", txt)
    blocks[f] = (m.group(1) if m else None, v.group(1) if v else None)
print(f"    registry block: {blocks}")
ok("C1 registry block == register state（無 BLOCK-DRIFT）",
   all(b[0] == "closed" and b[1] == "negative" for b in blocks.values()), blocks, "closed/negative")

# ══════════ C2 水位（close 之後被改嘅檔） ══════════
print("\nC2 水位：register close（2026-10-26T12:59:45）之後被改嘅檔")
close_ts = 0.0
import datetime as dt
close_ts = dt.datetime.fromisoformat("2026-10-26T12:59:45+08:00").timestamp()
changed = []
for sub in ("docs", "scripts", "data"):
    for f in sorted(os.listdir(os.path.join(M7, sub))):
        p = os.path.join(M7, sub, f)
        if os.path.isfile(p) and os.path.getmtime(p) > close_ts:
            changed.append(f"{sub}/{f}")
for f in ("01_problem_record.md", "SESSION_START.md", "README.md", "final_result.md"):  # 四件套
    p = os.path.join(M7, f)
    if os.path.isfile(p) and os.path.getmtime(p) > close_ts:
        changed.append(f)
print(f"    {len(changed)} 檔：{changed}")
WM = json.load(open(os.path.join(M7, "data", "watermark.json")))
wm_paths = {r["path"] for r in WM["rows"]}
uncovered = [c for c in changed if c not in wm_paths]
print(f"    水位（data/watermark.json，scripts/10 生成）列 {len(wm_paths)} 檔；未覆蓋 {len(uncovered)}：{uncovered[:8]}")
ok("C2 水位覆蓋全部 close 後變更檔", not uncovered, len(uncovered), 0,
   "水位屬 live field，由 scripts/10 重生；本檢查只核覆蓋性")
OUT["C2"] = dict(close_ts="2026-10-26T12:59:45+08:00", changed_after_close=changed,
                 watermark_listed=sorted(wm_paths), uncovered=uncovered)

# ══════════ C3 append-only（刪除行是否逐字重現） ══════════
print("\nC3 append-only：歷史檔刪除行是否【逐字重現】（換行正規化可接受；內容消失不可）")
gitdir = REPO if os.path.isdir(os.path.join(REPO, ".git")) else V11
commits = subprocess.run(["git", "-C", gitdir, "log", "--format=%H", "-10", "--",
                          "main_problems/m7_hadron_model_basis"], capture_output=True, text=True).stdout.split()
HISTPREFIX = ("01_problem_record", "SESSION_START")
rows, violations = [], []
for c in commits:
    out = subprocess.run(["git", "-C", gitdir, "show", "--format=", "--unified=0", c, "--",
                          "main_problems/m7_hadron_model_basis"], capture_output=True, text=True).stdout
    cur, dels, adds = None, [], []
    def is_live(d):   # 水位屬 live field（skill §1）：重生可刪舊行
        return bool(re.match(r"^\| \d+ \| `", d) or d.startswith("### §檔案閉合水位")
                    or d.startswith("> live field") or d.startswith("> 判「未變」")
                    or d.startswith("| # | 檔案 |") or d.startswith("|:-:|:-----"))

    def flush():
        if cur and dels:
            live = [d for d in dels if is_live(d)]
            lost = [d for d in dels if d.strip() and d not in adds and not is_live(d)]
            rows.append(dict(commit=c[:7], file=os.path.basename(cur), n_del=len(dels),
                             n_live=len(live), n_lost=len(lost)))
            if lost:
                violations.append(dict(commit=c[:7], file=os.path.basename(cur), lost=lost[:2]))
    for ln in out.splitlines():
        if ln.startswith("diff --git"):
            flush()
            cur = ln.split(" b/")[-1]
            base = os.path.basename(cur)
            if not (any(base.startswith(x) for x in HISTPREFIX) or "/docs/" in cur):
                cur = None
            dels, adds = [], []
        elif cur is not None:
            if ln.startswith("-") and not ln.startswith("---"):
                dels.append(ln[1:])
            elif ln.startswith("+") and not ln.startswith("+++"):
                adds.append(ln[1:])
    flush()
for r in rows:
    print(f"    {r['commit']} {r['file']:<34} 刪 {r['n_del']:>3} 行；其中水位 live-field {r['n_live']:>2}；未重現 {r['n_lost']}")
print(f"    違規（內容消失）：{violations or '無'}")
ok("C3 歷史檔零內容消失（append-only）", not violations, len(violations), 0,
   "刪除行若逐字出現於新增 ⟹ 屬尾行換行正規化（cat >> 撞無 EOL），非改寫")
OUT["C3"] = dict(rows=rows, violations=violations)

# ══════════ 總判 ══════════
print("\n" + "=" * 116)
print(f"結果：{NP[1]}/{NP[0]} PASS；findings = {len(OUT['findings'])}")
for f in OUT["findings"]:
    print(f"  🔴 {f['item']}：got {f['got']} | want {f['want']}  {f['note']}")
OUT["metrics"] = dict(b_hat=float(bhat), sig_b=float(sig_b), chi2=chi2, p=float(p3),
                      sum_gamma=float(s45), corr_pct=float(corr) * 100, mc_ratio=ratio,
                      cap=cap, g4M0=float(dev_xi), need_xi=float(need_xi),
                      uncovered=len(OUT["C2"]["uncovered"]))
OUT.update(npass=NP[1], ntotal=NP[0])
json.dump(OUT, open(os.path.join(M7, "data", "07_audit_full_chain.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/07_audit_full_chain.json")