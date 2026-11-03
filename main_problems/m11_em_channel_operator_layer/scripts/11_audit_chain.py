#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M11 scripts/11 — 全鏈審計（idcm-audit-verify §3.1 六項 ＋ §7.3 L2/L3 納入）
================================================================================
只讀（除寫 data/11_audit_chain.json ＋ data/run_11_audit_chain.log）。

檢查：
  ① 腳本集對稱：磁碟 scripts/*.py ↔ 文檔所列 ⇒ 兩側差集須空
  ② 逆向綁定：每個 headline 數字 → 有磁碟 producer（script）且其 log/JSON 內出現
  ③ 路徑衛生：scripts/*.py 內真寫入模式指向 /tmp/ ⇒ 0 命中
  ④ 產物綁定：data/*.json ↔ 產者腳本（NN 前綴）；腳本宣告之輸檔存在（無孤兒）
  ⑤ inline 幽靈：docs 內 ≥5 位有效數字，若不在任何 data/*.json 或 data/run_*.log，
     且不在 ALLOW（宣告之外部輸入常數）⇒ WARN 列出（須逐項申報）
  ⑥ driver：scripts/run_chain.sh 存在且列齊各輪
  E1／E2（§7.3）：l2_independent_verify.py、l3_doc_crosscheck.py 必須存在並納入
"""
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS, DATA, DOCS = ROOT / "scripts", ROOT / "data", ROOT / "docs"
LOG, RES = [], {}
NFAIL = NWARN = 0

HEADLINE = {  # ② headline 數字 → 期望 producer 前綴
    "C9 K 之 t = −0.3601": ("01", "0.3601"),
    "k₀ 3.815262": ("01", "3.815262"),
    "⟨O⟩_K 2.4213": ("01", "2.4213"),
    "1040 候選 0 PASS": ("02", "1040"),
    "18468 候選 0 PASS": ("03", "18468"),
    "R4 D 飽和 1.908984": ("04", "1.908984"),
    "k₀^crit D 2.5425": ("05", "2.5425"),
    "k₀^crit A 4.0730": ("05", "4.072960681856742"),
    "a^crit 0.26987": ("05", "0.26986741279064763"),
    "κ=0.300 帶": ("06", "0.300"),
    "em(0.79) 2.2773": ("07", "2.2773"),
    "G0 δM 0.01811": ("08", "0.01811"),
    "路線 C χ²/dof 3154": ("09", "3154"),
    "退化 1.9787": ("10", "1.9787"),
}
# §3.4 硬規則：每個對不上 data 嘅數字都要明示理由
ALLOW_EXPLICIT = {
    "0.3663": "R6 前 k₀=3.750839 之 K t 歷史列（docs/01 雙欄並列）",
    "0.1898": "同上 —— Σ⁰−Σ⁺ t 之舊 k₀ 版",
    "1.1437": "同上 —— Ξ⁻−Ξ⁰ t 之舊 k₀ 版",
    "0.5991": "同上 —— Σ⁻−Σ⁰ t 之舊 k₀ 版",
    "0.5422": "同上 —— π⁺−π⁰ t 之舊 k₀ 版",
    "0.9484": "R4 模型 C 之 S--S+ t（em/k₀_old；em 屬模型預言）",
    "0.7239": "R4 模型 A 之 S--S+ t（em/k₀_old）",
    "0.5666": "R4 模型 A/B 之 t（em/k₀_old）",
    "0.5163": "R4 模型 C 之 D+-D0 t（em/k₀_old）",
    "0.4306": "R4 模型 B 之 t（em/k₀_old）",
    "0.0245": "L2 首版單位陷阱之**假值**（已作廢；真值 0.018115，見 §十三 傳播表 #4）",
    "2.2217": "L2 首版 σ_O 加性傳播之**掩錯值**（已作廢；真值 2.108263，見 §14.4）",
    "5.9100": "同上 —— k₀^crit(C) 之掩錯值（已作廢；真值 5.335822，見 §14.4）",
    "4.0730": "k₀^crit(A) 之 4 位顯示（真值 4.07296068，見 §14.1）",
    "5.3358": "k₀^crit(C) 之 4 位顯示（真值 5.33582226）",
    "0.1038": "σ_O(Ξ) 正典之 4 位顯示（真值 0.103777；σ_Δm 不入）",
    "0.132744": "σ_O(Ξ) 補 σ_Δm 之值（L2 full 變體 sigO_XimmXi0）",
    "0.19666": "相對式錯值 |O|σ_Δm/(cX)（§14.6 拆解用）",
    "0.09055": "絕對式 σ_Δm/(cX)（σ_Δm=0.2318 時；§14.6）",
    "0.0798": "σ_X 地板 O·σ_X/X 之 Ξ 值（§14.6 達標分析）",
    "0.0435": "達標所需 σ_em 上限（§14.6）",
    "0.0386": "達標所需 σ_X 上限（§14.6）",
    "1.279": "§14.6 之比值 σ_O(full)/σ_O(正典) = 0.132744/0.103777（兩值皆在 L2）",
    "3899.73": "外部常數 PDG `ΔM²_K = −3899.73 MeV²`（docs/09 第 9 條替代渠道；非本記錄計算量）",
    "1.79": "M5 docs/56 所述 ±1.79%（其基準）；M11 實算 −1.760%（data/13 之 shift_pct），差屬取位",
    "0.001": "percentage 顯示（`0.001%` = 1e−5，docs/07 之 `k=2.2361≈√5` 偏差）；非 data 量",
    "0.2076": "M12 來件之導式錯值（誤乘 |O|；已於 §14.6 拆解）",
    "0.2318": "M12 來件之 σ_E（其表列值，非 M11 量）",
    "0.1363": "用 M12 之 σ_E=0.2318 重算之正確值（§14.6）",
}
ALLOW_PAT = re.compile(r"^(20\d\d|\d{4}\.\d{4,5}$|2\.2361|1\.166667|1\.186188|0\.9813|0\.5136262|0\.514|1[0-9]\.\d+|137\.035999|197\.3269804|0\.1973269804|"
                       r"2\.99792458|6\.02214076|1\.602176634|"      # 物理常數／年份
                       r"0\.37742216816391858835|522\.8283|2\.56|0\.094|"  # 校準與 X
                       r"4\.4492|1\.166158|1\.187051|"                # M5 對接量（上游）
                       r"0\.78|0\.79|0\.73|0\.6|0\.17|0\.10|0\.15|"   # FLAG ε 分支
                       r"0\.30|0\.0007884|4\.163|0\.2975|0\.305|2\.04)")  # κ 掃描輸出


def chk(tag, cond, detail="", warn=False):
    global NFAIL, NWARN
    ok = bool(cond)
    if not ok:
        if warn:
            NWARN += 1
        else:
            NFAIL += 1
    line = f"[{'PASS' if ok else ('WARN' if warn else 'FAIL')}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "warn": warn and not ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def listed_scripts():
    """由 01_problem_record／README／docs 抽出被引用之 scripts/NN"""
    names = set()
    for p in [ROOT / "01_problem_record.md", ROOT / "README.md", ROOT / "SESSION_START.md"] + sorted(DOCS.glob("*.md")):
        if not p.exists():
            continue
        txt = p.read_text(encoding="utf-8", errors="ignore")
        names |= set(re.findall(r"scripts/(\d{2})[a-z_0-9]*\.py", txt))
    return names


def main():
    hdr("① 腳本集對稱（磁碟 ↔ 文檔）")
    disk = {p.name[:2] for p in sorted(SCRIPTS.glob("[01][0-9]_*.py"))}
    disk |= {p.name[:2] for p in sorted(SCRIPTS.glob("2[0-9]_*.py"))}
    listed = listed_scripts()
    print(f"  磁碟 = {sorted(disk)}")
    print(f"  文檔 = {sorted(listed)}")
    chk("① 兩側差集空（磁碟有但文檔無）", not (disk - listed), f"{sorted(disk - listed)}")
    chk("① 兩側差集空（文檔有但磁碟無）", not (listed - disk), f"{sorted(listed - disk)}", warn=True)

    hdr("② 逆向綁定（headline 數字 → producer）")
    for label, (pre, token) in HEADLINE.items():
        prod = sorted(SCRIPTS.glob(f"{pre}_*.py"))
        if not prod:
            chk(f"② {label} → scripts/{pre}_*", False, "無 producer 檔")
            continue
        stem = prod[0].stem
        hay = ""
        for cand in (DATA / f"run_{stem}.log",):
            if cand.exists():
                hay += cand.read_text(encoding="utf-8", errors="ignore")
        for j in sorted(DATA.glob(f"{pre}_*.json")):
            hay += j.read_text(encoding="utf-8", errors="ignore")
        chk(f"② {label} → {stem}", token in hay, f"token={token!r} 命中={token in hay}")

    hdr("③ 路徑衛生（/tmp/ 真寫入模式）")
    pat = re.compile(r"(open\(\s*['\"]/tmp/|Path\(\s*['\"]/tmp/|>\s*/tmp/|write_text\(\s*['\"]?/tmp/)")
    hits = []
    for p in sorted(SCRIPTS.glob("*.py")):
        t = p.read_text(encoding="utf-8", errors="ignore")
        if pat.search(t):
            hits.append(p.name)
    chk("③ 無 /tmp/ 寫入", not hits, f"{hits}")

    hdr("④ 產物綁定（data ↔ script，無孤兒）")
    orphans, noprod = [], []
    for j in sorted(DATA.glob("[01][0-9]_*.json")) + sorted(DATA.glob("2[0-9]_*.json")):
        pre = j.name[:2]
        if not list(SCRIPTS.glob(f"{pre}_*.py")):
            orphans.append(j.name)
    for p in sorted(SCRIPTS.glob("[01][0-9]_*.py")) + sorted(SCRIPTS.glob("2[0-9]_*.py")):
        pre = p.name[:2]
        if not list(DATA.glob(f"{pre}_*.json")):
            noprod.append(p.name)
    chk("④ 無孤兒 JSON", not orphans, f"{orphans}")
    chk("④ 每支腳本有 JSON 產物", not noprod, f"{noprod}", warn=True)

    hdr("⑤ inline 幽靈（docs ≥5 位數字，不在 data 任何 json/log）")
    hay = ""
    for f in sorted(DATA.glob("*.json")) + sorted(DATA.glob("run_*.log")):
        hay += f.read_text(encoding="utf-8", errors="ignore")
    num = re.compile(r"(?<!\w\.)(\d+\.\d{4,})(?![\w.])")
    hay_nums = [float(x) for x in re.findall(r"\d+\.\d+", hay)]
    ghosts = {}
    for d in sorted(DOCS.glob("*.md")) + [ROOT / "01_problem_record.md"]:
        if not d.exists():
            continue
        for m in num.finditer(d.read_text(encoding="utf-8", errors="ignore").replace("\u2212", "-")):
            tok = m.group(1)
            if tok in hay or ALLOW_PAT.match(tok) or tok in ALLOW_EXPLICIT:
                continue
            k = len(tok.split(".")[1])
            v0 = float(tok)
            if any(abs(round(v, k) - v0) < 10 ** (-(k + 1)) for v in hay_nums):
                continue   # 四捨五入即可對上 data 內某值 ⟹ 非幽靈
            val = abs(v0)
            if val in (0.0,) or len(tok.split(".")[1]) < 4:
                continue
            ghosts.setdefault(tok, set()).add(d.name)
    for k, why in sorted(ALLOW_EXPLICIT.items()):
        print(f"    ✓ 已申報 {k:10s} — {why}")
    if ghosts:
        for k in sorted(ghosts, key=lambda x: -float(x)):
            print(f"    ~ {k:22s} {sorted(ghosts[k])}")
    chk("⑤ 零未申報幽靈", not ghosts, f"{len(ghosts)} 個（見上；須逐項入 ALLOW 或補 producer）", warn=True)

    hdr("⑥ driver（run_chain.sh）")
    drv = SCRIPTS / "run_chain.sh"
    ok = drv.exists()
    if ok:
        t = drv.read_text(encoding="utf-8")
        rounds = re.findall(r"run\s+\S+\s+\S+", t)
        ok = len(rounds) >= 11 and "OMP_NUM_THREADS=1" in t
        detail = f"{len(rounds)} 輪、OMP 鎖={'OMP_NUM_THREADS=1' in t}"
    else:
        detail = "缺 scripts/run_chain.sh"
    chk("⑥ driver 齊備", ok, detail)

    hdr("⑦ log 內容健全（無 traceback／路徑錯）")
    bad = []
    for lg in sorted(DATA.glob("run_*.log")):
        if "pre_patch" in lg.name:
            continue
        t = lg.read_text(encoding="utf-8", errors="ignore")
        for pat2 in ("Traceback (most recent call last)", "FileNotFoundError",
                     "No such file or directory", "ModuleNotFoundError", "SyntaxError"):
            if pat2 in t:
                bad.append(f"{lg.name}:{pat2}")
                break
    chk("⑦ 所有 run_*.log 無 traceback／路徑錯", not bad, f"{bad}")

    hdr("E §7.3 L2／L3 納入")
    for name, need in (("l2_independent_verify.json", "L2"), ("l3_doc_crosscheck.json", "L3")):
        p = DATA / name
        good = p.exists()
        if good:
            j = json.loads(p.read_text(encoding="utf-8"))
            summ = j.get("_summary", {})
            good = summ.get("fail", 1) == 0
            detail = f"checks={summ.get('checks')} fail={summ.get('fail')}"
        else:
            detail = "缺 JSON"
        chk(f"E {need} 綠燈", good, detail)

    hdr("匯總")
    n = len([k for k in RES])
    print(f"檢查項 = {n}，FAIL = {NFAIL}，WARN = {NWARN}")
    RES["_summary"] = {"checks": n, "fail": NFAIL, "warn": NWARN}
    (DATA / "11_audit_chain.json").write_text(json.dumps(RES, ensure_ascii=False, indent=2),
                                              encoding="utf-8")
    print(f"\n寫出: {DATA/'11_audit_chain.json'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())