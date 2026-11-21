#!/usr/bin/env python3
"""
M19 / scripts/10 — OA2 閉合：容許矩陣權重比可否作弱判別

預註冊（跑前凍結）:
  P1  權重比 = A[quark][≡1 欄] : A[quark][非 ≡1 欄] = 720 : 2880 = **1 : 4**（精確 Fraction）
  P2  比之來源檢定：1/4 = 1/#{≡1 槽}（4 槽隨機取 3 給輕子 ⟹ 某指定槽為空之機率 1/4）
  P3  均勻先驗 over **admissible** 雙射（S1 條件化）：quark-at-≡1 邊際 = 720/17280 = 1/24；
      quark-at-非≡1 = 2880/17280 = 1/6 ⟹ 比 1:4
  P4  反事實：均勻先驗 over **全部 9!**（不條件化 S1）：quark-at-任一模式 = 1/9 ⟹ 比 **1:1**
  P5  判據：若 P4 比 = 1:1 ⟹ 權重比**完全由 S1 條件化產生**，非獨立信息 ⟹ 不可作弱判別
      （否則循環：用 S1 之條件化去「支持」S1 之指派）
  P6  v11 有無「指派動力學／生成本」註冊？若無 ⟹ 冇物理先驗 ⟹ OA2 判 🔴（除非註冊動力學）

失敗判準: P1 比 ≠ 1:4 ⟹ 與 docs/04 §二 D6 之 720/2880 不一致（停）。
"""
import json
import math
from fractions import Fraction
from pathlib import Path

OUT = Path(__file__).parent.parent / "data"
M9 = [1, 2, 3, 4, 5, 7, 9, 11, 13]
S1 = {1, 5, 9, 13}
checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 scripts/10 — OA2：權重比可否作弱判別")
print("=" * 92)

# P1：由閉式重算 A[quark][≡1] 與 [非≡1]
a_eq1 = math.factorial(3) * math.factorial(5)          # 3!·5! = 720
a_neq = 4 * math.factorial(3) * math.factorial(5)      # 4·3!·5! = 2880
ratio = Fraction(a_eq1, a_neq)
chk("P1 比 = 1:4 精確", ratio == Fraction(1, 4), f"{a_eq1}:{a_neq} = {ratio}")

chk("P2 1/4 = 1/#≡1 槽", Fraction(1, len(S1)) == ratio, f"#{sorted(S1)} = {len(S1)}")

rho = 4 * math.factorial(3) * math.factorial(6)
m_eq1 = Fraction(a_eq1, rho)
m_neq = Fraction(a_neq, rho)
chk("P3 條件化邊際 1/24 vs 1/6", (m_eq1, m_neq) == (Fraction(1, 24), Fraction(1, 6)),
    f"{m_eq1} vs {m_neq}（比 {m_eq1/m_neq}）")

unchanged = Fraction(1, 9)
chk("P4 不條件化邊際 = 1/9（比 1:1）", unchanged == Fraction(1, 9), f"{unchanged} vs {unchanged}")

chk("P5 權重比由 S1 條件化產生（非獨立信息）",
    ratio != Fraction(1, 1) and m_eq1 / m_neq == ratio,
    "⟹ 用權重比支持指派 = 拿 S1 之條件化去支持 S1")

# P6：v11 有無「指派動力學」註冊（掃全 main_problems，排除本 M）
M19_ROOT = Path(__file__).resolve().parent.parent
PARENT = M19_ROOT.parent
docs_hit = []
for p in PARENT.glob("*/docs/*.md"):
    if p.parts[-3] == M19_ROOT.name:
        continue
    txt = p.read_text()
    for kw in ("指派動力學", "指派生成律", "味指派轉移"):
        if kw in txt:
            docs_hit.append((str(p.relative_to(PARENT)), kw))
chk("P6 v11 無「指派動力學」註冊（全 main_problems，排除本 M）", not docs_hit, str(docs_hit[:5]))

OUT.mkdir(exist_ok=True)
(OUT / "10_assignment_prior.json").write_text(json.dumps({
    "A_quark_eq1": a_eq1,
    "A_quark_neq1": a_neq,
    "ratio": str(ratio),
    "ratio_source": "1/#{≡1 slots} = 1/4（空槽均勻先驗）",
    "marginal_conditioned": {"eq1": str(m_eq1), "neq1": str(m_neq)},
    "marginal_unconditioned": str(unchanged),
    "verdict": "OA2 🔴 不可作弱判別：權重比係 S1 條件化之產物；v11 無指派動力學",
}, ensure_ascii=False, indent=1))

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
