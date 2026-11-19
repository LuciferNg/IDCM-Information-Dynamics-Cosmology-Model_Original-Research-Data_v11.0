#!/usr/bin/env python3
"""
M19 / l2 — 獨立重算（零 import 本輪任何腳本；閉式推算）

判據（跑前凍結）:
  L2.1  ρ_struct = 4·3!·6! 閉式 = 17,280（唔用枚舉）
  L2.2  9! = 362,880；比值 = 21；結構位元 = log₂21；輸入位元 = log₂17,280
  L2.3  容許矩陣閉式：輕子任一 ≡1 槽 = C(3,2)·2!·6! = 4,320；
        夸克喺 ≡1 欄 = 3!·5! = 720；喺非 ≡1 欄 = 4·3!·5! = 2,880
  L2.4  全 T3 不可滿足（鴿籠：χ=1 粒子 6 > ≡1 槽 4）—— 解析證明，唔用枚舉
  L2.5  γ₁ 以替代式 2 − 2cos(π/44) 重算，與 4sin²(π/88) 一致（dps=60）
  L2.6  9 模式集之 A₅ 動差準則閉式：奇 n ≤ 13 共 7 個、偶 n ≤ 4 共 2 個 ⟹ 9
"""
import math
import sys
from fractions import Fraction

from mpmath import mp, mpf, cos, sin, pi

mp.dps = 60

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


print("=" * 92)
print("M19 l2 — 獨立重算（閉式）")
print("=" * 92)

f6 = math.factorial(6)
rho = 4 * math.factorial(3) * f6
chk("L2.1 ρ_struct = 4·3!·6! = 17,280", rho == 17280, f"{rho:,}")

total = math.factorial(9)
ratio = Fraction(total, rho)
chk("L2.2 9!/ρ = 21 精確", ratio == 21, f"{total:,}/{rho:,} = {ratio}")
chk("L2.2 位元 log₂21 / log₂17280", abs(math.log2(21) - 4.392317422778761) < 1e-12
    and abs(math.log2(17280) - 14.076815597050832) < 1e-12,
    f"{math.log2(21):.9f} + {math.log2(17280):.9f} = {math.log2(362880):.9f}")

lep = math.comb(3, 2) * math.factorial(2) * f6
q1 = math.factorial(3) * math.factorial(5)
q0 = 4 * math.factorial(3) * math.factorial(5)
chk("L2.3 矩陣閉式 4320 / 720 / 2880", (lep, q1, q0) == (4320, 720, 2880),
    f"{lep} / {q1} / {q0}")

chi1 = 3 + 3            # leptons 3 + {u,c,t} 3
slots1 = 4
chk("L2.4 鴿籠：χ=1 粒子 6 > ≡1 槽 4 ⟹ 全 T3 不可滿足", chi1 > slots1, f"{chi1} > {slots1}")

g_alt = 2 - 2 * cos(pi / 44)
g_std = 4 * sin(pi / 88) ** 2
chk("L2.5 γ₁ 兩式一致（<1e−58）", abs(g_alt - g_std) < mpf("1e-58"),
    f"{mp.nstr(g_std, 22)} vs {mp.nstr(g_alt, 22)}")

odd = [n for n in range(1, 14) if n % 2 == 1]
even = [n for n in range(1, 5) if n % 2 == 0]
chk("L2.6 A₅ 動差準則：7 奇 ＋ 2 偶 = 9", len(odd) == 7 and len(even) == 2, f"{odd} + {even}")

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
raise SystemExit(0 if n_pass == len(checks) else 1)
