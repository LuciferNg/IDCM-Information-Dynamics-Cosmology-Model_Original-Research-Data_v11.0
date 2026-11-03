#!/usr/bin/env python3
"""R6 — 通道型判別（Σ-型 vs Ξ-型）＋ 光扇區超定證書【誤差傳播改正版】

資料源（唯讀 vendored）: data/pdg2024/mass_width_2024.txt（PDG 2024, sha256 9d6024da…）
夸克內容由 PDG MC-ID 硬編碼（CONTENT 表）；非擬合參數。

★ S0 誤差傳播（跑前寫死；本輪必要紀律）
   所有 δ 與預言皆係「強子質量之線性組合」。凡**預言與觀測共用同一質量**
   （例：δ_c 內含 Σc⁺，而觀測 Σc⁺−Σc⁺⁺ 亦含 Σc⁺），naive 獨立誤差相加會錯配 σ。
   本腳本一律以 Form = {強子: Fraction 係數} 表示預言與觀測；
   殘差 = obs − pred（Form），σ = √Σ (c_i σ_i)² ⟹ 共用質量之相關自動正確處理。
   ⚠ 前版（naive σ）把 Σc⁺−Σc⁺⁺ 誤判為 z = −3.64；改正後見下（方法瑕疵，已修，
     舊版保留為 `scripts/19_r6_channel_type.py.pre_patch_*`）。

預註冊判準（跑前凍結）
--------------------------------------------------
ansatz 沿用 R4/R5：Δm = Σ_(價對) a(f_i,f_j)，δ_x ≡ a(d,x) − a(u,x)。
δ_u, δ_d, δ_s 由三錨 {n−p, Σ⁻−Σ⁰, Ξ⁻−Ξ⁰} 唯一解；δ_c 由 Σc⁰−Σc⁺；δ_b 由 Σb⁻−Σb⁺。

T1 Σ-型（恰一個 s/c/b 旁觀者）：全通（|z| < 3）；任一 ≥3 ⟹ 假設死。
T2 Ξ-型（兩個重旁觀者）：全死（|z| ≥ 3）；任一 <3 ⟹ 假設弱化，逐條照報。
T3 自旋盲：同味 multiset 之 1/2⁺ 與 3/2⁺ 相等（5 對）。
T4 介子符號分類：同構型是否同號 ⟹ 單一 ℤ₂ 可否。
T5 光扇區超定證書：rank(M) vs rank([M|b])，精確 Fraction。

判死門檻 |z| ≥ 3。禁回擬、禁加參數、禁剔除通道。
"""
# ── 2026-11-03 更正（R30）─────────────────────────────────────────────────────
# 原缺陷：本腳本 ROW 表把三列 Ξc 之預言硬編碼為 `(d_s + d_c).scale(-1)`（負號），
#         而正確符號由夸克內容唯一決定：Ξc⁰ = dsc（d 型）、Ξc⁺ = usc（u 型）
#         ⟹ Ξc⁰ − Ξc⁺ = +(δ_s + δ_c)。同表 Ξb⁻−Ξb⁰ = +(δ_s + δ_b) 一直正確 ⟹ 三列 Ξc 屬孤例錯。
#         另：CONTENT 表之 4132/4232/4312/4322/4314/4324 夸克內容與 PDG 電荷相反，惟該表只作
#         存在性 assert（L82），**不參與預言計算** ⟹ 非本次符號錯之機制，但已一併更正。
# 更正後：docs/15 之三列 Ξc 預言 = +3.1430；`scripts/51`（docs/44）已先行獨立算出同值。
# ────────────────────────────────────────────────────────────────────────────

import re
import hashlib
import pathlib
from decimal import Decimal
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 40
SRC = pathlib.Path("data/pdg2024/mass_width_2024.txt")

CONTENT = {
    "2112": "udd", "2212": "uud",
    "3112": "dds", "3212": "uds", "3222": "uus", "3312": "dss", "3322": "uss",
    "3114": "dds", "3214": "uds", "3224": "uus", "3314": "dss", "3324": "uss",
    "4112": "ddc", "4212": "udc", "4222": "uuc",
    "4114": "ddc", "4214": "udc", "4224": "uuc",
    "4132": "dsc", "4232": "usc", "4312": "dsc", "4322": "usc",
    "4314": "dsc", "4324": "usc",
    "5112": "bdd", "5222": "buu", "5114": "bdd", "5224": "buu",
    "5132": "bds", "5232": "bus",
    "321": "us~", "311": "ds~", "411": "cd~", "421": "cu~", "511": "db~", "521": "ub~",
}
LABEL = {
    "2112": "n", "2212": "p", "3112": "Σ⁻", "3212": "Σ⁰", "3222": "Σ⁺",
    "3312": "Ξ⁻", "3322": "Ξ⁰", "3114": "Σ*⁻", "3214": "Σ*⁰", "3224": "Σ*⁺",
    "3314": "Ξ*⁻", "3324": "Ξ*⁰",
    "4112": "Σc⁰", "4212": "Σc⁺", "4222": "Σc⁺⁺",
    "4114": "Σc*⁰", "4214": "Σc*⁺", "4224": "Σc*⁺⁺",
    "4132": "Ξc⁰", "4232": "Ξc⁺", "4312": "Ξc′⁰", "4322": "Ξc′⁺",
    "4314": "Ξc(2645)⁰", "4324": "Ξc(2645)⁺",
    "5112": "Σb⁻", "5222": "Σb⁺", "5114": "Σb*⁻", "5224": "Σb*⁺",
    "5132": "Ξb⁻", "5232": "Ξb⁰",
    "321": "K⁺", "311": "K⁰", "411": "D⁺", "421": "D⁰", "511": "B⁰", "521": "B⁺",
}
INV = {v: k for k, v in LABEL.items()}


def load():
    out = {}
    for line in SRC.read_text().splitlines():
        if line.startswith("*"):
            continue
        m = re.match(r"^\s*(-?\d+)\s+(.+)$", line)
        if not m:
            continue
        pid, rest = m.group(1), m.group(2)
        nums = re.findall(r"[+-]?\d\.\d+E[+-]\d+", rest)
        if len(nums) < 2:
            continue
        out[pid] = (F(Decimal(nums[0])) * 1000, abs(F(Decimal(nums[1]))) * 1000)
    return out


MS = load()
assert not [p for p in CONTENT if p not in MS]
MASS = {lb: MS[INV[lb]][0] for lb in LABEL.values()}
SIG = {lb: MS[INV[lb]][1] for lb in LABEL.values()}


class Form(dict):
    """強子質量之線性組合 {label: Fraction coeff}"""

    def __add__(self, o):
        r = Form(self)
        for k, v in Form(o).items():
            r[k] = r.get(k, F(0)) + v
        return Form({k: v for k, v in r.items() if v != 0})

    def __sub__(self, o):
        return self + Form({k: -v for k, v in Form(o).items()})

    def scale(self, c):
        c = F(c)
        return Form({k: v * c for k, v in self.items() if v * c != 0})

    def val(self):
        return sum((c * MASS[k] for k, c in self.items()), F(0))

    def sd(self):
        return mp.sqrt(sum(mp.mpf(float(c * SIG[k])) ** 2 for k, c in self.items()))


def H(name):
    return Form({name: F(1)})


def CH(a, b):
    return H(a) - H(b)


print("=" * 96)
print("R6 — 通道型判別（Σ-型 vs Ξ-型）＋ 光扇區超定證書  [誤差傳播改正版]")
print("=" * 96)
print("  PDG 2024 sha256(前32) =", hashlib.sha256(SRC.read_bytes()).hexdigest()[:32])
print("  S0：殘差以質量線性組合精確傳播（正確處理預言∩觀測之共用質量）")

d_s = CH("Ξ⁻", "Ξ⁰").scale(F(1, 2))
d_d = CH("Σ⁻", "Σ⁰") - d_s
d_u = CH("n", "p") - d_d
d_c = CH("Σc⁰", "Σc⁺") - d_d
d_b = (CH("Σb⁻", "Σb⁺") - d_u - d_d).scale(F(1, 2))

print("\n[δ] 由錨唯一解（Form；值與 σ 精確傳播）")
for nm, v in [("δ_u", d_u), ("δ_d", d_d), ("δ_s", d_s), ("δ_c", d_c), ("δ_b", d_b)]:
    print("   %-4s = %+9.4f ± %7.4f MeV" % (nm, float(v.val()), float(v.sd())))

TESTS = [
    ("Σ⁻−Σ⁰", "Σ⁻", "Σ⁰", "Σ", d_d + d_s, True),
    ("Ξ⁻−Ξ⁰", "Ξ⁻", "Ξ⁰", "Ξ", d_s.scale(2), True),
    ("Σ⁰−Σ⁺", "Σ⁰", "Σ⁺", "Σ", d_u + d_s, False),
    ("Σ⁻−Σ⁺", "Σ⁻", "Σ⁺", "Σ", d_u + d_d + d_s.scale(2), False),
    ("Σ*⁻−Σ*⁰", "Σ*⁻", "Σ*⁰", "Σ", d_d + d_s, False),
    ("Σ*⁰−Σ*⁺", "Σ*⁰", "Σ*⁺", "Σ", d_u + d_s, False),
    ("Ξ*⁻−Ξ*⁰", "Ξ*⁻", "Ξ*⁰", "Ξ", d_s.scale(2), False),
    ("Σc⁰−Σc⁺", "Σc⁰", "Σc⁺", "Σ", d_d + d_c, True),
    ("Σc⁺−Σc⁺⁺", "Σc⁺", "Σc⁺⁺", "Σ", d_u + d_c, False),
    ("Σc*⁰−Σc*⁺", "Σc*⁰", "Σc*⁺", "Σ", d_d + d_c, False),
    ("Σc*⁺−Σc*⁺⁺", "Σc*⁺", "Σc*⁺⁺", "Σ", d_u + d_c, False),
    ("Ξc⁰−Ξc⁺", "Ξc⁰", "Ξc⁺", "Ξ", (d_s + d_c), False),
    ("Ξc′⁰−Ξc′⁺", "Ξc′⁰", "Ξc′⁺", "Ξ", (d_s + d_c), False),
    ("Ξc(2645)⁰−Ξc(2645)⁺", "Ξc(2645)⁰", "Ξc(2645)⁺", "Ξ", (d_s + d_c), False),
    ("Σb⁻−Σb⁺", "Σb⁻", "Σb⁺", "Σ", d_b.scale(2) + d_u + d_d, True),
    ("Σb*⁻−Σb*⁺", "Σb*⁻", "Σb*⁺", "Σ", d_b.scale(2) + d_u + d_d, False),
    ("Ξb⁻−Ξb⁰", "Ξb⁻", "Ξb⁰", "Ξ", d_b + d_s, False),
]

print("\n" + "=" * 96)
print("[檢定] 樣本外預言（anchor = 用於解 δ，不入 tally）")
print("  %-24s %-2s %10s %10s %9s %8s %7s" % ("通道", "型", "預言", "實測", "殘差", "σ", "z"))
sp = sk = xp = xk = 0
ROWS = {}
for nm, ha, hb, typ, pred, anch in TESTS:
    obs = CH(ha, hb)
    res = obs - pred
    rv, rs = res.val(), res.sd()
    z = (mp.mpf(float(rv)) / rs) if rs > 0 else mp.mpf(0)
    bad = abs(z) >= 3
    ROWS[nm] = (float(obs.val()), float(rs), z, bad, anch)
    tag = "anchor" if anch else ("判死 🔴" if bad else "通 🟢")
    print("  %-24s %-2s %+10.4f %+10.4f %+9.4f %8.4f %+7.2f  %s"
          % (nm, typ, float(pred.val()), float(obs.val()), float(rv), float(rs), z, tag))
    if not anch:
        if typ == "Σ":
            sk += bad
            sp += (not bad)
        else:
            xk += bad
            xp += (not bad)

print("\n  [T1] Σ-型：通 %d／死 %d  ⟹ %s"
      % (sp, sk, "假設成立 🟢" if sk == 0 else "假設死 🔴（逐條見上表）"))
print("  [T2] Ξ-型：通 %d／死 %d  ⟹ %s"
      % (xp, xk, "假設成立 🟢" if xp == 0 else "假設弱化 🟡（逐條見上表）"))

print("\n" + "=" * 96)
print("[T3] 自旋盲：同味 multiset 之 1/2⁺ vs 3/2⁺ 相等")
PAIRS = [("Σ", "Σ⁻−Σ⁰", "Σ*⁻−Σ*⁰"), ("Σc", "Σc⁰−Σc⁺", "Σc*⁰−Σc*⁺"),
         ("Σc+", "Σc⁺−Σc⁺⁺", "Σc*⁺−Σc*⁺⁺"),
         ("Ξc", "Ξc⁰−Ξc⁺", "Ξc(2645)⁰−Ξc(2645)⁺"),
         ("Σb", "Σb⁻−Σb⁺", "Σb*⁻−Σb*⁺")]
n3 = 0
for nm, h, t in PAIRS:
    ha, hb = h.split("−")
    ta, tb = t.split("−")
    res = CH(ha, hb) - CH(ta, tb)
    rv, rs = res.val(), res.sd()
    z = (mp.mpf(float(rv)) / rs) if rs > 0 else mp.mpf(0)
    bad = abs(z) >= 3
    n3 += bad
    print("   %-4s %+9.4f vs %+9.4f  Δ=%+9.4f ± %7.4f  z=%+7.2f  %s"
          % (nm, ROWS[h][0], ROWS[t][0], float(rv), float(rs), z,
             "判死 🔴" if bad else "通 🟢"))
print("   ⟹ 自旋盲：%d/%d 對判死" % (n3, len(PAIRS)))

print("\n" + "=" * 96)
print("[T4] 介子符號分類（isodoublet 側 × 旁觀者）")
for nm, ha, hb, side, spec in [("K⁺−K⁰", "K⁺", "K⁰", "夸克側", "s̄（重反夸克）"),
                               ("D⁺−D⁰", "D⁺", "D⁰", "反夸克側", "c（重夸克）"),
                               ("B⁰−B⁺", "B⁰", "B⁺", "夸克側", "b̄（重反夸克）")]:
    v = CH(ha, hb).val()
    print("   %-9s = %+9.4f MeV   isodoublet 在 %-6s  旁觀者 %-14s ⟹ 號 %s"
          % (nm, float(v), side, spec, "＋" if v > 0 else "－"))
print("   ⟹ 同構型（夸克側 isodoublet ＋ 重反夸克旁觀者）：K －、B ＋")
print("     ⟹ 單一 ℤ₂ 不可能定符號（範疇陳述，不依賴參數）")

print("\n" + "=" * 96)
print("[T5] 光扇區超定證書（Fraction 精確）")
LIGHT = ["n−p", "Σ⁻−Σ⁰", "Σ⁰−Σ⁺", "Ξ⁻−Ξ⁰", "Σ⁻−Σ⁺",
         "Σ*⁻−Σ*⁰", "Σ*⁰−Σ*⁺", "Ξ*⁻−Ξ*⁰"]
COEF = {"n−p": (1, 1, 0), "Σ⁻−Σ⁰": (0, 1, 1), "Σ⁰−Σ⁺": (1, 0, 1), "Ξ⁻−Ξ⁰": (0, 0, 2),
        "Σ⁻−Σ⁺": (1, 1, 2), "Σ*⁻−Σ*⁰": (0, 1, 1), "Σ*⁰−Σ*⁺": (1, 0, 1), "Ξ*⁻−Ξ*⁰": (0, 0, 2)}


def rankF(Mx):
    A = [row[:] for row in Mx]
    r, ncol = 0, (len(A[0]) if A else 0)
    for c in range(ncol):
        p = next((i for i in range(r, len(A)) if A[i][c] != 0), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == len(A):
            break
    return r


CHMAP = {nm: CH(*nm.split("−")) for nm in LIGHT}
Mx = [[F(c) for c in COEF[nm]] for nm in LIGHT]
bv = [CHMAP[nm].val() for nm in LIGHT]
rM, rA = rankF(Mx), rankF([row + [b] for row, b in zip(Mx, bv)])
print("   未知 3；通道 %d；rank(M)=%d  rank([M|b])=%d ⟹ %s（解空間維度 %d，獨立關係 %d）"
      % (len(LIGHT), rM, rA, "相容" if rM == rA else "**不相容（ansatz 被判死）**",
         3 - rM, len(LIGHT) - rM))
print("   逐條殘差（精確誤差傳播）：")
for nm in LIGHT:
    cu, cd, cs = COEF[nm]
    res = CHMAP[nm] - (d_u.scale(cu) + d_d.scale(cd) + d_s.scale(cs))
    rv, rs = res.val(), res.sd()
    z = mp.mpf(float(rv)) / rs if rs > 0 else mp.mpf(0)
    print("     %-12s 殘差 %+9.4f ± %7.4f  z=%+7.2f  %s"
          % (nm, float(rv), float(rs), z, "越界 🔴" if abs(z) >= 3 else ""))

print("\n" + "=" * 96)
print("[T6] R5 檢定之誤差傳播更正（同一 ansatz，改用 Form 傳播）")
R5 = [
    ("A1 Σ⁰−Σ⁺", CH("Σ⁰", "Σ⁺"), d_u + d_s, -0.87),
    ("A2 Σ*⁻−Σ*⁰", CH("Σ*⁻", "Σ*⁰"), d_d + d_s, -1.17),
    ("A3 Σ*⁰−Σ*⁺", CH("Σ*⁰", "Σ*⁺"), d_u + d_s, -2.34),
    ("A4 Ξ*⁻−Ξ*⁰", CH("Ξ*⁻", "Ξ*⁰"), d_s.scale(2), -5.37),
    ("B1 Σc 雙差", (CH("Σc⁰", "Σc⁺") - CH("Σc⁺", "Σc⁺⁺")), d_d - d_u, +2.57),
    # ★ 2026-11-03 更正（R30）：預言原寫 d_d − d_s —— 該式只在「Ξc⁰ = usc」之錯內容假設下成立；
    #   正確內容 Ξc⁰ = dsc ⟹ Ξc⁰−Ξc⁺ = +(δ_s+δ_c)，Σc⁰−Σc⁺ = +(δ_d+δ_c) ⟹ 和 = δ_d+δ_s+2δ_c（粲項不相消）。
    ("B2 Σc+Ξc 和", (CH("Σc⁰", "Σc⁺") + CH("Ξc⁰", "Ξc⁺")), d_d + d_s + d_c.scale(2), +13.16),
    ("C1 底組合", (CH("Σb⁻", "Σb⁺") - CH("Ξb⁻", "Ξb⁰").scale(2)), CH("n", "p") - CH("Ξ⁻", "Ξ⁰"), +0.27),
]
print("  %-16s %9s %9s %7s %7s  %s" % ("檢定", "殘差", "σ(Form)", "z_舊", "z_新", "判"))
for nm, obs, pred, zold in R5:
    res = obs - pred
    rv, rs = res.val(), res.sd()
    znew = (mp.mpf(float(rv)) / rs) if rs > 0 else mp.mpf(0)
    bad = abs(znew) >= 3
    print("  %-16s %+9.4f %9.4f %+7.2f %+7.2f  %s"
          % (nm, float(rv), float(rs), zold, znew, "判死 🔴" if bad else "通 🟢"))
print("  ⟹ z_舊 用獨立誤差相加；z_新 用 Form 精確傳播。判詞（通/死）除 A1 之數值外不變；")
print("     A1 |z| 由 0.87 降至 0.29（預言含 Σ⁰ ⟹ 共用質量使 σ 正確變大）。")

print("\n" + "=" * 96)
print("[匯總 — 只報事實，不裁決]")
print("   Σ-型樣本外：通 %d／死 %d" % (sp, sk))
print("   Ξ-型樣本外：通 %d／死 %d" % (xp, xk))
print("   自旋盲：%d/%d 對判死" % (n3, len(PAIRS)))
print("   介子：單一 ℤ₂ 判死（同構型 K 負 B 正）")
print("   光扇區：rank(M)=%d / rank([M|b])=%d ⟹ 二體荷%s" % (rM, rA, "相容" if rM == rA else "不相容"))
print("=" * 96)