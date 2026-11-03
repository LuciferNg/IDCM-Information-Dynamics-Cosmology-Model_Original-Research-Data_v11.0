#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
scripts/03_variants_and_lookelsewhere.py — M16 R3 + R4
輸入：data/run_02_scan_all_500.txt（階 ≤ 500 之全部 92,570 群之 irrep 維數重數）
輸出：
  §A 準則本身在全集之命中（應 = {A5}）
  §B 修飾詞敏感度：去「相異」/ 去「非平凡」後之命中集
  §C ★ look-elsewhere：以「Σ(DNT) = c」掃 c = 1..40，(∏(DNT) = |G|) ∧ (Σ = c) 之命中
  §D 準則族 (selector × feature × target) 之單條／聯合（兩條）鑑別力統計
全部以整數運算；bitmask 用 numpy packbits 建構。
"""
import os, sys, itertools, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
SCAN = os.path.join(DATA, "run_02_scan_all_500.txt")

NPASS = 0
NFAIL = 0
OUT = []


def chk(name, cond, detail=""):
    global NPASS, NFAIL
    if cond:
        NPASS += 1
        OUT.append(f"PASS {name}" + (f"  [{detail}]" if detail else ""))
    else:
        NFAIL += 1
        OUT.append(f"FAIL {name}" + (f"  [{detail}]" if detail else ""))


# ── 讀取 ────────────────────────────────────────────────────
LINE = re.compile(r"^(\d+)\s+(\d+)\s+\[(.*)\]\s*$")
rows = []
with open(SCAN) as f:
    for ln in f:
        m = LINE.match(ln.strip())
        if not m:
            continue
        n = int(m.group(1)); idx = int(m.group(2))
        pairs = re.findall(r"\[\s*(\d+)\s*,\s*(\d+)\s*\]", m.group(3))
        degs = []
        for d, k in pairs:
            degs.extend([int(d)] * int(k))
        rows.append((n, idx, degs))
OUT.append(f"§0 讀入群數 = {len(rows)}（階 ≤ 500；SmallGrp 1.5.4 / GAP 4.14.0）")
chk("0:count==92570", len(rows) == 92570, f"{len(rows)}")
chk("0:burnside", all(sum(d * d for d in degs) == n for n, _, degs in rows))

# ── 特徵表 ──────────────────────────────────────────────────
N = len(rows)
n_arr = np.array([r[0] for r in rows], dtype=np.int64)
DNT_prod = np.ones(N, dtype=np.int64)
DNT_sum = np.zeros(N, dtype=np.int64)
DNT_sumsq = np.zeros(N, dtype=np.int64)
DNT_cnt = np.zeros(N, dtype=np.int64)
NT_prod = np.ones(N, dtype=np.int64)
NT_sum = np.zeros(N, dtype=np.int64)
ALL_prod = np.ones(N, dtype=np.int64)
ALL_sum = np.zeros(N, dtype=np.int64)
ALL_sumsq = np.zeros(N, dtype=np.int64)
DNT_max = np.zeros(N, dtype=np.int64)
IS_A5 = np.zeros(N, dtype=bool)
CAP = 10 ** 7   # 飽和乘積上界：只需分辨 ≤ 500 之 |G| 與 ≤ 30 之常數


def sat_prod(vals):
    p = 1
    for v in vals:
        p *= v
        if p > CAP:
            return CAP
    return p


for j, (n, idx, degs) in enumerate(rows):
    if n == 60 and idx == 5:
        IS_A5[j] = True
    S = sorted({d for d in degs if d > 1})
    if S:
        DNT_prod[j] = sat_prod(S)
        DNT_sum[j] = sum(S)
        DNT_sumsq[j] = sum(d * d for d in S)
        DNT_cnt[j] = len(S)
        DNT_max[j] = S[-1]
    NT = [d for d in degs if d > 1]
    NT_prod[j] = sat_prod(NT)
    NT_sum[j] = sum(NT)
    ALL_prod[j] = sat_prod(degs)
    ALL_sum[j] = sum(degs)
    ALL_sumsq[j] = sum(d * d for d in degs)

OUT.append(f"§0 A5 索引={int(np.where(IS_A5)[0][0])}  |G|={n_arr[IS_A5][0]}  維數={rows[int(np.where(IS_A5)[0][0])][2]}")
chk("0:A5-identified", int(np.where(IS_A5)[0][0]) >= 0 and
    rows[int(np.where(IS_A5)[0][0])][:2] == (60, 5))


def mask_of(bool_arr):
    return int.from_bytes(np.packbits(bool_arr.astype(np.uint8), bitorder="little").tobytes(), "little")


A5MASK = mask_of(IS_A5)
TOTAL = (1 << N) - 1


def hits_of(m):
    cnt = m.bit_count()
    return cnt


def label(m):
    if m == A5MASK:
        return "== {A5} 唯一"
    if m & A5MASK:
        return f"⊇A5, |H|={m.bit_count()}"
    return f"∌A5, |H|={m.bit_count()}"


# ── §A 準則本身 ─────────────────────────────────────────────
OUT.append("§A 準則 C := [∏(相異非平凡) = |G|] ∧ [Σ(相異非平凡) = 12]")
m_prod = mask_of(DNT_prod == n_arr)
m_sum12 = mask_of(DNT_sum == 12)
mC = m_prod & m_sum12
OUT.append(f"  ∏=|G|      : {label(m_prod)}")
OUT.append(f"  Σ=12        : {label(m_sum12)}")
OUT.append(f"  兩者聯合     : {label(mC)}")
chk("A:crit-unique-A5", mC == A5MASK, f"|H|={mC.bit_count()}")

# ── §B 修飾詞敏感度（群層） ─────────────────────────────────
OUT.append("§B 修飾詞敏感度（群層，階 ≤ 500 全集）")
variants = {
    "V0 相異+非平凡（原準則）": (DNT_prod == n_arr) & (DNT_sum == 12),
    "V1 去「相異」：(∏ 非平凡重數積 = |G|) ∧ (Σ 非平凡重數和 = 12)": (NT_prod == n_arr) & (NT_sum == 12),
    "V2 去「非平凡」和條件：Σ(相異全部) = 12": (DNT_sum + 1 == 12),
    "V3 全積（重數）無和條件": (ALL_prod == n_arr),
}
for k, v in variants.items():
    mm = mask_of(v)
    OUT.append(f"  {k}: {label(mm)}")
    if k.startswith("V0"):
        chk("B:V0==A5", mm == A5MASK)
    if k.startswith("V1"):
        chk("B:V1-distinct-load-bearing", mm != A5MASK, f"|H|={mm.bit_count()}")
    if k.startswith("V2"):
        chk("B:V2-nontrivial-load-bearing", mm != A5MASK, f"|H|={mm.bit_count()}")

# ── §C ★ look-elsewhere：常數 c 掃描 ────────────────────────
OUT.append("§C ★ look-elsewhere：C_c := [∏(相異非平凡) = |G|] ∧ [Σ(相異非平凡) = c]，c = 1..40")
rowsC = []
for c in range(1, 41):
    mm = m_prod & mask_of(DNT_sum == c)
    rowsC.append((c, mm))
    OUT.append(f"  c={c:2d}: {label(mm)}")
uniq_c = [c for c, mm in rowsC if mm == A5MASK]
hit_c = [c for c, mm in rowsC if mm != 0]
OUT.append(f"  有命中之 c = {hit_c}")
OUT.append(f"  唯一命中 A₅ 之 c = {uniq_c}")
chk("C:c==12-unique", uniq_c == [12], f"{uniq_c}")
# 另一族：Σ(DNT) = c 單獨之鑑別力
rowsS = []
for c in range(1, 41):
    rowsS.append((c, mask_of(DNT_sum == c).bit_count()))
OUT.append("§C2 單條件 Σ(DNT) = c 之命中群數（無 ∏ 條件）：" + " ".join(f"{c}:{k}" for c, k in rowsS if k))

# ── §D 準則族：單條 / 兩條聯合 ──────────────────────────────
SEL = {
    "DNT": (DNT_prod, DNT_sum, DNT_sumsq, DNT_cnt),
    "NT":  (NT_prod, NT_sum, None, None),
    "DALL": (DNT_prod, DNT_sum + 1, DNT_sumsq + 1, DNT_cnt + 1),
    "ALL": (ALL_prod, ALL_sum, ALL_sumsq, np.ones(N, dtype=np.int64)),
}
FEAT = ["prod", "sum", "sumsq", "count"]
members = []
for sname, feats in SEL.items():
    for fi, fname in enumerate(FEAT):
        arr = feats[fi]
        if arr is None:
            continue
        members.append((f"{sname}.{fname}=", "|G|", mask_of(arr == n_arr)))
        for c in range(1, 31):
            members.append((f"{sname}.{fname}=", str(c), mask_of(arr == c)))
OUT.append(f"§D 準則族大小 |F| = {len(members)}（selector×feature×target）")
uniq_single = [m for m in members if m[2] == A5MASK]
OUT.append(f"  單條即唯一命中 A₅ 者：{len(uniq_single)} 條")
for a, b, _ in uniq_single:
    OUT.append(f"    {a}{b}")
sing_hit = [m for m in members if m[2] & A5MASK]
OUT.append(f"  含 A₅ 之單條：{len(sing_hit)} 條")

uniq_pair = []
npair = 0
for i in range(len(members)):
    mi = members[i][2]
    if not (mi & A5MASK):
        continue
    for j in range(i + 1, len(members)):
        mj = members[j][2]
        if not (mj & A5MASK):
            continue
        npair += 1
        inter = mi & mj
        if inter == A5MASK:
            uniq_pair.append((members[i][0] + members[i][1], members[j][0] + members[j][1]))
OUT.append(f"§D 兩條聯合：可測組合 {npair} 對；其中恰好唯一命中 A₅ 者 {len(uniq_pair)} 對")
for a, b in uniq_pair[:40]:
    OUT.append(f"    [{a}] ∧ [{b}]")
if len(uniq_pair) > 40:
    OUT.append(f"    …（其餘 {len(uniq_pair)-40} 對見 data/run_03_pairs.txt）")
with open(os.path.join(DATA, "run_03_pairs.txt"), "w") as f:
    for a, b in uniq_pair:
        f.write(f"[{a}] AND [{b}]\n")
LOOK_ELSEWHERE = len(uniq_pair) / npair if npair else float("nan")
OUT.append(f"★ look-elsewhere 率 = {len(uniq_pair)}/{npair} = {LOOK_ELSEWHERE:.4%}")
# 本 O 之準則屬其中一對
CANON = ("DNT.prod=|G|", "DNT.sum=12")
IS_IN = CANON in [(a, b) for a, b in uniq_pair] or CANON in [(b, a) for a, b in uniq_pair]
chk("D:canon-in-pair-family", IS_IN, f"{CANON}")

# ── 落檔 ───────────────────────────────────────────────────
print("\n".join(OUT))
print("─" * 60)
print(f"TOTAL PASS={NPASS} FAIL={NFAIL}")
sys.exit(0 if NFAIL == 0 else 1)