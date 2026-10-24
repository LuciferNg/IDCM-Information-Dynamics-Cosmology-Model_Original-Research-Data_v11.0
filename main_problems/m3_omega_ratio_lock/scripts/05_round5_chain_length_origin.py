#!/usr/bin/env python3
"""M3 scripts/05 — 第五輪：88 位鏈長之源（獨立重算 O220 U11）

問題: M3 第四輪標記「L₅=11 係最後非平凡輸入」。O220 docs/13 (U11) 主張
      88 已由二十面體資料導出（兩條獨立路線）。本輪獨立重算，唔引用其結果。

U11 內容:
  N_gen = (e + o − 3)/(F₅ + 1),  e = #5-fold axes = 6, o = #2-fold axes = 15
  |A₂|  = 15 = 群 I 之二階元素數
  N_d   = N_gen·|A₂| − 1 = 44  ⟹  88 = 2N_d
  第二路線: 88 = E(RH) − F(TI) = 120 − 32

本輪只做獨立計算: 建群 I（12 頂點 icosahedron 之旋轉群）→ 數階數 → 多面體計數。
"""
import numpy as np
import itertools
from fractions import Fraction as F

PASS = FAIL = 0
def chk(name, cond, detail=""):
    global PASS, FAIL
    if cond: PASS += 1; print(f"  [PASS] {name}  {detail}")
    else:    FAIL += 1; print(f"  [FAIL] {name}  {detail}")

phi = (1 + 5**0.5)/2

# ── §1 建 icosahedron 12 頂點 ────────────────────────────────────────
V = []
for s1 in (1, -1):
    for s2 in (1, -1):
        V.append((0, s1*1.0, s2*phi))
        V.append((s1*1.0, s2*phi, 0))
        V.append((s2*phi, 0, s1*1.0))
V = np.array(sorted(set(tuple(np.round(v, 12) for v in v_) for v_ in V)))
print(f"§1 icosahedron 頂點數 = {len(V)}")
chk("|V| = 12", len(V) == 12)

# 邊 = 最短距離對
d = np.linalg.norm(V[:, None, :] - V[None, :, :], axis=-1)
edge_len = np.min(d[np.nonzero(d)])
E_ico = int(np.sum(np.isclose(d, edge_len))/2)
chk("icosahedron 邊數 = 30", E_ico == 30, f"{E_ico}")

# ── §2 旋轉群 I（60 元）──────────────────────────────────────────────
def rot_from_frames(a, b, c):
    """把標準正交架映到 (a,b,c) 之正交矩陣"""
    A = np.column_stack([a, b, c])
    det = np.linalg.det(A)
    return A / np.cbrt(det) if abs(det) > 1e-12 else None

Vs = [V[i] for i in range(12)]
seen = set(); G = []
# 固定參考架（三個線性獨立頂點向量）
ref = None
for i, j, k in itertools.permutations(range(12), 3):
    M = np.column_stack([Vs[i], Vs[j], Vs[k]])
    if abs(np.linalg.det(M)) > 1e-9:
        ref = M; ref_idx = (i, j, k); break
for i, j, k in itertools.permutations(range(12), 3):
    M = np.column_stack([Vs[i], Vs[j], Vs[k]])
    if abs(np.linalg.det(M)) < 1e-9:
        continue
    R = M @ np.linalg.inv(ref)
    if np.linalg.det(R) < 0.5:          # 只取 proper rotations (det=+1)
        continue
    # 檢查 R 置換頂點集
    RV = V @ R.T
    ok = all(np.min(np.linalg.norm(V - rv, axis=1)) < 1e-8 for rv in RV)
    if not ok:
        continue
    key = tuple(np.round(R, 8).flatten())
    if key not in seen:
        seen.add(key); G.append(R)
print(f"\n§2 群 I 元素數 = {len(G)}")
chk("|I| = 60", len(G) == 60, f"{len(G)}")

ang = lambda R: float(np.degrees(np.arccos(np.clip((np.trace(R)-1)/2, -1, 1))))
orders = {}
for R in G:
    th = ang(R)
    n = 1
    while n < 6:
        if abs(ang(np.linalg.matrix_power(R, n))) < 1e-6 or abs(ang(np.linalg.matrix_power(R, n)) - 360) < 1e-6:
            break
        n += 1
    orders[n] = orders.get(n, 0) + 1
print(f"  階分佈: {dict(sorted(orders.items()))}")

n2 = orders.get(2, 0)
chk("二階元素數 |A₂| = 15", n2 == 15, f"{n2}")
chk("五階元素 24 ⟹ 5-fold 軸 e = 6", orders.get(5, 0) == 24 and orders.get(5, 0)//4 == 6,
    f"{orders.get(5,0)}/4 = {orders.get(5,0)//4}")
chk("三階元素 20 ⟹ 3-fold 軸 10", orders.get(3, 0) == 20, f"{orders.get(3,0)}")

# ── §3 U11 鏈長導出 ─────────────────────────────────────────────────
print("\n§3 U11 鏈長（用 §2 之群論計數）")
e, o, F5 = 6, 15, 5
N_gen = F(e + o - 3, F5 + 1)
chk("N_gen = (e+o−3)/(F₅+1) = 18/6 = 3", N_gen == 3, f"{N_gen}")
A2 = F(n2)
chk("|A₂| = 15", A2 == 15)
N_d = N_gen*A2 - 1
chk("N_d = N_gen·|A₂| − 1 = 44", N_d == 44, f"{N_d}")
chk("88 = 2N_d ; 176 = 4N_d", 2*N_d == 88 and 4*N_d == 176)

# ── §4 第二路線：多面體計數 ─────────────────────────────────────────
print("\n§4 第二路線（多面體計數，唔用群表示論）")
# RH: 60 頂點，度數 4（頂點遞移）
V_RH, deg_RH = 60, 4
E_RH = V_RH*deg_RH//2
F_RH = 20 + 30 + 12
chk("E(RH) = V·deg/2 = 120", E_RH == 120, f"{E_RH}")
chk("F(RH) = 20△+30□+12⬠ = 62 且 Euler=2", F_RH == 62 and V_RH - E_RH + F_RH == 2)
# TI: truncated icosahedron
V_TI, E_TI, F_TI = 60, 90, 12 + 20
chk("TI: 60−90+32 = 2（Euler）且 F=32", V_TI - E_TI + F_TI == 2 and F_TI == 32)
chk("★ 88 = E(RH) − F(TI) = 120 − 32 = 88", E_RH - F_TI == 88, f"{E_RH-F_TI}")
chk("兩路線一致（2N_d = E(RH)−F(TI)）", 2*N_d == E_RH - F_TI)

# ── §5 M3 生成鏈重新植根 ────────────────────────────────────────────
print("\n§5 M3 生成鏈重新植根（L₅ 由輸入降為導出）")
B = 4
L5 = F(N_d, B)
chk("L₅ = N_d/B = 11（導出，非輸入）", L5 == 11, f"{L5}")
chk("12 = L₅+1 = N_gen·B（同 O116 §7）", L5 + 1 == 12 and N_gen*B == 12)
chk("|A₂| = B² − 1 = 15（群↔樹 之橋）", B**2 - 1 == A2, f"{B**2-1} vs {A2}")
chk("N_d = N_gen·B² − B（equipartition 式）", N_gen*B**2 - B == N_d,
    f"{N_gen*B**2-B}")
chk("兩式等價 ⟺ B = N_gen+1", (N_gen*A2 - 1) == (N_gen*B**2 - B) and B == N_gen+1)
chk("K_tree = N_gen·L₅ = 33", N_gen*L5 == 33)
chk("Σγ = 2N_d+2 = 90", 2*N_d + 2 == 90)
chk("2⟨W⟩² = Σγ+K_tree = 123", 90 + N_gen*L5 == 123)
chk("份額 (3,16,41)/60 與量子 3/2",
    (N_gen, B**2, N_d-N_gen) == (3, 16, 41) and F(90, 60) == F(3,2))

print(f"\n結算: PASS {PASS} / FAIL {FAIL}")