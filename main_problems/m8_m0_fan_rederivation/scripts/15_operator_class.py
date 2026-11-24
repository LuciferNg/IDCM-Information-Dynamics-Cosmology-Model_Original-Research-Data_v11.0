#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
M8 scripts/15 — 門 (C) 算子類判定（**判定輪，唔係推導輪**）
================================================================================
問題：門 (C) 需要「粗粒化算子 C ∉ Comm(K)，具唯一 RG 特徵值 r ≠ 1，N 次到 IR」。
本輪只判三件事：

  B2（定理，電腦核）  **任何** 抽樣式粗粒化（surjective linear decimation ＋ 其 pseudo-inverse）
                     都係 projection ⟹ 譜 ⊆ {0,1} ⟹ **冇 r ≠ 1**。
                     對 block size 2/4/8 與均勻/非均勻權重逐一驗。
  B3（判定）          `[P,K] ≠ 0` 可否同時成立（唔對易 ≠ 有標度）—— 實測 Frobenius 範數。
  B1（定理，電腦示範）「存在非交換 C 具唯一 r」對**任何** r 都可滿足 ⟹ 該公理對目標零約束
                     ⟹ 零預言內容（`dof` = 1 連續 + 1 整數）。
                     ⚠ 此項直接判死「新 C」版本，唔需要動 v10 O210。

執行: /usr/bin/python3 scripts/15_operator_class.py
輸出: data/15_operator_class.json、data/run_15.log
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

LOG = []
RES = {}
NFAIL = 0


def chk(tag, cond, detail=""):
    global NFAIL
    ok = bool(cond)
    if not ok:
        NFAIL += 1
    line = f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  |  {detail}" if detail else "")
    print(line)
    LOG.append(line)
    RES[tag] = {"pass": ok, "detail": detail}
    return ok


def hdr(s):
    print("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)
    LOG.append("\n" + "=" * 78 + f"\n{s}\n" + "=" * 78)


def chain_K(N=88, periodic=False):
    """88 站鏈 Laplacian（2 對角、−1 鄰接；open 或 periodic）"""
    K = 2.0 * np.eye(N)
    for i in range(N - 1):
        K[i, i + 1] = K[i + 1, i] = -1.0
    if periodic:
        K[0, N - 1] = K[N - 1, 0] = -1.0
    return K


def decimation(N=88, block=2, weights=None):
    """A: N → N/block（區塊平均）；P = Aᵀ(A Aᵀ)⁻¹A（pseudo-inverse 重建 ⟹ projection）"""
    m = N // block
    A = np.zeros((m, N))
    for b in range(m):
        w = np.ones(block) if weights is None else np.asarray(weights, dtype=float)
        A[b, b * block:(b + 1) * block] = w / w.sum()
    G = A @ A.T
    P = A.T @ np.linalg.inv(G) @ A
    return A, P


def main():
    hdr("§0 定義（唔涉任何目標值；純結構）")
    print("  K  = 88 站鏈 Laplacian（2 對角、−1 鄰接）")
    print("  A  = 區塊平均抽樣 N→N/b；P = Aᵀ(A Aᵀ)⁻¹A（廣義重建 ⟹ 投影）")
    print("  C  = 門 (C) 所需之一般算子（∉ Comm(K)，唯一特徵值 r ≠ 1）")

    hdr("§1 B2：抽樣式粗粒化之譜 ⊆ {0,1}（塊 2/4/8 × 均勻/非均勻）")
    cases = [(2, None), (2, [1, 3]), (4, None), (4, [1, 2, 3, 5]), (8, None), (8, [1, 1, 2, 2, 3, 3, 5, 7])]
    for block, w in cases:
        A, P = decimation(88, block, w)
        ev = np.linalg.eigvalsh((P + P.T) / 2)
        n1 = int(np.sum(np.abs(ev - 1) < 1e-9))
        n0 = int(np.sum(np.abs(ev) < 1e-9))
        tag = f"B2.块{block}{'·均' if w is None else '·非均'}"
        chk(f"{tag}：譜 ⊆ {{0,1}}（1×{n1}, 0×{n0}，其餘 {88-n1-n0}）",
            n0 + n1 == 88, f"1×{n1}, 0×{n0}, 最大|λ|={np.abs(ev).max():.3g}")
        RES.setdefault("_spec", {})[tag] = {"ones": n1, "zeros": n0}

    hdr("§2 B3：`[P,K] ≠ 0` 與「有標度」無關（唔對易 ≠ r ≠ 1）")
    K = chain_K(88)
    A, P = decimation(88, 2, None)
    comm = P @ K - K @ P
    fro = float(np.linalg.norm(comm))
    evP = np.linalg.eigvalsh((P + P.T) / 2)
    print(f"  ‖[P,K]‖_F = {fro:.6f}")
    print(f"  spec(P) 唯一非平凡值 = {sorted(set(np.round(evP, 9)))}")
    chk("B3.1 自然 decimation 之 [P,K] ≠ 0", fro > 1e-6, f"‖·‖_F = {fro:.6f}")
    chk("B3.2 同時其譜只有 {0,1} ⟹ 非交換 ≠ 有標度",
        bool(np.all(np.isclose(np.abs(evP), 0) | np.isclose(np.abs(evP), 1))),
        f"max|λ| = {np.abs(evP).max():.6f}")
    RES["_comm"] = {"fro": fro, "spec_uniq": sorted(set(np.round(evP, 9).tolist()))}

    hdr("§3 B1：示範「∃ 非交換 C 具唯一 r」對任何 r 皆可滿足 ⟹ 零約束")
    rng = np.random.default_rng(20261124)
    targets = ["1e5", "3.7", "2.2852930605e-21"]      # 含 μ 本身
    for ts in targets:
        r = float(ts)
        C0 = np.diag([r] + [1.0] * 87)                 # 唯一 r ≠ 1，其餘 1
        S = np.eye(88) + 0.05 * rng.standard_normal((88, 88))
        C = S @ C0 @ np.linalg.inv(S)
        ev = np.linalg.eigvals(C)
        uniq_r = ev[np.abs(ev - 1) > 1e-9]
        comm = C @ K - K @ C
        fro = float(np.linalg.norm(comm))
        print(f"  r = {ts:>18}：非 1 特徵值數 = {len(uniq_r)}（值 {np.round(uniq_r,12)}）｜‖[C,K]‖_F = {fro:.3g}")
        chk(f"B1.{ts} 構造出唯一 r 且 C ∉ Comm(K)", len(uniq_r) == 1 and fro > 1e-3,
            f"非1值 {len(uniq_r)}，‖[C,K]‖_F={fro:.3g}")
        RES.setdefault("_demo", {})[ts] = {"n_nonunit": int(len(uniq_r)), "comm_fro": fro}

    hdr("§4 判詞")
    print("  B2：抽樣式 decimation ⟹ 投影 ⟹ 譜 {0,1} ⟹ **冇 r ≠ 1**（全部塊大小/權重）")
    print("  B1：一般 C 可對任何 r 構造 ⟹ 該公理對 μ 零約束 ⟹ 零預言 ⟹ 不可證偽")
    print("  ⟹ 門 (C) 兩讀法皆死；且 **唔需要** 判定是否撞 v10 O210（死因更早）")
    chk("V4.1 判詞成立（B1 ∧ B2 皆 PASS）", NFAIL == 0, f"FAIL = {NFAIL}")

    hdr("§5 匯總")
    nchk = len([k for k in RES if not k.startswith("_")])
    print(f"檢查項 = {nchk}，FAIL = {NFAIL}")
    RES["_summary"] = {"checks": nchk, "fail": NFAIL}
    (DATA / "15_operator_class.json").write_text(
        json.dumps(RES, ensure_ascii=False, indent=2), encoding="utf-8")
    (DATA / "run_15.log").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\n寫出: {DATA / '15_operator_class.json'}\n寫出: {DATA / 'run_15.log'}")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
