#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""S4：E₈ 分塊外存之「機制驗證」＋ 規模設計數字
   只驗證「分塊累加 == 稠密計算」（float64 bit-identical），不涉及 E₈ 算子語義
"""
import json
import os

import numpy as np

rng = np.random.default_rng(20261127)
N_IN, N_OUT, N_CH = 1 << 12, 1 << 11, 17        # 小規模：輸入 4096、輸出 2048、17 通道

# 稠密參考：out = Σ_ch M_ch @ in
In = rng.standard_normal(N_IN)
Ms = [rng.standard_normal((N_OUT, N_IN)) for _ in range(N_CH)]
dense = np.zeros(N_OUT)
for M in Ms:
    dense += M @ In

# 分塊：輸出分 8 塊，累加器為 memmap（模擬 29.43 GB 級外存累加器）
p = "/tmp/s4_acc.bin"
blk = N_OUT // 8
acc = np.memmap(p, dtype=np.float64, mode="w+", shape=(N_OUT,))
for b in range(8):
    lo, hi = b * blk, (b + 1) * blk
    part = np.zeros(hi - lo)
    for M in Ms:                                  # 通道切片（模擬 rank 表重取）
        part += M[lo:hi] @ In
    acc[lo:hi] = part
    acc.flush()
chunked = np.array(acc)

bit_equal = bool(np.array_equal(chunked.view(np.uint64), dense.view(np.uint64)))
print("分塊 vs 稠密：bit-identical =", bit_equal, "｜max|diff| =", float(np.abs(chunked - dense).max()))

# ---- E₈ 規模設計數字（精確整數）----
STATES = 3679075400
out = {
    "mechanism_check": {"n_out": N_OUT, "n_channels": N_CH, "blocks": 8,
                        "bit_identical": bit_equal,
                        "max_abs_diff": float(np.abs(chunked - dense).max()),
                        "note": "僅驗證累加機制（float64 結合律無損），不驗 E₈ 語義"},
    "e8_sizing": {
        "states_V7_9": STATES,
        "bytes_float32": STATES * 4,
        "bytes_float64": STATES * 8,
        "gib_float32": STATES * 4 / 2**30,
        "gib_float64": STATES * 8 / 2**30,
        "chunk_states_2p24": 1 << 24,
        "chunk_gib_float64": (1 << 24) * 8 / 2**30,
        "n_chunks": -(-STATES // (1 << 24)),
        "disk_free_gib_ref": 119,
        "vram_ref_gib": 1.25,
    },
}
with open("/home/wsl/idcm_ckpt/e8/mrB/s4_chunk_validate.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
os.remove(p)
e = out["e8_sizing"]
print(f"E₈ 規模：{STATES} 態｜float32 {e['gib_float32']:.3f} GiB｜float64 {e['gib_float64']:.3f} GiB")
print(f"分塊 2^24 態 = {e['chunk_gib_float64']:.3f} GiB/塊 ⟹ {e['n_chunks']} 塊（外存累加器 {e['gib_float64']:.2f} GiB ≪ 119 GiB）")
