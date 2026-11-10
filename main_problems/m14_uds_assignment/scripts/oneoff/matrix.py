import json, os
import numpy as np

CK = "/home/wsl/.hermes/cache/scratch/e7ckpt2"
CB = "/home/wsl/.hermes/cache/scratch/e7ckptB"
GO = "/home/wsl/.hermes/cache/scratch/gpuout"
OUT = "/home/wsl/Projects/IDCM-Original/IDCM-v11.0--Reality-Architecture/main_problems/m14_uds_assignment/data/chain_matrix.json"


def load(p, keys=("E2", "E3", "E4", "E5", "E6", "E7")):
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    return {k: d.get(k) for k in keys}


A_cpu = json.load(open(f"{CK}/chain_result_f64.json"))
B_cpu = json.load(open(f"{CB}/chain_resultB.json")) if os.path.exists(f"{CB}/chain_resultB.json") else None
A_gpu = json.load(open(f"{GO}/gpu_chain_A.json"))
B_gpu = json.load(open(f"{GO}/gpu_chain_B.json")) if os.path.exists(f"{GO}/gpu_chain_B.json") else None

keys = ["E1", "E2", "E3", "E4", "E5", "E6", "E7"]
tab = {}
for name, src in (("A_cpu", A_cpu), ("A_gpu", A_gpu), ("B_cpu", B_cpu), ("B_gpu", B_gpu)):
    if src is None:
        continue
    tab[name] = {k: (float(src[k]) if src.get(k) is not None else None) for k in keys}
    tab[name]["E1"] = 0.004871899480351504   # E_1 ≡ γ_n（定義值，兩路徑相同）

rows = []
for k in keys:
    r = {"E": k}
    for name in ("A_cpu", "A_gpu", "B_cpu", "B_gpu"):
        r[name] = tab.get(name, {}).get(k)
    for pair in (("A_cpu", "A_gpu"), ("B_cpu", "B_gpu")):
        a, b = r.get(pair[0]), r.get(pair[1])
        if a is not None and b is not None and a != 0:
            r[f"rel_{pair[0][0]}"] = abs(b - a) / abs(a)
        else:
            r[f"rel_{pair[0][0]}"] = None
    rows.append(r)

out = {"rows": rows, "sources": {k: v for k, v in tab.items()}}
json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)

print(f"{'E':>4} {'A_CPU':>20} {'A_GPU':>20} {'rel(A)':>10}   {'B_CPU':>20} {'B_GPU':>20} {'rel(B)':>10}")
for r in rows:
    def f(x):
        return "—" if x is None else (f"{x:+.12f}" if abs(x) < 1e3 else f"{x:+.9e}")
    def g(x):
        return "—" if x is None else f"{x:.2e}"
    print(f"{r['E']:>4} {f(r['A_cpu']):>20} {f(r['A_gpu']):>20} {g(r['rel_A']):>10}   "
          f"{f(r['B_cpu']):>20} {f(r['B_gpu']):>20} {g(r['rel_B']):>10}")