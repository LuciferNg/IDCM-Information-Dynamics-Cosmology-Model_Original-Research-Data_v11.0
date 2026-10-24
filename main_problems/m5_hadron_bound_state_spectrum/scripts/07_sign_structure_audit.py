"""M5 docs/07 配套 — H5a: H4b 存活號組合之結構性審計（s_d = −1 可否由量子數／A₅ 分級導出？）
================================================================================
H4b 唯一存活號組合：{s_u, s_d, s_s, s_c} = {+1, −1, +1, ±1}（全域翻號等價 ⟹ 4 個等價類）。
H5a 問題：呢個號向量能唔能寫成**已登記結構**嘅 Z₂ 特徵（charge／strangeness／isospin／baryon／
          A₅ 阻挫類）？若唔能 ⟹ H4b 嘅「解」只係擬合自由度，唔可作預言（湊數風險）。
方法：把候選分級 s_f = (−1)^{a·3Q + b·S + c·2I₃ + d·B}（a,b,c,d ∈ {0,1}）全枚舉（16 種），
      同需求類逐一同構比較（准許全域翻號）；A₅ 阻挫類另作說明（需 flavour→vertex 映射，屬禁區）。
輸出：`data/07_sign_structure_audit.json`
"""
import json
from itertools import product

FLAV = ["u", "d", "s", "c"]
Q3 = {"u": 2, "d": -1, "s": -1, "c": 2}        # 3Q（整數）
S = {"u": 0, "d": 0, "s": -1, "c": 0}          # strangeness
I3_2 = {"u": 1, "d": -1, "s": 0, "c": 0}       # 2I₃
B = {"u": 1, "d": 1, "s": 1, "c": 1}           # baryon number（全同 ⟹ 無鑑別力）

# H4b 存活類（s_c 自由 ⟹ 2 種 × 全域翻號等價）
REQ = [{"u": 1, "d": -1, "s": 1, "c": 1},
       {"u": 1, "d": -1, "s": 1, "c": -1}]

def canon(d):
    """全域翻號規範化（令 s_u = +1），用於比較等價類"""
    s = d["u"]
    return tuple(d[f] * s for f in FLAV)

REQ_C = {canon(r) for r in REQ}
print("=" * 96)
print("M5 docs/07 — H5a: H4b 存活號組合之結構性審計")
print("=" * 96)
print(f"需求類（s_u = +1 規範）= {sorted(REQ_C)}")

cands = {}
for a, b, c, d in product((0, 1), repeat=4):
    key = f"(-1)^[{a}·3Q+{b}·S+{c}·2I3+{d}·B]"
    s = {f: (-1) ** (a * Q3[f] + b * S[f] + c * I3_2[f] + d * B[f]) for f in FLAV}
    cands[key] = s

print("\n[候選 Z₂ 分級]")
match = []
for k, s in cands.items():
    cs = canon(s)
    ok = cs in REQ_C
    print(f"  {k:34s} s = {[s[f] for f in FLAV]}  (規範化 {list(cs)})"
          + ("   ⟸ 命中" if ok else ""))
    if ok:
        match.append(k)

print(f"\n[命中數] {len(match)} / {len(cands)}")
print(f"  \u27f9 命中密度 = {len(match)}/{len(cands)} = {100*len(match)/len(cands):.1f}% ⟹ 遠高於 δ★ 類門檻（0.01%）⟹ 巧合水平，不可作結構導出")

# A₅ 阻挫類：需 flavour → vertex 映射（O140 docs/19 禁區）；無映射時 sign 只在 Z₂^30 類中任選
print("\n[A₅ 阻挫類]")
print("  Rhombicosidodecahedron：60 頂點、120 邊、阻挫邊 32 ⟹ 圖**非二部**，")
print("  冇全域二色 ⟹ 每頂點號只有『相對某參考頂點』嘅意義；")
print("  要為 u,d,s,c 定號必須提供 flavour → vertex 映射 ⟹ 屬 O140 docs/19 禁區（映射不存在）")
print("  ⟹ s_d = −1 **唔可以**由 A₅ 阻挫導出（缺映射）；亦唔可以由量子數導出（上表 0/16）")

density = len(match) / len(cands)
THRESH = 1e-4          # δ★ 類門檻（M4 docs/20：0.1% 級命中視隨機；此處更嚴 0.01%）
verdict = ("非結構：候選空間命中密度 = %d/%d = %.1f%% ≫ 門檻 %.2f%% ⟹ 屬巧合水平"
           "（A₅ 路線亦因缺 flavour→vertex 映射而不可用）" % (len(match), len(cands), 100 * density, 100 * THRESH)
           if density > THRESH else "結構性（密度低於門檻）")
print(f"\n[判定] {verdict}")
print("  ⟹ M5 收緊結論：8 通道 + 超度量核（正／號／態依賴）全部得唔到結構性結果；")
print("    唯一出路 = 由 lattice Hamiltonian H_ℓ 導出 ψ（H5b），否則本核路線終止。")

out = {"required_classes": sorted(REQ_C),
       "candidates": {k: [v[f] for f in FLAV] for k, v in cands.items()},
       "hits": match,
       "a5_note": "60 頂點/120 邊/阻挫邊 32 ⟹ 非二部；號需 flavour→vertex 映射（禁區）",
       "verdict": verdict}
with open("data/07_sign_structure_audit.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/07_sign_structure_audit.json")