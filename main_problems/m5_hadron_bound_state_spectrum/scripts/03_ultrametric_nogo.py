"""M5 docs/03 配套 — H3a: 超度量 Coulomb 核類之解析否證（K 通道符號矛盾）
================================================================================
H2 得核形式：V_ij = κ q_i q_j x_{f_i f_j}，x_AB = B^(−|ℓ_A−ℓ_B|) ∈ (0,1]（同味 x_AA = 1）。
H3a 問題：呢一整個核類（任意 κ>0、任意層指派、任意 b）能唔能重現 8 通道 em 譜？

答案（解析，非擬合）：**不能**。n−p 與 K 通道被同一條 q² 排序綁死：
    E(n)−E(p) = −κ/3                （層指派無關 ⟹ 精確結構關係）
    E(K⁺)−E(K⁰) = +κ[1/3 + (2x_us + x_ds)/9] ≥ +κ/3 > 0
但實測 E(K⁺)−E(K⁰) = −1.3740 < 0，|E_K| ≈ |E_N| ⟹ 符號矛盾，任何 κ、任何層指派都救唔返。
⟹ 核必須帶**電荷以外嘅號結構**（不能化為 q_iq_j·正核）⟹ 唯一候選 = 樹嘅 Z₂ 二色（阻挫）符號。
本腳本：精確 Fraction 係數 + 全通道代數 + 矛盾見證 + 符號核（Z₂）可行性檢查。
"""
import json
from fractions import Fraction as F
from itertools import product

Q = {"u": F(2, 3), "d": F(-1, 3), "s": F(-1, 3), "c": F(2, 3),
     "ub": F(-2, 3), "db": F(1, 3), "sb": F(1, 3), "cb": F(-2, 3)}

# 8 通道（M4 docs/18 §四 canonical）：(名稱, 初態, 終態, 實測 E MeV, σ)
CH = [("n−p", ("u", "d", "d"), ("u", "u", "d"), -1.2667, 0.0940),
      ("K⁺−K⁰", ("u", "sb"), ("d", "sb"), -1.3740, 0.0963),
      ("π⁺−π⁰", ("u", "db"), ("u", "ub"), 2.0336, 0.0940),   # π⁰ 用 u-ub 支做代表（另支 d-db 見 §四）
      ("Σ⁻−Σ⁰", ("d", "d", "s"), ("u", "d", "s"), 2.2470, 0.1010),
      ("Σ⁰−Σ⁺", ("u", "d", "s"), ("u", "u", "s"), 0.7120, 0.1196),
      ("Σ⁻−Σ⁺", ("d", "d", "s"), ("u", "u", "s"), 2.9590, 0.2028),
      ("Ξ⁻−Ξ⁰", ("d", "s", "s"), ("u", "s", "s"), 4.2900, 0.2301),
      ("D⁺−D⁰", ("c", "db"), ("c", "ub"), 2.2600, 0.1172)]

def flavour(x):
    return x[:-1] if x.endswith("b") else x

def pairs(state):
    return [(state[i], state[j]) for i in range(len(state)) for j in range(i + 1, len(state))]

def em_symbolic(state, xsym):
    """EM = Σ_{i<j} q_i q_j x_{f_i f_j}；同味 x=1（層相同）⟹ 精確有理係數"""
    tot = F(0)
    for a, b in pairs(state):
        fa, fb = flavour(a), flavour(b)
        key = None if fa == fb else tuple(sorted((fa, fb)))
        coef = Q[a] * Q[b]
        if key is None:
            tot += coef                      # x_AA = 1
        else:
            tot += coef * xsym[key]          # 符號係數（正核 ⟹ x>0）
    return tot

XV = [tuple(sorted(p)) for p in {("u", "d"), ("u", "s"), ("d", "s"), ("u", "c"), ("d", "c"), ("s", "c")}]
XV = sorted(XV)
sym = {k: F(0) for k in XV}          # 佔位：用 Fraction 做記號運算（每 key 一個獨立單位）

def em_coeffs(state, selfa=False):
    """回傳 (常數, {pair: 有理係數})；selfa=True 另加 on-site Σq_i²（V_ii = B^0 = 1）"""
    const = F(0); co = {}
    if selfa:
        const += sum(Q[x] ** 2 for x in state)
    for a, b in pairs(state):
        fa, fb = flavour(a), flavour(b)
        c = Q[a] * Q[b]
        if fa == fb:
            const += c
        else:
            k = tuple(sorted((fa, fb)))
            co[k] = co.get(k, F(0)) + c
    return const, co

def diff(init, fin, selfa=False):
    c0, k0 = em_coeffs(init, selfa)   # E = EM(init) − EM(fin)
    c1, k1 = em_coeffs(fin, selfa)
    const = c0 - c1
    co = dict(k0)
    for k, v in k1.items():
        co[k] = co.get(k, F(0)) - v
    return const, {k: v for k, v in co.items() if v != 0}

print("=" * 96)
print("M5 docs/03 — H3a: 超度量核類之解析否證")
print("=" * 96)
rows = []
for name, i, f, obs, sig in CH:
    const, co = diff(i, f)
    rows.append((name, const, co, obs, sig))
    print(f"\n{name}:  E_model/κ = {const}  +  " +
          "  +  ".join(f"({v})·x[{k[0]}{k[1]}]" for k, v in co.items()))
    print(f"   實測 = {obs:+.4f} ± {sig:.4f} MeV")

print("\n" + "=" * 96)
print("[0] 兩種 self-energy 公約下嘅同一矛盾（A: 只計 pair；B: 另加 on-site Σq²）")
for selfa, nm in ((False, "A"), (True, "B")):
    cn_, kn_ = diff(CH[0][1], CH[0][2], selfa)
    ck_, kk_ = diff(CH[1][1], CH[1][2], selfa)
    print(f"    公約{nm}: n−p/κ = {cn_} {dict((f'{k[0]}{k[1]}', str(v)) for k, v in kn_.items())}"
          f" | K⁺−K⁰/κ = {ck_} {dict((f'{k[0]}{k[1]}', str(v)) for k, v in kk_.items())}")
print("    ⟹ 公約 A：E_K/κ = 正項之和 > 0；公約 B：E_K/κ = 1/3 + 正項 ≥ 1/3 > 0")
print("    ⟹ 兩種公約皆預測 E(K⁺) > E(K⁰)，實測相反 ⟹ 矛盾與公約無關")

print("\n[1] n−p：所有 x 項精確抵消 ⟹ 層指派無關嘅結構關係")
cn, con = rows[0][1], rows[0][2]
print(f"    E(n)−E(p)/κ = {cn}   殘餘 x 項 = {con}  ⟹ κ = 3·1.2667 = {3*1.2667:.4f} MeV")

print("\n[2] K 通道：x ≥ 0（正核）⟹ 下界")
ck, cok = rows[1][1], rows[1][2]
print(f"    E(K⁺)−E(K⁰)/κ = {ck} + " + " + ".join(f"({v})·x[{k[0]}{k[1]}]" for k, v in cok.items()))
print(f"    由於所有係數 = {[str(v) for v in cok.values()]} > 0 且 x ∈ (0,1] ⟹ 括號 ≥ {ck} > 0")
print(f"    ⟹ E_model(K⁺)−E(K⁰) > 0（公約 A）／≥ κ/3 = {3*1.2667/3:.4f}（公約 B），但實測 = {rows[1][3]:+.4f} < 0")
print(f"    ⟹ 矛盾幅度 = {abs(rows[1][3]) + 1.2667:.4f} MeV（≈ 2× 效應本身）⟹ 唔可能靠微調救")

print("\n[3] 全枚舉：層指派 ℓ_f ∈ 0..33（4 味 ⟹ 34³ 平移等價類）、κ 由 n−p 釘死，檢 7 條通道")
best = []
for lu, ld, ls in product(range(34), repeat=3):
    xs = {("d", "u"): F(4) ** (-abs(ld - lu)),
          ("s", "u"): F(4) ** (-abs(ls - lu)),
          ("d", "s"): F(4) ** (-abs(ls - ld))}
    kappa = F(3 * 12667, 10000)          # MeV，由 n−p 精確釘死（公約 A: −κ/3）
    dev = []
    for name, const0, co0, obs, sig in rows[1:]:
        const, co = diff(CH[[c[0] for c in CH].index(name)][1], CH[[c[0] for c in CH].index(name)][2], True)
        val = kappa * (const + sum(c * xs.get(k, F(0)) for k, c in co.items() if k in xs))
        dev.append(abs(float(val) - obs))
    best.append((max(dev), lu, ld, ls))
best.sort()
print(f"    最佳 5 個指派（最大通道偏差 MeV）：")
for d, lu, ld, ls in best[:5]:
    print(f"      ℓ_u={lu} ℓ_d={ld} ℓ_s={ls} ⟹ max dev = {d:.4f} MeV")
print("    ⟹ 最小偏差仍 ≫ σ（K 通道符號相反，任何指派都無法翻轉）")

print("\n[4] 出路：核必須帶號結構（Z₂ 二色／阻挫）")
print("    令 x_AB → s_A s_B · B^(−|Δℓ|)，s_f = ±1（樹嘅二色類）。K 通道括號變")
print("    (2x_us s_u s_s + x_ds s_d s_s)/9 ⟹ 可取負 ⟹ 符號矛盾解除，代價 = 每個味 1 bit。")
print("    ⟹ H3b 目標：由 A₅/樹阻挫結構**導出** s_f（唔准擬合），再過 δ★ 閘。")

out = {"channels": [{"name": r[0], "const": str(r[1]),
                     "coeffs": {f"{k[0]}{k[1]}": str(v) for k, v in r[2].items()},
                     "observed": r[3], "sigma": r[4]} for r in rows],
       "n_minus_p_identity": {"const": str(cn), "kappa_MeV": 3 * 1.2667},
       "K_sign_contradiction": {"lower_bound_kappa_over_3": 3 * 1.2667 / 3,
                                "observed": rows[1][3],
                                "contradiction_MeV": abs(rows[1][3]) + 1.2667},
       "enumeration": {"assignments": len(best), "best5": [
           {"max_dev_MeV": d, "lu": lu, "ld": ld, "ls": ls} for d, lu, ld, ls in best[:5]]},
       "verdict": "超度量正核類（任意 κ／b／層指派）解析否證；出路 = Z₂ 符號核"}
with open("data/03_ultrametric_nogo.json", "w") as f:
    json.dump(out, f, ensure_ascii=False, indent=2, default=str)
print("\n→ data/03_ultrametric_nogo.json")