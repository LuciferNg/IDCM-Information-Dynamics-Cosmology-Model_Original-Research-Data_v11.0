"""M7 scripts/21 — M2 身份層傳播回核：互換律之字面讀法 ⟹ 非可分族被**結構性**排除
==========================================================================================================
M2（2026-11-23 04:18）裁決：採字面讀法 `E_sync(σ) := m(1−σ)/M₀`，`λ≠0` 非可分族被交換律排除；
例 `σ_n = 1 iff n≡0 mod 3` 偏差 `668.31876 MeV`；`T_M2′` 前提升為「交換律結構級」。

本檔獨立重算（**閉式路線**，與 scripts/16 之雙迴圈路線不同）：
  §1 `A_σ = ((Σσγ)² − Σσ²γ²)/(2Σγ)`（精確代數恆等）⟹ 偏差 = `M₀(A_σ + A_{1−σ})`
  §2 例（mod 3）偏差值與 M2 之 `668.31876` 比對
  §3 **定理**：`λ≠0` ⟹ **所有**遮罩皆違反（`dev=0` 需兩側各 ≤1 模，而合計 = `N_d` = 44 ⟹ 不可能）⟹ M2 之「排除」係結構，非 MC
輸出 data/21_m2_reading_verify.json
"""
import json
import os

from mpmath import mp, mpf, sin, pi, nstr

mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
M7 = os.path.abspath(os.path.join(HERE, ".."))
OUT = {"module": "M7", "script": "21_m2_reading_verify.py", "checks": []}
NP = [0, 0]


def ok(name, cond, got, want, note=""):
    NP[0] += 1
    NP[1] += 1 if cond else 0
    print(f"  {'✅' if cond else '🔴'} {name:<52} {got}" + (f" | want {want}" if want else "") + (f"  {note}" if note else ""))
    OUT["checks"].append(dict(item=name, pass_=bool(cond), got=str(got), want=str(want), note=note))


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


NM = 44
G = [mpf(0)] + [4 * sin(n * pi / 88) ** 2 for n in range(1, NM + 1)]
SG = sum(G[1:])
M0 = mpf("27.900915")
W = sum(G[n] * G[m] / SG for n in range(1, NM + 1) for m in range(n + 1, NM + 1))

print("=" * 100)
print("M7 scripts/21 — M2 互換律字面讀法（閉式路線）")
print("=" * 100)

hdr("§1 閉式：A_σ = ((Σσγ)² − Σσγ²)/(2Σγ)")
def A_closed(mask):
    s = sum(G[n] for n in range(1, NM + 1) if mask(n))
    q = sum(G[n] ** 2 for n in range(1, NM + 1) if mask(n))
    return (s * s - q) / (2 * SG)


def A_loop(mask):
    idx = [n for n in range(1, NM + 1) if mask(n)]
    return sum(G[a] * G[b] / SG for i, a in enumerate(idx) for b in idx[i + 1:])


m3 = lambda n: n % 3 == 0
d_cl, d_lp = A_closed(m3), A_loop(m3)
print(f"    mod-3 遮罩：A_σ(閉式) = {nstr(d_cl, 12)}；A_σ(雙迴圈) = {nstr(d_lp, 12)}")
ok("兩路線一致（<1e−45）", abs(d_cl - d_lp) < mpf("1e-45"), nstr(abs(d_cl - d_lp), 3), "<1e−45")

hdr("§2 偏差 = M₀(A_σ + A_{1−σ})（λ=−1）vs M2 之 668.31876 MeV")
dev = M0 * (A_closed(m3) + A_closed(lambda n: not m3(n)))
print(f"    dev = {nstr(dev, 12)} MeV；M2 報 668.31876 MeV")
ok("偏差與 M2 一致（<1e−4 MeV）", abs(dev - mpf("668.31876")) < mpf("1e-4"), nstr(dev, 10), "668.31876")
ok("量級 O(10²) MeV ⟹ λ 族被強烈排除", dev > 100, nstr(dev, 6) + " MeV", ">100 MeV")

hdr("§3 定理：λ≠0 ⟹ 凡有 ≥1 對佔用模之遮罩皆違反（結構性，非統計）")
print("    A_σ ≥ 0；A_σ = 0 ⟺ σ 佔用模 ≤ 1。dev = M₀|λ|(A_σ + A_{1−σ}) = 0 ⟺ 兩側皆 ≤1 模")
print("    但 σ 與 1−σ 之佔用模數合計 = N_d = 44 ⟹ 兩側皆 ≤1 不可能 ⟹ **dev > 0 對所有遮罩**")
print("    ⟹ λ≠0 族被交換律**普遍**排除（不止非可分族）")
cnt_zero, cnt_pos, worst = 0, 0, mpf(0)
for nocc in range(0, NM + 1):          # 以佔用模數族分類（完備掃描 0..44）
    mask = lambda n, k=nocc: n <= k
    a1, a2 = A_closed(mask), A_closed(lambda n, k=nocc: n > k)
    dv = M0 * (a1 + a2)
    if dv == 0:
        cnt_zero += 1
    else:
        cnt_pos += 1
    worst = max(worst, dv)
print(f"    0..44 佔用模數掃描：dev=0 之族 = {cnt_zero}（k≤1 兩端）；dev>0 = {cnt_pos}；max dev = {nstr(worst, 8)} MeV")
ok("dev = 0 之族 = 0（N_d=44 ⟹ 兩側皆≤1 不可能）", cnt_zero == 0, cnt_zero, 0,
   "⟹ λ≠0 **普遍**被排除 ⟹ 結構性，非 200 組 MC 之經驗結果")
ok("λ_min = −Σγ/W 與 M2 C3 一致", abs(-SG / W - mpf("-2.0694941")) < mpf("5e-8"), nstr(-SG / W, 9), "-2.0694941")
ok("W = 43.48888889（M2 值）", abs(W - mpf("43.4888888889")) < mpf("1e-9"), nstr(W, 12), "43.4888888889")

print("\n" + "=" * 100)
print(f"結果：{NP[1]}/{NP[0]} PASS")
OUT.update(npass=NP[1], ntotal=NP[0],
           metrics=dict(A_m3=float(d_cl), dev_mod3=float(dev), zero_families=int(cnt_zero),
                        max_dev=float(worst), lambda_min=float(-SG / W), W=float(W)),
           verdict="literal reading: exchange law kills every non-separable family structurally "
                   "(dev=0 only at k<=1 trivial ends), not merely 200/200 MC")
json.dump(OUT, open(os.path.join(M7, "data", "21_m2_reading_verify.json"), "w"), ensure_ascii=False, indent=1)
print("JSON → data/21_m2_reading_verify.json")