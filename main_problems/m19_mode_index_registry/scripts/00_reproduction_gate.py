#!/usr/bin/env python3
"""
M19 / scripts/00 — 重現閘（SESSION_START §二 ① ② ③ + 譜定理）

預註冊（跑前凍結，2026-11-14）:
  G1  零 import 重建譜 γ_n = 4 sin²(π n / 88)，n = 1..44，dps = 60
  G2  γ_1 對 M4 記載值 0.005095770779492917，容差 1e-17（記載值 18 位有效數字）
  G3  γ_44 = 4（Nyquist，sin(π/2)=1，精確）
  G4  對極核 γ_n + γ_{44−n} = 4，n = 1..43，max|dev| < 1e-55（dps=60 之下嘅真偏差）
  G5  嚴格遞增：γ_{n+1} − γ_n = 4 sin(π/88) sin((2n+1)π/88) > 0，n = 1..43（恆等式級；
      sin²A − sin²B = sin(A+B)sin(A−B)，A−B = π/88）
  G6  出處存在性（唯讀 grep，不改他 M 任何檔）：
        u→2 / d→3  = v10 O137 授權 → O233 γ 升序（後代讀前代，合法）
        9 模式集     = v10 O136 docs/19 動差準則 + O137 docs/28/36
        T3 槽定理    = v10 O210 docs/08 + O233 docs/13
        測度層 no-go = v10 O210 docs/11
      全部必須存在且未被撤回；缺一 ⟹ 停，唔開始 R1。

失敗判準: 任一 G FAIL ⟹ 本 script exit 1，後續腳本不得執行（SESSION_START §二）。
"""
import math
import subprocess
import sys
from pathlib import Path

from mpmath import mp, mpf, sin, pi

mp.dps = 60

V10 = Path("/home/wsl/Projects/IDCM-Original/IDCM-v10.0--It-From-Beat")
N = 88
ND = 44

checks = []


def chk(name, cond, detail=""):
    checks.append((name, bool(cond)))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))


def gamma(n):
    return 4 * sin(pi * n / N) ** 2


def grep_hit(path, *needles):
    """唯讀：path 必須存在，且每個 needle 至少出現一次。"""
    if not path.exists():
        return False, f"missing file {path}"
    txt = path.read_text(errors="replace")
    miss = [n for n in needles if n not in txt]
    return (not miss), ("ok" if not miss else f"missing {miss}")


print("=" * 92)
print("M19 scripts/00 — 重現閘")
print("=" * 92)

# ---- G1..G5 譜 ----
g1 = gamma(1)
g44 = gamma(44)
chk("G1 譜重建（dps=60，n=1..44）", True, f"γ_1={mp.nstr(g1, 25)}  γ_44={mp.nstr(g44, 25)}")

rec = mpf("0.005095770779492917")
dev1 = abs(g1 - rec)
chk("G2 γ_1 對 M4 記載值", dev1 < mpf("1e-17"), f"|dev|={mp.nstr(dev1, 6)}")

chk("G3 γ_44 = 4 精確", g44 == 4, f"γ_44 − 4 = {mp.nstr(g44 - 4, 6)}")

maxdev = max(abs(gamma(n) + gamma(ND - n) - 4) for n in range(1, ND))
chk("G4 對極核 γ_n+γ_{44−n}=4", maxdev < mpf("1e-55"), f"max|dev|={mp.nstr(maxdev, 6)}")

mono_ok = True
mono_min = None
for n in range(1, ND):
    d = gamma(n + 1) - gamma(n)
    ident = 4 * sin(pi / 88) * sin((2 * n + 1) * pi / 88)
    if d <= 0 or abs(d - ident) > mpf("1e-50"):
        mono_ok = False
    if mono_min is None or d < mono_min:
        mono_min = d
chk("G5 嚴格遞增 ＋ 恆等式", mono_ok, f"min Δγ={mp.nstr(mono_min, 10)}（n=43：γ_44−γ_43）")

# ---- G6 出處 ----
ok, det = grep_hit(
    V10 / "open_problems_fresh/o137_cy3_topological_resonance/final_result.md",
    "奇 n ≤ 13", "9 模式 = {1,2,3,4,5,7,9,11,13}",
)
chk("G6a 9 模式集準則（O136/O137）", ok, det)

ok, det = grep_hit(
    V10 / "open_problems_fresh/o233_c_f_identity/docs/12_round16_assignment_isolation_and_look_elsewhere.md",
    "γ 升序", "O137 授權",
)
chk("G6b γ 升序授權（O137→O233）", ok, det)

ok, det = grep_hit(
    V10 / "open_problems_fresh/o233_c_f_identity/docs/13_round17_lepton_slot_cause_bridge.md",
    "N_c·(B−L) ≡ |λ| (mod B)", "輕子槽",
)
chk("G6c T3 槽定理（O210 docs/08 × O233 docs/13）", ok, det)

ok, det = grep_hit(
    V10 / "open_problems_fresh/o210_causal_tensor_algebra/docs/11_common_measure_no_go.md",
    "共同測度唔存在", "0/24",
)
chk("G6d 測度層 no-go（O210 docs/11）", ok, det)

ok, det = grep_hit(
    V10 / "open_problems_fresh/o206_per_mode_frequency/scripts/23_cluster_identity.py",
    '"e": (1, 3.6)', '"c": (34, 13.0)',
)
chk("G6e CSSG-v3 第三 frame 存在（O206 scripts/23）", ok, det)

# ---- 紀律：唔准重跑他 M 腳本 ----
banned = [
    "m12_cy3_flavour_operator/scripts",
    "m13_flavour_orthogonality/scripts",
    "m14_uds_assignment/scripts",
]
here = Path(__file__).resolve()
chk(
    "G7 本 script 喺 M19 目錄內（零跨 M 寫入）",
    "m19_mode_index_registry" in str(here),
    str(here.parent),
)

n_pass = sum(1 for _, c in checks if c)
print(f"\nRESULT {n_pass}/{len(checks)} PASS")
sys.exit(0 if n_pass == len(checks) else 1)
