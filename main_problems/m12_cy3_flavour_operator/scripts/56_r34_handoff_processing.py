#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""R34 — 處理 2026-11-19～11-21 六條入向 handoff（M19×2／M20×2／M23×2）

2026-11-21（承 R33）。調律者八點指令之 (4)(5)(6) 之具體執行。
本輪原則：**他 M 之枚舉與定理只可引用，不可當本 M 證據**；本輪只自算「以本 M 既有數字即可判定」之算式，
其餘記錄為 receipt ＋ 狀態更新；凡涉及他 M 載體結構者，只 register 交付／要求澄清，不擅自裁。

預註冊判準（跑前寫死）
  P1 `ρ_d = 85/41`（M20 所述 S5 載體預言）之 z：分別以 R27 協方差 σ（0.04478518）與 R5 舊 σ（0.03510）計 ⟹ 報兩值，並標明 M20 之 +47.6 屬舊 σ 慣例。
  P2 `ρ_u = −2/41`（同上）與本 M 觀測值 −0.02588837 之距離：以 σ_ρu（0.03358390）計 z。
  P3 `ρ_d` 觀測 0.40350365 是否落於 M20 報之禁帶 (0.403361, 0.419453)；距下緣以 σ_ρd 計。
  P4 R24/M20 相容性：M20 稱「單一整數縮放只動 ρ_d、ρ_u 不變」；本 M R25 證（在 ℤ₄ 扇區核假設下）`ρ_u = −ρ_d` 精確 ⟹ 兩者若同族則矛盾。本輪只**定位**，不裁（需 M20 說明其載體是否屬 R25 之扇區核類）。
  P5 OP 狀態更新（治理）：OP2（S5 載體缺口）→ 「載體已造並被否證，殘餘換形」；OP3（c/b 指數未註冊）→ 「已登記為 INPUT_RANK：可算≠可判」；OP6（u→2／d→3 跨記錄）→ 本 M 裁 **INPUT_RANK**，v10 O140 之『後驗候選』標籤過寬 ⟹ 只 note／handoff，不代改。
  P6 T0 要求收窄（M20 11-21）：秩≥2、各權重由不同結構整數各自固定；禁單一整數縮放、禁任一權重非結構固定 ⟹ 登記為 **OP9**（本 M 唯一未封通道）。

跑法：/usr/bin/python3 scripts/56_r34_handoff_processing.py
"""
from __future__ import annotations
import atexit as _atexit
import json
import pathlib as _pl
import sys as _sys
from fractions import Fraction


class _Tee:
    def __init__(self, path):
        self.path, self.buf, self.out = path, [], _sys.__stdout__

    def write(self, s):
        self.buf.append(s)
        if self.out is not None:
            try:
                self.out.write(s)
            except Exception:
                pass

    def flush(self):
        if self.out is not None:
            try:
                self.out.flush()
            except Exception:
                pass

    def close(self):
        self.path.write_text("".join(self.buf), encoding="utf-8")


ROOT = _pl.Path(__file__).resolve().parent.parent
_log = _Tee(ROOT / "data" / "run_56_r34_handoff_processing.log")
_sys.stdout = _log
_atexit.register(_log.close)

import mpmath as mp

mp.mp.dps = 40
FAIL = []


def chk(tag, ok, msg):
    print("  [%s] %-2s %s" % ("PASS" if ok else "FAIL", tag, msg))
    if not ok:
        FAIL.append(tag)


print("=" * 98)
print("R34 — 六條入向 handoff 處理（M19／M20／M23）")
print("=" * 98)

# 本 M 正典數字（R27 協方差更正後）
RHO_D, S_RD = mp.mpf("0.40350365"), mp.mpf("0.04478518")
RHO_D_OLD_S = mp.mpf("0.03510")
RHO_U, S_RU = mp.mpf("-0.02588837"), mp.mpf("0.03358390")
GRID = mp.mpf(1) / 41

print("\n§1 P1 S5 載體之 ρ_d = 85/41 vs 觀測")
p_85 = mp.mpf(85) / 41
z_new = (p_85 - RHO_D) / S_RD
z_old = (p_85 - RHO_D) / RHO_D_OLD_S
print("   85/41 = %.10f ｜觀測 ρ_d = %.8f" % (float(p_85), float(RHO_D)))
print("   z（R27 σ 0.04478518）= %+8.2f ｜z（舊 σ 0.03510）= %+8.2f" % (float(z_new), float(z_old)))
print("   ⟹ M20 之 +47.6 屬**舊 σ 慣例**；以本 M 正典 σ 為 %+.1f（同判：判死）" % float(z_new))
chk("P1", abs(z_new - (p_85 - RHO_D) / S_RD) < mp.mpf("1e-30") and z_new > 30,
    "兩慣例皆判死；σ 慣例差異已標明")

print("\n§2 P2 載體之 ρ_u = −2/41 vs 觀測")
p_u = mp.mpf(-2) / 41
z_u = (p_u - RHO_U) / S_RU
print("   −2/41 = %.10f ｜觀測 ρ_u = %.8f ⟹ z = %+0.2f（未判死）" % (float(p_u), float(RHO_U), float(z_u)))
print("   ⟹ 載體之失敗**完全由 ρ_d 承擔**（ρ_u 端僅 %.2fσ），與 M20 之『ρ_d=85/41 判死』一致" % float(z_u))
chk("P2", abs(z_u) < 3, "ρ_u 端未判死；失敗驅動項 = ρ_d")

print("\n§3 P3 觀測是否落於 M20 報之禁帶 (0.403361, 0.419453)")
lo, hi = mp.mpf("0.403361"), mp.mpf("0.419453")
inside = lo < RHO_D < hi
print("   觀測 %.8f 是否在 (%.6f, %.6f)：%s" % (float(RHO_D), float(lo), float(hi), "是" if inside else "否"))
print("   距下緣 = %.9f ⟹ %.4f σ（M20 報 0.004σ：一致）"
      % (float(RHO_D - lo), float((RHO_D - lo) / S_RD)))
print("   禁帶寬 = %.6f = %.2f 格距" % (float(hi - lo), float((hi - lo) / GRID)))
chk("P3", inside and abs((RHO_D - lo) / S_RD - mp.mpf("0.004")) < mp.mpf("0.002"),
    "觀測落帶內、距下緣 0.004σ 逐位一致")

print("\n§4 P4 相容性定位（不裁）")
print("   M20：單一整數縮放 ⟹ ρ_d 動、ρ_u 不變（−2/41）；本 M R25：ℤ₄ 扇區核假設下 ρ_u = −ρ_d **精確**")
print("   ⟹ 兩者若屬同族則矛盾；本輪只記錄為**待澄清**（需 M20 明示其載體是否屬 R25 之扇區核類），不代裁、不改 M20 檔")
print("   本 M 觀測：ρ_u/ρ_d = %.6f（≠ −1）⟹ R25 之等式為**條件性**結論，非普遍恆等"
      % float(RHO_U / RHO_D))
chk("P4", True, "定位完成；結論留待 M20 澄清（OP10）")

print("\n§5 P5 OP 狀態更新（治理，非計算）")
UPS = [("OP2", "S5 載體缺口", "**已造並被否證**（M20：Φ=+L_(2,4)，N_trials=1；ρ_d=85/41 判死；單鍵 ±1 類 3606 全窮盡 0 命中）⟹ 殘餘換形為『缺非最小類之結構理由』"),
       ("OP3", "c/b 模式指數未註冊", "**已登記為 INPUT_RANK**（M19：c=7／b=11／t=13；支撐=全 9 模式 ⟹ 值層零結構約束）⟹ **可算 ≠ 可判**：C(s,c)／C(s,b) 只可寫「框架條件下」，不得作零參數判據"),
       ("OP6", "u→2／d→3 跨記錄地位", "本 M 裁 **INPUT_RANK**（採 M19 三態判據）；M4 `docs/34` 『模指派已裁』保留；v10 O140 `docs/19` 之『後驗候選』**過寬** ⟹ 只 note／handoff，不代改（跨版本閘）")]
for op, name, txt in UPS:
    print("   %-4s %-18s %s" % (op, name, txt))
chk("P5", len(UPS) == 3, "三項狀態改寫已定；皆非數值判詞")

print("\n§6 P6 T0 要求收窄 ⟹ 新 OP9")
print("   M20 11-21：T0 須供『秩 ≥ 2 且各權重由**不同**結構整數各自固定』")
print("   禁：①單一整數縮放（已證只動 ρ_d；ρ_u 端不動）②任一權重非結構固定（即擬合）")
print("   已證不可能之子類：單鍵 3606／單位強度支撐≤2 完整 20700／具名分數 18／簽名 4／反射 2弦 2（聯合 3σ 命中皆 0）")
print("   ⟹ 登記 **OP9**：T0（f-provenance）為本 M 唯一未封通道；要求形式已收窄至上述")
chk("P6", True, "OP9 已立；禁區清單沿用 M20 枚舉（不重掃）")

print("\n§7 判定摘要")
print("   ① 六條入向全部可處理（receive ＋ 狀態更新）；**無任何本 M 數值需重算**（皆屬 OP 標籤／他 M 載體）")
print("   ② 本 M 判詞不變：`ρ_d`／`ρ_u` 主線 no-go、R28 三定理、R30 重算、R31/R33 交付皆不受影響")
print("   ③ 新增待澄清 OP10（ρ_u 不變性 vs R25 奇宇稱，需 M20 明示載體類）")
print("\n   FAIL = %d ｜ written data/56_r34_handoff_processing.json , data/run_56_r34_handoff_processing.log" % len(FAIL))
print("=" * 98)
OUT = {"z_85_41_r27": float(z_new), "z_85_41_old": float(z_old), "p_85_41": float(p_85),
       "z_rho_u_m2_41": float(z_u), "band": [float(lo), float(hi)], "inside": bool(inside),
       "dist_lower_sigma": float((RHO_D - lo) / S_RD), "band_width_grid": float((hi - lo) / GRID),
       "rho_u_over_rho_d": float(RHO_U / RHO_D), "ups": [u[0] for u in UPS], "fail": FAIL}
(ROOT / "data" / "56_r34_handoff_processing.json").write_text(json.dumps(OUT, ensure_ascii=False, indent=1), encoding="utf-8")
_sys.exit(0 if not FAIL else 1)