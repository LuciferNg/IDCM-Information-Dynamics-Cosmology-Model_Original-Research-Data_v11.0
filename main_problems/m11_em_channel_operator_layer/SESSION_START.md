# M11 SESSION_START — em 通道算子層（U5 攻堅）

## 起點
M4 `docs/19` 診斷：`em` 通道之瓶頸係**概念**（缺一個結構層），唔係精度。
M4 `docs/37 §一 B4` 於 2026-10-28 版本紀律（禁止 v11→v10 反向推送）下懸空 ⟹ M11 承接。

## 本 session 已完成
1. `E_ch` 8 通道由 PDG 2024 獨立重算（與 M4 `docs/18` 一致，max|Δ| = 3e-5）。
2. **F1**：`E_ch` 唔係純 em（K 通道強 IB 殘餘 −3.639 > em +2.265）⟹ 目標更正為「IB 算子」。
3. **F2**：QN 差矩陣 rank = 2（`Fraction` 精確）；QN 盲區 spread 5.664 MeV ⟹ 獨立再證 O124。
4. **F3**：pairwise 環核（M5 `g(d)`，零參數）純 em 側 −11.1%；混合對象 2 通道反號（缺強 IB 項）。
5. **F4**：`R/X ∈ [−1.32, +1.36]`，spread 2.68 ⟹ 單一 `⟨O⟩` 不可行（獨立再證 M4 `docs/13`）。
6. **F5**：算子 = `em[環核] + strong-IB[f(ch)]` 兩項；攻堅面收窄。

## 下一步（先立判準）
N1 取 lattice 純 em 全套 → N2 求 `f(ch)` 結構形式（零新參數，盲測 7 通道）→ N3 全滅則 no-go。

## 紀律
- 一律 `/usr/bin/python3`、`PYTHONDONTWRITEBYTECODE=1`；rank 用 `Fraction`。
- 禁以 `E_ch` 當純 em；禁反解 mode 指派；禁再引單一 `⟨O⟩`（見 `docs/01 §六`）。
