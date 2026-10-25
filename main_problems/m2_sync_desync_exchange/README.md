# M2 — 同步態 ↔ 脫同步態互換公式（質能公式）

## 狀態: 🟡 partial — 首輪推導完成 (2026-10-24)

互換律 `m/M₀ + E_sync = Σγ = 90`（零新常數）；P1✅ P2✅ P3🟡（上界✅/值🔴）P4🟡；絕對標度 🔴。59/59 checks PASS。詳見 `docs/01_exchange_law.md`。

## 檔案

| 檔案 | 描述 | 狀態 |
|:----|:-----|:----:|
| `01_problem_record.md` | 問題陳述 + 已確立結構 + 攻堅路線 | 已註冊 |
| `SESSION_START.md` | 起手式 + 接手起點 | 已註冊 |
| `docs/01_exchange_law.md` | 首輪推導：帳本/算子/測試 + 評級 | 🟡 |
| `scripts/01_ledger_exact.py` | 帳本精確核對 (16/16) | ✅ |
| `scripts/02_exchange_operator.py` | 互換算子 + 質量泛函 (11/11) | ✅ |
| `scripts/03_null_rigidity.py` | 尖截斷/MC null/剛性/Planck (7/7) | ✅ |
| `scripts/04_uniqueness_scaleblind.py` | p-族唯一性/標度盲/beat 恆等式 (7/7) | ✅ |
| `scripts/05_beat_duality.py` | 振幅讀法互換律＋對偶 (6/6) | ✅ |
| `scripts/06_rate_bound.py` | 速率上界＋值 no-go (5/5) | ✅ |
| `scripts/07_phase_slip_link.py` | O7 相位滑移鏈接 (7/7) | ✅ |
| `data/01-07_*.json` | 七腳本輸出 | ✅ |

## 依賴

O121 (sync sector E_sync=61.5) + O237 (質量=脫相位) + O204 (⟨W⟩²=61.5) + O212 (dof=1) + O73 (c 第一性原理) + O215 (Ω 預算) + M1 (質量檔案 schema)

## 註冊

pid `M2`（created 2026-10-15）
