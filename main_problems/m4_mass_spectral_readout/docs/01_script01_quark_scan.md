# M4 docs/01 — `scripts/01_quark_scan.py`：前設問題與結論報告

> 腳本：`scripts/01_quark_scan.py`（2026-10-21 建立，/tmp 遷入）
> 執行 log：`data/run_01_quark_scan.log`（sha256 `c7e8b76ec2dfdb136a74cc1823ce5455117db753d8f8706fc12357c16241179b`，76 行）
> 首跑 crash 記錄：`data/run_01_quark_scan.pre_fix.log`（sha256 `a69f2894391843817f066c1687ef412779c90b637d2dc5cb0f0fa1951fac670a`）
> Canonical run（2026-10-24 05:34，**假設從未執行過，全支重跑**）：`PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 scripts/01_quark_scan.py` → **exit 0，102 ms，76 行**
> 狀態：🔴 單點譜比路線（本腳本部分）全滅

## 一、前設問題（腳本 docstring 原文）

```
夸克差值閉合測試 01 — 質量數據 + 譜比掃描
相同方法 (R7 定案): 質量 = 譜讀數差值/比, 唔 fit 絕對標度。
夸克質量 (PDG 2024, MS̄):
  m_u = 2.16 MeV, m_d = 4.67 MeV, m_s = 93.4 MeV  @ 2 GeV
  m_c = 1.27 GeV @ m_c ; m_b = 4.18 GeV @ m_b ; m_t = 172.69 GeV @ pole
問: 質量比 (m_i/m_j) 有冇 γ_n 譜結構? (同 R7 A/B/C/D 族)
    同步邊界讀數 = γ_n? 1/γ_n? γ_n²? Σ/γ?
```

## 二、事前判準（掃描前寫死）

- 族 A：`γ_a/γ_b`，`a,b ∈ 1..44`（K_chain L=88 譜，1936 值）
- 族 B：跨比 look-elsewhere 密度 —— 全部 15 比、30 個方向，逐個報最優 dev
- 族 C：`m_q = γ_n² · c_k`，`c_k ∈ S ∪ {γ_k} ∪ {1/γ_k}`，44 mode × 99 候選 = 4356 組合/夸克
- MC null：200 隨機夸克質量（log10 ∈ [log10 2e-3, log10 173]）
- 判定：`|dev| < 0.1%` 為命中門檻；命中數要同池密度比

## 三、原始輸出要點（引自 log，未改）

| 段 | 結果 |
|:--|:--|
| [0] | log10 質量比矩陣（15 比）|
| [A] | `m_c/m_s = 13.5974` vs `γ20/γ5 = 13.6032` **dev 0.0424%**；`m_s/m_d = 20` vs γ35/γ6 = 19.9245 dev 0.3773%；`m_d/m_u = 2.16204` vs γ33/γ19 dev 0.2728% |
| [B] | `\|dev\|<0.1%` 嘅比數 **2/30**（兩者 = `m_s/m_c` 與 `m_c/m_s`，同一 ratio）；`\|dev\|<1%` **14/30** |
| [C] | 最優組合：m_s `γ23²·γ2` dev 0.0653%、m_c `γ22²·γ8` dev **0.0022%**、m_b `γ43²·1/γ38` dev 0.0330%；`<0.1%` 組合數 0–2 /夸克 |
| [C] MC | 隨機夸克 `best<0.1%` = **44.0%**，`best<1%` = **96.0%** |

## 四、結論報告

1. **族 C 有 0.0022% 級「命中」，但完全冇預言力**：每夸克要一個獨立 `A_q = m/γ²`，6 個夸克 = 6 個自由參數 ⟹ 零自由度 ⇒ 純 fit。MC null 44%/96% 證明「池大到任何點都撞得到」。
2. **族 A 唯一 <0.1% 命中 = `m_c/m_s ≈ γ20/γ5`**（0.042%），但 30 個方向中 14 個 <1% ⟹ 池密度高；正式判定交 02（MC null）與 09（統一 protocol）。
3. **判定**：本腳本**唔提供任何可升格閉合**；作為閘（gate）用途成立 —— 它證實「夸克單點質量比 vs γ池」係高密度池。
4. **[C] 段 bug（已修，附註）**：/tmp 原版用「A_q 與 γ 向量同 index 配對再格式化」⟹ `TypeError` 全支 crash（見 `pre_fix.log`）。改為 well-posed 自由 index 掃描後，**[0]/[A]/[B] 60 行與原版 byte-identical** ⟹ 只有 [C] 語義被修正，其餘結果未變。
5. **數據潔癖聲明**：本腳本質量 = PDG 2024 **引用值**（m_d ±10.3%）；任何 <0.1% 級結論在此精度下**不成立**，需 lattice（見 11）。

## 五、與其他文檔／腳本關係

- → `02` 追擊唯一 0.1% 命中（m_c/m_s）
- → `03` 測「讀數×譜強度」形式
- → `09` 統一 density protocol 給最終 PN（本腳本 [B] 密度 → 09 的 PN 值）
- → `11` lattice 重判（PDG 引用值誤差主導嘅兩個命中）
- 記錄於 `docs/01–08`（閘重建）、`docs/09_script09_density_gate.md`（統一判決）

---

## 六、閘核對 — 原 docs/r01 §一「閘結果（期望值 vs 實跑）」

> 閘 = 遷入 scripts 逐支重跑，同 SESSION_START §一 期望值對照。本支結果：

| 期望（SESSION_START §一） | 實跑（2026-10-24 04:56） | 判 |
|:--|:--|:--|
| `m_c/m_s` vs γ20/γ5 dev 0.0424%；池 density 高 | [A] `13.5974` → `13.6032` **dev 0.0424%**；[B] `\|dev\|<0.1%` 比數 **2/30** | 🟢 |

**閘判**：8 支全部跑通；7 支逐位吻合；01 首跑 crash（腳本 bug，非物理，見 `docs/01` §七）；05 MC 數字與記錄不符（見 `docs/05` §七）。

---

## 附錄（原 docs/r01）：R1 推導鏈重建（8 支遷入 scripts 逐支驗證）(2026-10-24)

> 全文照錄（原標題層級降一級；一符一字未改）。內含 8 支閘表（§一）、01 [C] TypeError 取證（§二）、05 MC protocol 不符（§三）、本輪新增知識（§四）。

> 狀態：🟢 閘通過（8/8 跑通）＋ 🟡 兩項誠實修正（01 crash、05 MC protocol 不一致）
> 腳本：`scripts/01_quark_scan.py` … `scripts/08_span_12_judge.py`（原 8 支，未改計算內容）
> 執行：`/usr/bin/python3`、`PYTHONDONTWRITEBYTECODE=1`；每支 stdout 落 `data/run_NN_*.log`
> 環境：2026-10-24 04:56 HKT，WSL Debian，numpy + mpmath + sympy

---

### 一、閘結果（期望值 vs 實跑）

| Script | SESSION_START §一 期望 | 實跑（本輪） | 判 |
|:--|:--|:--|:--:|
| 01_quark_scan | m_c/m_s γ20/γ5 dev 0.0424%；池 density 高 | [A] `m_c/m_s = 13.5974` → γ20/γ5=13.6032 **dev 0.0424%**；[B] `\|dev\|<0.1%` 比數 **2/30** | 🟢 |
| 02_mc_mcs_judge | 撞池 97.9%＝後驗判死 | MC 5000 族：**97.9%** | 🟢 逐位對 |
| 03_reading_strength | m=γ²v 全部 >6%；γv 1.3–1.7% 事後揀 | γ²·v 最小 dev **6.42%**（m_b）、最大 195.7%；γ·v：c 1.29%、t 1.67% | 🟢 |
| 04_lepton_anchor | μ/e 1.15%、τ/μ 1.04%、τ/e 77%；MC 16.8% | 1.1524%／1.0394%／77.4256%；MC 3000：**16.8%** | 🟢 |
| 05_hadron_dm | n−p/m_e γ29/γ16 0.0195%（density 殺） | Δm(n−p)/m_e=2.53092 → γ29/γ16=2.530426 **dev 0.0195%**；MC null **43.3%**（記錄寫 8.7%，見 §三） | 🟡 |
| 06_hop_probe | s/d=20 池零偏差 | `s/d = 20.000000` → hop 結構池 20 **dev 0.000%** | 🟢 |
| 07_hit_followup | 誤差帶 [17.4,22.6] | **[17.409, 22.591]**（m_d ±10.28%） | 🟢 |
| 08_span_12_judge | R=11.994978 dev 0.042%；MC 0.166% | **R=11.994978, dev 0.0419%**；MC 50000：**0.166%** | 🟢 逐位對 |

**閘判**：8 支全部跑通；7 支數值逐位吻合 §一 期望；01 首跑 crash（腳本 bug，非物理）→ 修復後 [0]/[A]/[B] 60 行 **byte-identical**（`diff` 已驗）。05 嘅 MC 數字與記錄不符 → 見 §三。

---

### 二、F1 — scripts/01 [C] 段 TypeError（腳本 bug，2026-10-24 修）

**症狀**（首跑 exit=1）：`TypeError: unsupported format string passed to numpy.ndarray.__format__`
於 `scripts/01_quark_scan.py:118`。

**根因**：/tmp 原版寫法：
```python
Aq = m/g2                    # 44 值向量（逐 mode 嘅 A_q）
best = np.argmin(np.abs(Aq - inv_g))   # 同 index 配對：γ_n² vs 1/γ_n
print(f"... {b1[0]:.4g} (dev {abs(b1[0]/Aq-1)*100:.2f}%)")
```
① `Aq` 係向量 ⟹ `abs(b1[0]/Aq-1)*100` 係向量 ⟹ `:.2f` 直接 TypeError；
② 更根本：`argmin(|A_q(n) − 1/γ_n|)` 係**同 index 硬配**，冇物理意義（A_q 用 mode n、候選又用同一 n）。

**修法（well-posed 化，非降低標準）**：改為自由 index 掃描
`∃(n,k): γ_n²·c_k = m_q`，`c_k ∈ S ∪ {γ_k} ∪ {1/γ_k}`（99 候選 × 44 mode = 4356/夸克），
並加 MC density（200 隨機夸克質量）。
- 結果：m_s 最好 γ23²·γ2 dev 0.065%；m_c γ22²·γ8 dev 0.0022%；**best<0.1% 組合數 0–2**
- **density：隨機夸克質量 best<0.1% = 44.0%、best<1% = 96.0%** ⟹ 池密度隨機級
  ⟹ 「A_q 有譜身份」🔴（同 03 §E「冇預言力」結論一致，本輪首次量化）

**取證**：`data/run_01_quark_scan.pre_fix.log`（崩潰原文，保存）／`data/run_01_quark_scan.log`（修後）；
backup：`~/.hermes/backups/scripts__01_quark_scan.py.v22pre.20261024_045821`。

---

### 三、F2 — scripts/05 MC null protocol 與記錄不符（8.7% vs 43.3%）

| 項 | 記錄（01_problem_record §2.4） | 遷入腳本實跑 | docstring 聲稱 |
|:--|:--|:--|:--|
| MC null（γ比池命中） | **8.7%** | **43.3%** | 「隨機質量差對 (0.5..6 MeV log-均勻)」 |

- 腳本碼用域 `log10(0.1)…log10(10)`、thr 0.1%；docstring 寫 0.5–6 MeV ⟹ **碼／文檔不一致**
- 網格掃（lo×hi×thr×單比/族級，3000 MC × 168 組）：域 0.1–10/thr 0.015% 給 8.7%；
  域 0.5–6/thr 0.1% 給 55.2% ⟹ **8.7% 唔可以由遷入腳本唯一復現**（係另一個版本嘅 protocol）
- 三支腳本 MC 準則互不可比：02 族級 thr 0.1%（質量域）→97.9%；04 單比 thr 0.05%→16.8%；
  05 單比 thr 0.1%（差比域）→43.3%
- **結論**：05 嘅「density 殺」裁定**不變**（43.3% > 8.7% ⟹ 更強否定），但記錄數字係**不可復現嘅舊 protocol**。
  統一 protocol 見 `docs/09_density_unified.md`（`scripts/09_density_gate.py`）：
  以 **命中自身 dev** 做 thr、以受測比值域做 null 域 ⟹ n−p/m_e 命中 PN = **30.98%** 🔴

---

### 四、本輪新增知識（唔改舊結論）

1. 01 [C] 係全支唯一 crash 點 ⟹ SESSION_START §一「8/8 跑通」在 /tmp 原版**未曾達成**（記錄膨脹已修）
2. 「γ比池密度」需要用**單一 protocol** 定義，否則 sector 之間數字唔可比（已立 `scripts/09`）
3. m₃ 標度歸一化：`m₃(v=246)=53.5858 meV`（R7 記錄 53.586 ✓）；`m₃(v_SM=246.2196)=53.6815 meV`（+0.179%）
   ⟹ 引用 R7 數字時**必須標 v 版本**（P6 docs/12 已改用 v_SM）

---

*append-only；2026-10-24 R1 建立。*
