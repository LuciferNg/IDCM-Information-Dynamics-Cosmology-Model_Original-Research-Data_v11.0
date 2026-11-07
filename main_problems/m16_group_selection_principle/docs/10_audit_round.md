# docs/10 — 全鏈審計記錄（2026-11-07）

**觸發**：用戶要求「未審計腳本／未檢查文檔全部審計，確保全鏈一致無漂移」。
**結果**：L1 11/11 ｜ L2 18/18 ｜ L3 56/56 ｜ audit 31/31（FAIL 0）

---

## §1 審計範圍（磁碟實況）

```
scripts/  (15) 00_gap_named_degrees.g ｜ 01_arithmetic_reduction.py ｜ 02_gap_full_scan.g
               02b_chunk_driver.py ｜ 02c_order_counts.g ｜ 03_variants_and_lookelsewhere.py
               04_gap_alternative_principles.g ｜ 05_candidate_principles.py
               06_minimum_cost_principle.py ｜ 07_variational_axiom_attack.py
               08_mf1prime_mf2_attack.py ｜ 09_mf1prime_dissolution.py
               gap_run.sh ｜ gap_setup.sh ｜ run_02_parallel.sh
               審計基建：l2_independent_verify.py ｜ l3_doc_crosscheck.py
                        10_audit_verify.py ｜ 11_l1_reproducibility.py ｜ run_chain.sh
docs/     (9)  01–09（09 = 開放問題權威清單）
data/     全部 run 記錄＋92,570 群掃描原始檔（md5 3dfd5d214f279707fba4390eb42b00e8）
```

## §2 四閘定義與結果

```
L1  scripts/11_l1_reproducibility.py
    每支核心腳本跑兩次，stdout 逐位元組相同、rc=0。
    02 全量掃描例外（~35 min）：以「重生成 md5 ≡ 入庫 md5」作等價證明 —— 本次重跑
    92,570 行、md5 3dfd5d21… 完全一致（11 分片、workers done fail=0）。
    ⟹ 11/11 PASS

L2  scripts/l2_independent_verify.py
    零 import L1、零讀 docs；INPUTS = 鏈常數（B、N_d、L₅、L₄、N_gen、Σγ、⟨W⟩²）＋
    Planck 2018 測量值 ＋ GAP 原始掃描檔（自建 parser）、mpmath dps=40 ＋ 精確 Fraction。
    18 項檢查（含 6 條 identity：3·4·5 = 60 = N_d+B²；γ_3+γ_5 = γ_max = 2γ_4；
    Σd = B²；Σγ = 2(N_d+1)；⟨W⟩²/Σγ = 41/60；Burnside |G| = Σd²）。
    ⟹ 18/18 PASS（內部派生 69 keys）

L3  scripts/l3_doc_crosscheck.py
    期望值只來自 L2；8 份文檔逐 claim 核；細整數一律帶 anchor（R6）；負控制通過。
    ⟹ 56/56 PASS，逐檔命中率 100%（docs/01–08）

audit  scripts/10_audit_verify.py
    ① 腳本集對稱（磁碟 ↔ 宣告；文檔引用全部存在）
    ② 逆向綁定（headline 數字 → producer script 之 log 實際含該值）
    ③ 路徑衛生（scripts 內零 /tmp、.hermes/cache、/var/tmp；GAP 走 data/gap/root）
    ④ 產物綁定（13 個宣告產物全在；無孤兒 log）
    ⑤⑦ 反向數值覆蓋（docs ≥6 位有效數字 → data 池，rel ≤ 1e-3）
    ⑥ driver（run_chain.sh 覆蓋全部核心＋L2/L3/audit，含傳遞覆蓋）
    A8a 舊值掃描 ｜ A10 時序閘 ｜ A12 commit 完整性 ｜ G1 GAP 部署
    ⟹ 31/31 PASS
```

## §3 審計發現之真缺陷（3 項，全部已修）

```
D1【路徑衛生，最嚴重】scripts/gap_run.sh 之 GAPROOT 指向
   /home/wsl/.hermes/cache/scratch/gapdebs/root —— 該目錄 24h 閒置即被 prune。
   實測：本輪審計開始時該目錄**已被清空** ⟹ 00/02/02c/04 四支 GAP 腳本全部不可重現
   （「昨日跑得通」變今日 FileNotFoundError），而所有內容閘仍會照樣綠。
   修：新增 scripts/gap_setup.sh（apt-get download ＋ dpkg-deb -x 到 data/gap/root，
       版本釘 gap 4.14.0-3 / smallgrp 1.5.4-2）；gap_run.sh 改指向記錄內持久位；
       data/gap/ 入 .gitignore（二進位不入庫、可重建）；audit ③ 與 G1 加閘。
   驗：重裝後 GAP 4.14.0 / SmallGrp 1.5.4；全量掃描復現 md5 完全一致。

D2【舊值掃描自我阻擋】docs/06 §「新殘項」與 docs/07 §2.5 之「1/6 = 16.7%」代價
   係 MF1′ 之舊判詞；MF1′ 已於 docs/08 撤回，但原文無標註 ⟹ 讀者會當成活值。
   修：兩處加 append-only 標註（原文不動），指明「已失效，見 docs/08」。

D3【文檔環境段過時】docs/02 §「環境」寫「deb 解包至 scratch」＋「長片 1,705,468 ms」。
   修：append-only 附註（改 data/gap/root；執行時間屬環境量，不入 L2／L3 池）。

另：audit 自身三處誤報已修（①c 跨記錄引用誤判、⑥ 傳遞覆蓋缺 l2/l3 前綴、
    A12a 排除清單漏 L1/L2/L3 輸出）—— 屬審計器 bug，非記錄缺陷，已記錄於此以免後人重踩。
```

## §4 漂移檢查（無漂移）

```
① 掃描數據：92,570 行 / md5 3dfd5d21…（重跑後逐位元組相同）
② 群層統計：17（∏=|G|）｜267（Σ=12）｜1（合取）｜5（去相異）｜586（去非平凡）
            ｜11（DNT={3,4,5}）｜24（DNT ⊇ {3,4,5}）｜5（perfect ≤500）｜13（階 60）
③ 高精度值：docs 內所有 ≥6 位有效數字全部可追溯至 data 池（audit ⑤⑦ PASS）
④ 時序：全部產物 mtime ≥ producer（audit A10 PASS）
⑤ 文檔 claim：56 條全部命中（L3），含 Planck Ω_b / k / z 值與 γ-Δ 精確表
```

## §5 對下游之影響

```
本 O 內：零漂移 ⟹ 無下游需重算（四數全綠即證）
本 O 外（只 handoff，不改他檔 §4.6a）：
  M15  docs/24 L74、final_result L95 之「M16 只開不推」→ 已 handoff（附終態更新）
  M17  README:27、01_problem_record §CQ5「🔴 阻塞（待 M16）」→ 已 handoff（阻塞解除）
  M3   原 R5 handoff 仍有效；補一條 MF2 閉合說明（perfect 非公理）
  上游 O171 三項（OB1–OB3）、基礎層（OB4）→ 版本閘：v11→v10 不推，見 docs/09 §2
```

## §6 重跑方式

```bash
bash scripts/gap_setup.sh      # 一次性（冪等）
bash scripts/run_chain.sh      # L1 → L2 → L3 → audit（約 4 分鐘）
bash scripts/run_chain.sh --heavy   # 另加重跑 92,570 群全量掃描（~35 分鐘）
```