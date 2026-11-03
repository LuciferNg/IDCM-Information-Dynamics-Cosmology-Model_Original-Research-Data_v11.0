# scripts/02_gap_full_scan.g — M16 R2：全域掃描（階 ≤ 500 之全部群）
# 執行（於 m16 目錄）：bash scripts/gap_run.sh scripts/02_gap_full_scan.g
# 輸出：data/run_02_scan_all_500.txt（每行：n  i  [[deg,mult],...]）
#       標準輸出：每個階之進度 + 總數 + GAP 版本
# 註：本機 GAP 未綁 ARGV，故路徑寫死。

outfile := "data/run_02_scan_all_500.txt";
tot := 0;
str := OutputTextFile(outfile, false);
SetPrintFormattingStatus(str, false);
t0 := Runtime();
for n in [1..500] do
  if SmallGroupsAvailable(n) then
    k := NumberSmallGroups(n);
    for i in [1..k] do
      g := SmallGroup(n,i);
      PrintTo(str, n, " ", i, " ", CharacterDegrees(g), "\n");
      tot := tot + 1;
    od;
    Print(n, " ", k, " cum=", tot, " ms=", Runtime()-t0, "\n");
  else
    Print(n, " NOT-AVAILABLE\n");
  fi;
od;
CloseStream(str);
Print("TOTAL ", tot, " ms=", Runtime()-t0, "\n");
Print("GAP_VERSION ", GAPInfo.Version, "\n");
Print("SMALLGRP ", InstalledPackageVersion("smallgrp"), "\n");
QUIT;