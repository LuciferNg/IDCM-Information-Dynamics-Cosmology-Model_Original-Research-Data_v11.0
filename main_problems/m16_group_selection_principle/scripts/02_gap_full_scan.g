# scripts/02_gap_full_scan.g — M16 R2：全域掃描（階 ≤ 500 之全部群）
# 執行：bash scripts/gap_run.sh scripts/02_gap_full_scan.g <outfile>
# 輸出每行：n  i  [[deg,mult],...]   （CharacterDegrees，GAP 內建）
# 用途：事後以 Python（scripts/03）計算準則及其變體之命中率。

tot := 0;
str := OutputTextFile(ARGV[1], false);
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
Print("TOTAL ", tot, " ms=", Runtime()-t0, "\n");
Print("GAP_VERSION ", GAPInfo.Version, "\n");
CloseStream(str);
QUIT;