# scripts/02c_order_counts.g — 各階之群數（負載平衡用；資料來源 SmallGrp 1.5.4）
for n in [1..500] do
  if SmallGroupsAvailable(n) then Print(n, " ", NumberSmallGroups(n), "\n"); fi;
od;
QUIT;
