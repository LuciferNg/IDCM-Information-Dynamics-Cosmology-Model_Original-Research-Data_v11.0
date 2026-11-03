# scripts/04_gap_alternative_principles.g — M16：替代選擇原理之檢定（GAP 4.14.0）
# 執行（於 m16 目錄）：bash scripts/gap_run.sh scripts/04_gap_alternative_principles.g
# §A 階 ≤ 500 之全部「非阿貝爾單群」（最小非阿貝爾單群原理之候選集）
# §B 階 60 之 13 群：簡單性 / 完備性(G=G') / 線性特徵標數 k=|G:G'|
# §C 階 ≤ 500 之全部「完備群」中，階 60 者（完備性 + 階 60 之替代原理）

Print("GAP ", GAPInfo.Version, " smallgrp ", InstalledPackageVersion("smallgrp"), "\n");
Print("--- §A non-abelian simple groups, order <= 500 ---\n");
cnt := 0;
for n in [1..500] do
  if n mod 2 = 0 and not IsPrimePowerInt(n) and SmallGroupsAvailable(n) then
    for i in [1..NumberSmallGroups(n)] do
      g := SmallGroup(n,i);
      if IsSimple(g) and not IsAbelian(g) then
        Print("  SIMPLE ", n, "/", i, " ", StructureDescription(g), "\n");
        cnt := cnt + 1;
      fi;
    od;
  fi;
od;
Print("non-abelian simple (order<=500) count = ", cnt, "\n");

Print("--- §B order 60 : all 13 groups ---\n");
for i in [1..NumberSmallGroups(60)] do
  g := SmallGroup(60,i);
  degs := List(Irr(CharacterTable(g)), x -> x[1]);
  k := Number(degs, d -> d = 1);
  Print("60/", i, " ", StructureDescription(g),
        " simple=", IsSimple(g), " perfect=", IsPerfect(g),
        " |G'|=", Size(DerivedSubgroup(g)), " k(linear)=", k,
        " DNT=", Set(Filtered(degs, d -> d > 1)), "\n");
od;

Print("--- §C perfect groups, order<=120 (count only) ---\n");
pc := 0;
for n in [1..120] do
  if SmallGroupsAvailable(n) then
    pc := pc + Number(AllSmallGroups(Size, n, IsPerfect, true));
  fi;
od;
Print("  perfect groups order<=120 : ", pc, "\n");
QUIT;