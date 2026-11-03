# scripts/00_gap_named_degrees.g — 具名群之 irrep 維數（GAP 4.14.0 權威讀出）
# 執行（於 m16 目錄）：bash scripts/gap_run.sh scripts/00_gap_named_degrees.g
Print("GAP ", GAPInfo.Version, " smallgrp ", InstalledPackageVersion("smallgrp"), "\n");
for nm in [["C1",Group(())],["S3",SymmetricGroup(3)],["A4",AlternatingGroup(4)],
           ["S4",SymmetricGroup(4)],["A5",AlternatingGroup(5)],["S5",SymmetricGroup(5)],
           ["PSL(2,7)",PSL(2,7)],["A6",AlternatingGroup(6)],["S6",SymmetricGroup(6)],
           ["A7",AlternatingGroup(7)]] do
  g := nm[2];
  Print(nm[1], " |G|=", Size(g), " ncls=", NrConjugacyClasses(g),
        " degrees=", List(Irr(CharacterTable(g)), x -> x[1]), "\n");
od;
Print("--- order 60 : all 13 groups ---\n");
for i in [1..NumberSmallGroups(60)] do
  g := SmallGroup(60,i);
  Print("60/", i, " ", StructureDescription(g), " degrees=",
        List(Irr(CharacterTable(g)), x -> x[1]), "\n");
od;
QUIT;
