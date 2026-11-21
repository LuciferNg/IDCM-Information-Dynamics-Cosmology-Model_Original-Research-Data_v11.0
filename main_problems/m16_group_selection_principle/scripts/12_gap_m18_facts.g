# scripts/12_gap_m18_facts.g — M18 三向張力裁決所需之群論事實（GAP 4.14.0 / SmallGrp 1.5.4）
# 輸出：stdout（由 scripts/12_m18_faithfulness_adjudication.py 解析 → data/run_12.log）
Print("GAP ", GAPInfo.Version, " smallgrp ", InstalledPackageVersion("smallgrp"), "\n");
A5 := AlternatingGroup(5);
Print("simple=", IsSimpleGroup(A5), " order=", Size(A5), "\n");
Print("classes=", List(ConjugacyClasses(A5), Size), "\n");
# §A 置換作用：|Hom(A5, S_k)|
for k in [1..6] do
  Print("hom k=", k, " count=", Length(AllHomomorphisms(A5, SymmetricGroup(k))), "\n");
od;
# §B 線性作用：不可約維數 ＋ 忠實性（kernel = 1）
irr := Irr(A5);
for i in [1..Length(irr)] do
  Print("irr deg=", irr[i][1], " faithful=", Size(KernelOfCharacter(irr[i])) = 1,
        " char=", List([1..Length(irr[i])], j -> irr[i][j]), "\n");
od;
Print("min_faithful_linear_deg=",
      Minimum(List(Filtered([1..Length(irr)], i -> Size(KernelOfCharacter(irr[i])) = 1),
                   i -> irr[i][1])), "\n");
QUIT;
