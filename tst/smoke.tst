# Standard GAP package smoke check; run with TestPackage("fermionAHSS").
gap> pages := koAHSS(CyclicGroup(2), 0, 0, 1, 5);;
gap> Assert(0, Length(pages) = 5);
gap> expected := [[[0],[],[2]],[[],[]],[[2]],[],[]];;
gap> Assert(0, ForAll(pages, table -> table = expected));
gap> Assert(0, koAHSS(CyclicGroup(2), 0, 0, 1) = expected);
gap> Assert(0, koAHSS(CyclicGroup(2), 0, 0, -1) = [[[0]],[],[],[],[]]);
gap> R := ResolutionFiniteGroup(CyclicGroup(2),4);;
gap> B := koAHSSHAPSpace(R,koAHSSNaturalOperations());;
gap> Assert(0, koAHSS(CyclicGroup(2),[1],[1],1,5) = koAHSSpages(B,[1],[1],1,5));
