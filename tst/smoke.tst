# Finite-group checks with known answers (previous package, same formulas).
gap> pages := koAHSS_batch(CyclicGroup(2), 0, 0, 1, 5);;
gap> Assert(0, Length(pages) = 5 and pages[5] = [[[0],[],[2]],[[],[]],[[2]],[],[]]);
gap> Assert(0, koAHSS_batch(CyclicGroup(2), 0, 0, -1) = [[[0]],[],[],[],[]]);
gap> Assert(0, koAHSS(CyclicGroup(2), 0, 0, 1).lines = [[[2],[],[2],fail,fail]]);
gap> Assert(0, koFull(CyclicGroup(2), 0, 0, 3).invariants = [0, 8]);
gap> c4 := koFull_batch(CyclicGroup(4), [1], 0, 6);;
gap> Assert(0, c4.invariants = [[],[2],[2],[2,4],[2],[4],[],[4]]);
gap> Assert(0, ForAll(c4.degreeResults, d -> d.heavyMeasurements = 0));
gap> Assert(0, koFull(CyclicGroup(3), 0, 0, 5).invariants = [9]);
gap> Assert(0, koFull(CyclicGroup(9), 0, 0, 5).invariants = [3, 27]);
gap> d8 := koFull(ResolutionFiniteGroup(DihedralGroup(8), 9), [1,1], [1,0,1], 3);;
gap> Assert(0, d8.invariants = [2, 2, 8] and d8.degreeResult.certificateLevel = "light-R");
gap> Assert(0, koFull(CyclicGroup(4), 0, 0, 5).invariants = [4]);
gap> Assert(0, koFull(CyclicGroup(4), 0, [1], 5).invariants = [2, 32]);
gap> z3z := ResolutionDirectProduct(ResolutionFiniteGroup(CyclicGroup(3), 9), ResolutionAbelianGroup([0], 9));;
gap> Assert(0, koFull(z3z, 0, 0, 6).invariants = [9]);
