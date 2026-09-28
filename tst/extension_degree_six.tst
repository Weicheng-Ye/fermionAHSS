# Opt-in degree-six comparison; not part of TestPackage (about five minutes):
# gap -q --quitonbreak -c 'LoadPackage("fermionAHSS"); Assert(0,Test(Filename(DirectoriesPackageLibrary("fermionAHSS","tst"),"extension_degree_six.tst"),rec(compareFunction:="uptowhitespace"))); QUIT;'
gap> # Degree six on the one-cell C2 bar, whose comparison is an isomorphism:
gap> # the native model and the complete-bar section model agree on states
gap> # with A zero. The C generator doubles into D, and the B generator
gap> # doubles into C with an integral D carry.
gap> transferR:=ResolutionFiniteGroup(CyclicGroup(2),9);;
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS([1],0,8);;
gap> transferModel:=KOAHSS_ExtensionTransferredModel(transferBackend,6);;
gap> transferBar:=KOAHSS_ExtensionBarModel(transferBackend,6);;
gap> Assert(0,transferModel.supports(6) and transferBar.supports(6));
gap> transferC:=transferModel.zero(6);; transferC.C:=[1];;
gap> transferB:=transferModel.zero(6);; transferB.B:=[1];;
gap> Assert(0,transferModel.d(6,transferC)=transferModel.zero(7) and transferBar.d(6,transferC)=transferModel.zero(7));
gap> Assert(0,transferModel.d(6,transferB)=transferModel.zero(7) and transferBar.d(6,transferB)=transferModel.zero(7));
gap> # The native product omits the pure-C normalization, an integral
gap> # coboundary, so its D layer agrees with the reference model's up to a
gap> # D-coboundary; the lower layers agree exactly.
gap> transferSameClass:=function(x,y) return ForAll(["A","B","C"],f->x.(f)=y.(f)) and koAHSSSolveIntegerSystem(transferBar.matrix(6,true),x.D-y.D)<>fail; end;;
gap> transferSquare:=transferModel.xtimes(6,transferC,transferC);;
gap> Assert(0,transferSquare.C=[0] and transferSquare.D mod 2=1);
gap> Assert(0,transferSameClass(transferSquare,transferBar.xtimes(6,transferC,transferC)));
gap> transferSquare:=transferModel.xtimes(6,transferB,transferB);;
gap> Assert(0,transferSquare.B=[0] and transferSquare.C=[1] and transferSquare.D mod 2=1);
gap> Assert(0,transferSameClass(transferSquare,transferBar.xtimes(6,transferB,transferB)));
gap> Assert(0,transferModel.divideLeft(6,transferB,transferSquare)=transferB);
gap> transferGauge:=transferModel.zero(5);; transferGauge.C:=[1];; transferGauge.D:=[2];;
gap> Assert(0,transferSameClass(transferModel.act(6,transferGauge,transferC),transferBar.xtimes(6,transferBar.d(5,transferGauge),transferC)));
gap> transferModel.close();; transferBar.close();;
