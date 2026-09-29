# Opt-in degree-six comparison; not part of TestPackage (about six minutes):
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
gap> Assert(0,transferSquare.C=[0] and transferSquare.D[1] mod 2=1);
gap> Assert(0,transferSameClass(transferSquare,transferBar.xtimes(6,transferC,transferC)));
gap> transferSquare:=transferModel.xtimes(6,transferB,transferB);;
gap> Assert(0,transferSquare.B=[0] and transferSquare.C=[1] and transferSquare.D[1] mod 2=1);
gap> Assert(0,transferSameClass(transferSquare,transferBar.xtimes(6,transferB,transferB)));
gap> Assert(0,transferModel.divideLeft(6,transferB,transferSquare)=transferB);
gap> transferGauge:=transferModel.zero(5);; transferGauge.C:=[1];; transferGauge.D:=[2];;
gap> Assert(0,transferSameClass(transferModel.act(6,transferGauge,transferC),transferBar.xtimes(6,transferBar.d(5,transferGauge),transferC)));
gap> transferModel.close();; transferBar.close();;
gap> # The three-local model in degree six, checked by suspension: the
gap> # degree-six window of Z/3 x Z is the suspension of the degree-five
gap> # window of Z/3 (A=Z/3 at (3,0), D=Z/3 at (7,-4)), so the group is Z/9.
gap> susR:=ResolutionDirectProduct(ResolutionFiniteGroup(CyclicGroup(3),9),ResolutionAbelianGroup([0],9));;
gap> susFull:=koFull(susR,0,0,6);;
gap> Assert(0,susFull.invariants=[9]);
gap> susA:=First(susFull.degreeResult.extensionVectors,v->v.layer="A");;
gap> Assert(0,susA.result.witness.model="three-local" and susA.result.witness.prime=3);
gap> Assert(0,susA.result.witness.measuredLayers=["D"] and susA.result.lowerCoordinates[1] mod 3<>0);
gap> Assert(0,susFull.degreeResult.layers.A.fullLifts[1].model="three-local");
