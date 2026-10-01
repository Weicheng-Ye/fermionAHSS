# Light rows (doc/extensions.md, "Light rows"): the relations that neither
# the target-layer shortcut nor a primary operation settles are read from
# residues of transported defining data, with no flat lift, reflected product
# or gauge search. FERMIONAHSS_LIGHT_RELATIONS=0, or the override below,
# measures them in the transferred model instead; the groups agree.
gap> lightRow := function(result,layer) return First(result.extensionVectors,v->v.layer=layer); end;;
gap> # C4 with s=[1], degrees -1 to 6.
gap> lightC4 := koFull_batch(CyclicGroup(4),[1],0,6);;
gap> Assert(0,lightC4.invariants=[[],[2],[2],[2,4],[2],[4],[],[4]]);
gap> # Degree two: B over D in the strong absorption regime (D=Z/2), the
gap> # row (L,0) of 2b=c with no evaluation and no worker.
gap> lightTwo := lightC4.degreeResults[4];;
gap> Assert(0,lightTwo.certificateLevel="light-R" and not IsBound(lightTwo.modelId));
gap> Assert(0,lightRow(lightTwo,"B").result.witness.light.precision="ETA");
gap> # Degree four: A over C by Theorem C (omega=0), 2a=c.
gap> Assert(0,lightRow(lightC4.degreeResults[6],"A").result.witness.model="light-R");
gap> # Degree six: B over D, the page form with a D gauge, 2b=d.
gap> lightSix := lightRow(lightC4.degreeResults[8],"B");;
gap> Assert(0,lightSix.result.lowerCoordinates=[1] and lightSix.result.witness.light.precision="DS");
gap> # D8 in degree three (s=[1,1], omega=[1,0,1] on this resolution): two B
gap> # generators with dependent leading C parts, D elementary; the pivot is
gap> # absorbed and the other row is the formal sum of a kernel atom and the
gap> # pivot. The model gives the same group.
gap> lightR := ResolutionFiniteGroup(DihedralGroup(8),9);;
gap> lightD8 := koFull(lightR,[1,1],[1,0,1],3);;
gap> Assert(0,lightD8.invariants=[2,2,8] and lightD8.degreeResult.certificateLevel="light-R");
gap> Assert(0,Set(List(Filtered(lightD8.degreeResult.extensionVectors,v->v.layer="B"),
>     v->v.result.witness.light.precision))=["ETA"]);
gap> # The prime three in degree five: 3a = 2Y with rho_3 Y = (rho_3 A)^3, Z/9.
gap> lightZ3 := koFull(CyclicGroup(3),0,0,5);;
gap> Assert(0,lightZ3.invariants=[9] and lightRow(lightZ3.degreeResult,"A").result.witness.model="light-R");
gap> # The switch restores the measurement in the model.
gap> KOAHSS_EXTENSION_RELATION_OVERRIDE.light := false;;
gap> heavySix := koFull(CyclicGroup(4),[1],0,6);;
gap> Unbind(KOAHSS_EXTENSION_RELATION_OVERRIDE.light);
gap> Assert(0,heavySix.invariants=[4] and heavySix.degreeResult.certificateLevel="transfer-R");
gap> Assert(0,lightRow(heavySix.degreeResult,"B").result.witness.model="complete");
gap> Assert(0,ForAll(lightC4.degreeResults,d->d.heavyMeasurements=0 and not IsBound(d.lightFallbacks)));
gap> Assert(0,lightD8.degreeResult.heavyMeasurements=0 and not IsBound(lightD8.degreeResult.lightFallbacks));
gap> Assert(0,lightZ3.degreeResult.heavyMeasurements=0 and not IsBound(lightZ3.degreeResult.lightFallbacks));
gap> Assert(0,heavySix.degreeResult.heavyMeasurements>0);
gap> Read(Filename(DirectoriesPackageLibrary("fermionAHSS","tst"),"extension_light_frame.g"));
gap> Read(Filename(DirectoriesPackageLibrary("fermionAHSS","tst"),"extension_light_failures.g"));
