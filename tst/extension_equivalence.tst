# Exact ordered comparison; this synthetic model isolates affine choices
# without relying on a classification table or on cochain associativity.
gap> gaugeModel := rec(modelId:="exact-gauge-comparison-test");;
gap> gaugeModel.dimension := function(n) if n<0 then return 0; fi; return 1; end;;
gap> gaugeModel.zero := k->rec(A:=List([1..gaugeModel.dimension(k-3)],i->0),B:=List([1..gaugeModel.dimension(k-2)],i->0),C:=List([1..gaugeModel.dimension(k-1)],i->0),D:=List([1..gaugeModel.dimension(k+1)],i->0));;
gap> gaugeModel.coboundary := function(n,v,signed) if signed and n=0 then return -2*v; elif signed and n=4 then return 2*v; fi; return List([1..gaugeModel.dimension(n+1)],i->0); end;;
gap> gaugeModel.d := function(k,x) if k=4 then return gaugeModel.zero(5); fi; Assert(0,k=3); return rec(A:=-2*x.A,B:=[0],C:=ShallowCopy(x.B),D:=x.C+2*x.D); end;;
gap> gaugeModel.xtimes := function(k,x,y) return rec(A:=x.A+y.A,B:=List(x.B+y.B,z->z mod 2),C:=List(x.C+y.C,z->z mod 2),D:=x.D+y.D+[x.C[1]*y.C[1]]); end;;
gap> gaugeTarget := rec(A:=[2],B:=[0],C:=[0],D:=[2]);;
gap> gaugeCanonical := rec(A:=[0],B:=[0],C:=[1],D:=[0]);;
gap> gaugeProof := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,gaugeTarget,gaugeCanonical);;
gap> Assert(0,gaugeProof.status="computed" and gaugeProof.equalityVerified and gaugeProof.boundaryFlatnessVerified);
gap> Assert(0,gaugeProof.gauge=rec(A:=[-1],B:=[1],C:=[1],D:=[0]));
gap> Assert(0,gaugeProof.boundary=rec(A:=[2],B:=[0],C:=[1],D:=[1]));
gap> Assert(0,gaugeProof.product=gaugeTarget and gaugeProof.leadingEquation.chosenPrimitive=[-1]);
gap> Assert(0,gaugeProof.definingSystemWitness.comparisonEquationVerified and not IsBound(gaugeProof.definingSystemWitness.flatnessVerified));
gap> Assert(0,not IsMutable(gaugeProof.gauge));
gap> # In k3 there is no gauge A: vary the leading binary gauge B instead.
gap> gaugeThree := ShallowCopy(gaugeModel);;
gap> gaugeThree.coboundary := function(n,v,signed) if signed and n=3 then return 2*v; fi; return List([1..gaugeModel.dimension(n+1)],i->0); end;;
gap> gaugeThree.d := function(k,x) if k=3 then return gaugeModel.zero(4); fi; Assert(0,k=2); return rec(A:=[0],B:=[0],C:=ShallowCopy(x.B),D:=x.C+2*x.D); end;;
gap> gaugeThreeProof := KOAHSS_ExtensionGaugeCompare(gaugeThree,3,rec(A:=[0],B:=[0],C:=[1],D:=[1]),gaugeThree.zero(3));;
gap> Assert(0,gaugeThreeProof.status="computed" and gaugeThreeProof.gauge=rec(A:=[],B:=[1],C:=[1],D:=[0]));
gap> Assert(0,gaugeThreeProof.search.leadingChoices=2 and gaugeThreeProof.leadingEquation.modulus=2);
gap> # An integral leading-kernel adjustment is tried after the particular.
gap> gaugeKernel := ShallowCopy(gaugeModel);;
gap> gaugeKernel.coboundary := function(n,v,signed) if signed and n=4 then return 2*v; fi; return List([1..gaugeModel.dimension(n+1)],i->0); end;;
gap> gaugeKernel.d := function(k,x) if k=4 then return gaugeModel.zero(5); fi; Assert(0,k=3); return rec(A:=[0],B:=[0],C:=List(x.A,z->z mod 2),D:=x.C+2*x.D); end;;
gap> gaugeKernelTarget := rec(A:=[0],B:=[0],C:=[1],D:=[1]);;
gap> gaugeKernelProof := KOAHSS_ExtensionGaugeCompare(gaugeKernel,4,gaugeKernelTarget,gaugeKernel.zero(4));;
gap> Assert(0,gaugeKernelProof.status="computed" and gaugeKernelProof.gauge.A=[1] and gaugeKernelProof.winningStage="ABCD" and Last(gaugeKernelProof.attemptedStages).leadingChoices=2);
gap> gaugeLimited := KOAHSS_ExtensionGaugeCompare(gaugeKernel,4,gaugeKernelTarget,gaugeKernel.zero(4),rec(maxLeadingChoices:=1));;
gap> Assert(0,gaugeLimited.status="unresolved" and gaugeLimited.resourceLimit and not gaugeLimited.equalityVerified);
gap> gaugeRadiusZero := KOAHSS_ExtensionGaugeCompare(gaugeKernel,4,gaugeKernelTarget,gaugeKernel.zero(4),rec(integerRadius:=0));;
gap> Assert(0,gaugeRadiusZero.status="unresolved" and not gaugeRadiusZero.resourceLimit and not IsBound(gaugeRadiusZero.gauge));
gap> gaugeNoPrimitive := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,rec(A:=[1],B:=[0],C:=[0],D:=[0]),gaugeModel.zero(4));;
gap> Assert(0,gaugeNoPrimitive.status="unresolved" and not gaugeNoPrimitive.equalityVerified);
gap> # The smallest sufficient gauge support wins, with every higher
gap> # component held at the same literal zero throughout its stage.
gap> gaugeDOnly := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,rec(A:=[0],B:=[0],C:=[0],D:=[2]),gaugeModel.zero(4));;
gap> Assert(0,gaugeDOnly.status="computed" and gaugeDOnly.winningStage="D" and gaugeDOnly.gauge=rec(A:=[0],B:=[0],C:=[0],D:=[1]));
gap> Assert(0,List(gaugeDOnly.attemptedStages,x->x.stage)=["D"] and gaugeDOnly.attemptedStages[1].fixedZeroLayers=["A","B","C"]);
gap> gaugeCD := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,rec(A:=[0],B:=[0],C:=[0],D:=[1]),gaugeModel.zero(4));;
gap> Assert(0,gaugeCD.status="computed" and gaugeCD.winningStage="CD" and gaugeCD.gauge=rec(A:=[0],B:=[0],C:=[1],D:=[0]));
gap> Assert(0,List(gaugeCD.attemptedStages,x->x.stage)=["D","CD"] and gaugeCD.attemptedStages[1].searchComplete);
gap> Assert(0,gaugeThreeProof.winningStage="BCD" and List(gaugeThreeProof.attemptedStages,x->x.stage)=["D","CD","BCD"]);
gap> Assert(0,gaugeProof.winningStage="ABCD" and List(gaugeProof.attemptedStages,x->x.stage)=["D","CD","BCD","ABCD"]);
gap> # Stage selection lets the reducer hold its lower normal form fixed
gap> # while trying the same support across candidate lower coordinates.
gap> gaugeOnlyCD := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,rec(A:=[0],B:=[0],C:=[0],D:=[1]),gaugeModel.zero(4),rec(stages:=["CD"]));;
gap> Assert(0,gaugeOnlyCD.status="computed" and gaugeOnlyCD.winningStage="CD" and Length(gaugeOnlyCD.attemptedStages)=1);
gap> gaugeTooSmall := KOAHSS_ExtensionGaugeCompare(gaugeModel,4,rec(A:=[0],B:=[0],C:=[0],D:=[1]),gaugeModel.zero(4),rec(stages:=["D"]));;
gap> Assert(0,gaugeTooSmall.status="unresolved" and gaugeTooSmall.searchComplete and not gaugeTooSmall.resourceLimit);
gap> # Priority vectors never discard unpreferred full-kernel directions.
gap> gaugePriority := ShallowCopy(gaugeKernel);;
gap> gaugePriority.gaugeKernelRepresentatives := function(n,signed,family) return [List([1..gaugePriority.dimension(n)],i->0)]; end;;
gap> gaugePriorityProof := KOAHSS_ExtensionGaugeCompare(gaugePriority,4,gaugeKernelTarget,gaugeKernel.zero(4));;
gap> Assert(0,gaugePriorityProof.status="computed" and gaugePriorityProof.gauge.A=[1] and gaugePriorityProof.leadingEquation.kernelSearchSource="complete-kernel");
gap> gaugePriority.gaugeKernelRepresentatives := function(n,signed,family) return family.homogeneousGenerators; end;;
gap> gaugePriorityProof := KOAHSS_ExtensionGaugeCompare(gaugePriority,4,gaugeKernelTarget,gaugeKernel.zero(4));;
gap> Assert(0,gaugePriorityProof.status="computed" and gaugePriorityProof.leadingEquation.kernelSearchSource="preferred");
gap> # A nonzero preferred direction is tried and fails. The overlapping
gap> # complete-kernel search must retain its other independent direction.
gap> gaugeNonzeroPriority := rec(modelId:="nonzero-priority-fallback");;
gap> gaugeNonzeroPriority.dimension := function(n) if n<0 then return 0; fi; return 2; end;;
gap> gaugeNonzeroPriority.zero := k->rec(A:=List([1..gaugeNonzeroPriority.dimension(k-3)],i->0),B:=List([1..gaugeNonzeroPriority.dimension(k-2)],i->0),C:=List([1..gaugeNonzeroPriority.dimension(k-1)],i->0),D:=List([1..gaugeNonzeroPriority.dimension(k+1)],i->0));;
gap> gaugeNonzeroPriority.coboundary := function(n,v,signed) return List([1..gaugeNonzeroPriority.dimension(n+1)],i->0); end;;
gap> gaugeNonzeroPriority.d := function(k,x) if k=3 then return gaugeNonzeroPriority.zero(4); fi; Assert(0,k=2); return rec(A:=[0,0],B:=[0,0],C:=ShallowCopy(x.B),D:=[0,0]); end;;
gap> gaugeNonzeroPriority.xtimes := function(k,x,y) return rec(A:=x.A+y.A,B:=List(x.B+y.B,z->z mod 2),C:=List(x.C+y.C,z->z mod 2),D:=x.D+y.D); end;;
gap> gaugeNonzeroPriority.gaugeKernelRepresentatives := function(n,signed,family) return [[1,0]]; end;;
gap> gaugeNonzeroPriorityTarget := rec(A:=[0,0],B:=[0,0],C:=[0,1],D:=[0,0]);;
gap> gaugeNonzeroPriorityProof := KOAHSS_ExtensionGaugeCompare(gaugeNonzeroPriority,3,gaugeNonzeroPriorityTarget,gaugeNonzeroPriority.zero(3));;
gap> Assert(0,gaugeNonzeroPriorityProof.status="computed" and gaugeNonzeroPriorityProof.winningStage="BCD" and gaugeNonzeroPriorityProof.equalityVerified);
gap> Assert(0,gaugeNonzeroPriorityProof.gauge.B=[0,1] and gaugeNonzeroPriorityProof.leadingEquation.kernelSearchSource="complete-kernel");
gap> Assert(0,Last(gaugeNonzeroPriorityProof.attemptedStages).preferredKernelDirections=1 and Last(gaugeNonzeroPriorityProof.attemptedStages).leadingChoices=3);
gap> # The same allowance is shared by all stages, including unsuccessful
gap> # lower searches; failure never fabricates an equality certificate.
gap> gaugeTiny := KOAHSS_ExtensionGaugeCompare(gaugeKernel,4,gaugeKernelTarget,gaugeKernel.zero(4),rec(maxChoices:=1));;
gap> Assert(0,gaugeTiny.status="unresolved" and gaugeTiny.resourceLimit and gaugeTiny.search.differentialEvaluations<=1 and not IsBound(gaugeTiny.gauge));
gap> # Production cochain fixtures: verify the full ordered equality for the
gap> # signed k4 A-doubling relation, including its integral boundary carry.
gap> # These are complete tuples; no expected abelian-group invariant is used.
gap> gaugeActualFour := function()
> local ahss,model,target,canonical,proof;
> ahss:=koAHSS(CyclicGroup(2),[1],[1],4,rec(details:=true));
> model:=KOAHSS_ExtensionBarModel(ahss._context.backend,4);
> target:=rec(A:=[2],B:=[0],C:=[0],D:=[-42]);
> canonical:=rec(A:=[0],B:=[1],C:=[1],D:=[-18]);
> proof:=KOAHSS_ExtensionGaugeCompare(model,4,target,canonical);
> model.close();
> return proof;
> end;;
gap> gaugeActualFourProof := gaugeActualFour();;
gap> Assert(0,gaugeActualFourProof.status="computed" and gaugeActualFourProof.equalityVerified);
gap> Assert(0,gaugeActualFourProof.gauge=rec(A:=[-1],B:=[0],C:=[0],D:=[0]));
gap> Assert(0,gaugeActualFourProof.boundary=rec(A:=[2],B:=[1],C:=[1],D:=[-29]));
gap> Assert(0,gaugeActualFourProof.product=rec(A:=[2],B:=[0],C:=[0],D:=[-42]));
gap> # With no A gauge in k3, the same API compares marked C/D tuples.
gap> gaugeActualThree := function()
> local ahss,model,tuple,proof;
> ahss:=koAHSS(CyclicGroup(2),0,0,3,rec(details:=true));
> model:=KOAHSS_ExtensionBarModel(ahss._context.backend,3);
> tuple:=rec(A:=[0],B:=[0],C:=[1],D:=[1]);
> proof:=KOAHSS_ExtensionGaugeCompare(model,3,tuple,tuple);
> model.close();
> return proof;
> end;;
gap> gaugeActualThreeProof := gaugeActualThree();;
gap> Assert(0,gaugeActualThreeProof.status="computed" and gaugeActualThreeProof.equalityVerified);
gap> Assert(0,gaugeActualThreeProof.gauge=rec(A:=[],B:=[0],C:=[0],D:=[0]));
