# The E6 D quotient removes an incoming higher-layer image which is not an
# ordinary integral coboundary. Its surviving generator must keep its mark.
gap> stagedReductionModel := rec(modelId:="higher-gauge-reduction-test");;
gap> stagedReductionModel.dimension := function(n) if n<0 then return 0; elif n>=4 then return 2; fi; return 1; end;;
gap> stagedReductionModel.zero := k->rec(A:=List([1..stagedReductionModel.dimension(k-3)],i->0),B:=List([1..stagedReductionModel.dimension(k-2)],i->0),C:=List([1..stagedReductionModel.dimension(k-1)],i->0),D:=List([1..stagedReductionModel.dimension(k+1)],i->0));;
gap> stagedReductionModel.coboundary := function(n,v,signed) if n=4 and signed then return 2*v; fi; return List([1..stagedReductionModel.dimension(n+1)],i->0); end;;
gap> stagedReductionModel.d := function(k,x) local y; y:=stagedReductionModel.zero(k+1); if k=3 then y.D:=2*x.D+[x.C[1],0]; fi; return y; end;;
gap> stagedReductionModel.xtimes := function(k,x,y) return rec(A:=x.A+y.A,B:=List(x.B+y.B,a->a mod 2),C:=List(x.C+y.C,a->a mod 2),D:=x.D+y.D); end;;
gap> stagedD := stagedReductionModel.zero(4);; stagedD.D:=[0,1];; MakeImmutable(stagedD);;
gap> stagedLayers := rec(D:=rec(orders:=[2],fullLifts:=[rec(state:=stagedD)]));;
gap> stagedLower := rec(generatorCount:=1,generatorIds:=["D:1"],presentationId:="marked-test-D",layers:=[rec(name:="D",orders:=[2],startColumn:=1)]);;
gap> stagedTarget := stagedReductionModel.zero(4);; stagedTarget.D:=[1,1];;
gap> Assert(0,koAHSSSolveIntegerSystem([[0,1],[2,0],[0,2]],stagedTarget.D)=fail);
gap> stagedReduction := KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedTarget,stagedLayers,stagedLower);;
gap> Assert(0,stagedReduction.status="computed" and stagedReduction.lowerCoordinates=[1]);
gap> Assert(0,stagedReduction.lowerPresentationId=stagedLower.presentationId);
gap> Assert(0,stagedReduction.canonicalComparison.winningStage="CD");
gap> Assert(0,List(stagedReduction.gaugeSearchAttempts,r->r.stage)=["D","D","CD","CD"]);
gap> Assert(0,stagedReduction.canonicalComparison.gauge.C=[1] and stagedReduction.canonicalComparison.gauge.D=[0,0]);
gap> Assert(0,stagedReduction.canonicalComparison.equalityVerified and stagedReduction.canonicalComparison.boundaryFlatnessVerified);
gap> Assert(0,IsIdenticalObj(stagedLayers.D.fullLifts[1].state,stagedD) and stagedD.D=[0,1]);
gap> # Known page-projection coordinates still require an actual gauge proof.
gap> stagedProjected := KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedTarget,stagedLayers,stagedLower,[1]);;
gap> Assert(0,stagedProjected.status="computed" and stagedProjected.canonicalComparison.winningStage="CD");
gap> stagedWrong := KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedTarget,stagedLayers,stagedLower,[0]);;
gap> Assert(0,stagedWrong.status="unresolved" and not IsBound(stagedWrong.lowerCoordinates));
gap> # After recording a C generator with its D carry, the remaining D
gap> # projects to zero in E6 but is not an ordinary coboundary. Compare the
gap> # original full target, not just that residual, with the canonical word.
gap> stagedC := stagedReductionModel.zero(4);; stagedC.C:=[1];; stagedC.D:=[0,1];; MakeImmutable(stagedC);;
gap> stagedLayers.C:=rec(orders:=[2],fullLifts:=[rec(state:=stagedC)]);;
gap> stagedLowerC := rec(generatorCount:=2,generatorIds:=["D:1","C:1"],presentationId:="marked-test-DC",layers:=[rec(name:="D",orders:=[2],startColumn:=1),rec(name:="C",orders:=[2],startColumn:=2)]);;
gap> stagedOriginal := stagedReductionModel.zero(4);; stagedOriginal.C:=[1];; stagedOriginal.D:=[1,1];;
gap> stagedCarry := KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedOriginal,stagedLayers,stagedLowerC,[0,1]);;
gap> Assert(0,stagedCarry.status="computed" and stagedCarry.lowerCoordinates=[0,1]);
gap> Assert(0,stagedCarry.canonicalComparison.target=stagedOriginal and stagedCarry.canonicalLowerProduct=stagedC);
gap> Assert(0,stagedCarry.canonicalComparison.boundary.D=[1,0] and stagedCarry.canonicalComparison.winningStage="CD");
gap> # An unavailable coordinate search must not invent free coordinates.
gap> stagedFree := StructuralCopy(stagedLower);; stagedFree.layers[1].orders:=[0];;
gap> Assert(0,KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedTarget,stagedLayers,stagedFree).status="unresolved");
gap> stagedTooLarge := rec(generatorCount:=1,layers:=[rec(name:="D",orders:=[64],startColumn:=1)]);;
gap> Assert(0,KOAHSS_ExtensionGaugeReduce(stagedReductionModel,4,stagedTarget,stagedLayers,stagedTooLarge).status="unresolved");
gap> # Incomplete oracle evidence survives through the abstract wrapper.
gap> stagedUnresolved := koAHSSExtensionFromLayers(rec(D:=[2],C:=[2]),function(l,i,m,h) return rec(status:="unresolved",reason:="test unresolved comparison",gaugeSearchAttempts:=[rec(stage:="CD")]); end);;
gap> Assert(0,stagedUnresolved.pendingRelation.gaugeSearchAttempts[1].stage="CD");
