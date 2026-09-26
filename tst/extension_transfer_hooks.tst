# Native curvature cannot be applied to gauges as though it were a boundary.
# This triangular toy model has C^2 carrying once into D, and delta_D=2.
gap> transferHookActions:=0;; transferHookDivisions:=0;;
gap> transferHookModel:=rec(modelId:="test-native-transfer",certificateLevel:="transfer-R",dimension:=n->Maximum(0,Minimum(1,n+1)));;
gap> transferHookModel.zero:=k->rec(A:=List([1..transferHookModel.dimension(k-3)],i->0),B:=List([1..transferHookModel.dimension(k-2)],i->0),C:=List([1..transferHookModel.dimension(k-1)],i->0),D:=[0]);;
gap> transferHookModel.matrix:=function(n,signed) if n<0 then return []; elif n=3 then return [[2]]; else return [[0]]; fi; end;;
gap> transferHookModel.coboundary:=function(n,v,signed) if n<0 then return [0]; else return v*transferHookModel.matrix(n,signed); fi; end;;
gap> transferHookModel.lift:=function(n,v,signed) return ShallowCopy(v); end;;
gap> transferHookModel.project:=transferHookModel.lift;;
gap> transferHookModel.d:=function(k,x) if k<>3 then Error("native curvature was incorrectly used on a gauge"); fi; return transferHookModel.zero(k+1); end;;
gap> transferHookModel.xtimes:=function(k,x,y) return rec(A:=x.A+y.A,B:=List(x.B+y.B,z->z mod 2),C:=List(x.C+y.C,z->z mod 2),D:=x.D+y.D+[QuoInt(x.C[1]+y.C[1],2)]); end;;
gap> transferHookModel.act:=function(k,e,c) local out; Assert(0,k=3); transferHookActions:=transferHookActions+1; out:=rec(A:=ShallowCopy(c.A),B:=ShallowCopy(c.B),C:=ShallowCopy(c.C),D:=c.D+2*e.D); return out; end;;
gap> transferHookModel.divideLeft:=function(k,left,total) local right; Assert(0,k=3); transferHookDivisions:=transferHookDivisions+1; right:=rec(A:=total.A-left.A,B:=List(total.B-left.B,z->z mod 2),C:=List(total.C-left.C,z->z mod 2),D:=total.D-left.D); right.D:=right.D-[QuoInt(left.C[1]+right.C[1],2)]; return right; end;;
gap> transferHookCanonical:=transferHookModel.zero(3);; transferHookCanonical.C:=[1];;
gap> transferHookTarget:=StructuralCopy(transferHookCanonical);; transferHookTarget.D:=[2];;
gap> transferHookProof:=KOAHSS_ExtensionGaugeCompare(transferHookModel,3,transferHookTarget,transferHookCanonical);;
gap> Assert(0,transferHookProof.status="computed" and transferHookProof.gauge.D=[1]);
gap> Assert(0,transferHookProof.certificateLevel="transfer-R" and transferHookProof.actionFlatnessVerified);
gap> Assert(0,not IsBound(transferHookProof.boundary) and transferHookProof.product=transferHookTarget);
gap> Assert(0,transferHookProof.search.equivalenceScope="resolution gauges");
gap> transferHookRight:=KOAHSS_ExtensionDivideLeft(transferHookModel,3,transferHookCanonical,transferHookTarget);;
gap> Assert(0,transferHookDivisions=1 and transferHookModel.xtimes(3,transferHookCanonical,transferHookRight)=transferHookTarget);
gap> transferHookLayers:=rec(A:=rec(name:="A",p:=0,orders:=[],generators:=[],cochains:=[]),B:=rec(name:="B",p:=1,orders:=[],generators:=[],cochains:=[]),C:=rec(name:="C",p:=2,orders:=[2],generators:=[1],cochains:=[[1]]),D:=rec(name:="D",p:=4,orders:=[2],generators:=[1],cochains:=[[1]]));;
gap> transferHookOracle:=KOAHSS_ExtensionHigherOracle(rec(),3,transferHookLayers,transferHookModel);;
gap> transferHookResult:=koAHSSExtensionFromLayers(transferHookLayers,transferHookOracle);;
gap> Assert(0,transferHookResult.status="computed" and transferHookResult.invariants=[4]);
gap> Assert(0,transferHookActions>4 and transferHookDivisions>1);
gap> transferHookRelation:=First(transferHookResult.extensionVectors,v->v.layer="C").result.witness.reduction;;
gap> Assert(0,transferHookRelation.canonicalComparison.certificateLevel="transfer-R");
gap> Assert(0,ForAll(transferHookRelation.reductionSteps,step->step.boundary.equation="state = act(gauge, canonical)"));
gap> # Public options and supplied resolutions preserve the same backend basis.
gap> transferHookResolution:=ResolutionFiniteGroup(CyclicGroup(2),5);;
gap> transferHookFull:=koFull(transferHookResolution,0,0,2,rec(extensionModel:="transfer"));;
gap> Assert(0,transferHookFull.invariants=[[0],[],[2,2],[2,2]]);
gap> Assert(0,IsIdenticalObj(transferHookFull.ahss._context.resolution,transferHookResolution));
gap> Assert(0,transferHookFull.extensionModel="transfer" and ForAll(transferHookFull.degreeResults,r->r.modelSelection.selected="low"));
gap> koFull(CyclicGroup(2),0,0,1,rec(extensionModel:="unknown"));
Error, koFull: extensionModel must be bar or transfer
gap> koFull(transferHookFull.ahss,rec(extensionModel:=1));
Error, koFull: extensionModel must be bar or transfer
gap> koFull(transferHookFull.ahss,rec(other:=true));
Error, koFull: options must be a record containing only extensionModel
gap> # A setup refusal must select the reference model and retain its reason.
gap> transferSavedFactory:=KOAHSS_ExtensionTransferredModel;;
gap> MakeReadWriteGlobal("KOAHSS_ExtensionTransferredModel");
gap> KOAHSS_ExtensionTransferredModel:=function(backend,k) return rec(status:="unresolved",reason:="controlled native setup refusal"); end;;
gap> transferFallbackFull:=koFull(CyclicGroup(2),0,0,3,rec(extensionModel:="transfer"));;
gap> KOAHSS_ExtensionTransferredModel:=transferSavedFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionTransferredModel");
gap> Assert(0,transferFallbackFull.status="computed");
gap> Assert(0,transferFallbackFull.degreeResults[5].modelSelection.selected="bar" and transferFallbackFull.degreeResults[5].modelSelection.fallback);
gap> Assert(0,transferFallbackFull.degreeResults[5].transferAttempt.reason="controlled native setup refusal");
gap> # A native finite audit alone cannot silently certify gauge completeness.
gap> transferBarWrapper:=function(backend,k)
> local wrapped;
> wrapped:=KOAHSS_ExtensionBarModel(backend,k);
> if wrapped.status="computed" then
>     wrapped.act:=function(degree,gauge,canonical) return wrapped.xtimes(degree,wrapped.d(degree-1,gauge),canonical); end;
>     wrapped.classEquivalenceVerified:=false;
> fi;
> return wrapped;
> end;;
gap> MakeReadWriteGlobal("KOAHSS_ExtensionTransferredModel");
gap> KOAHSS_ExtensionTransferredModel:=transferBarWrapper;;
gap> transferUncertifiedFull:=koFull(transferFallbackFull.ahss,rec(extensionModel:="transfer"));;
gap> KOAHSS_ExtensionTransferredModel:=transferSavedFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionTransferredModel");
gap> Assert(0,transferUncertifiedFull.invariants=transferFallbackFull.invariants);
gap> Assert(0,transferUncertifiedFull.degreeResults[5].modelSelection.selected="bar");
gap> Assert(0,transferUncertifiedFull.degreeResults[5].transferAttempt.pendingLayer="transfer-certification");
gap> Assert(0,transferUncertifiedFull.degreeResults[5].transferAttempt.candidatePresentation.status="computed");
gap> # The identity wrapper really uses the same complete-bar equations.
gap> MakeReadWriteGlobal("KOAHSS_ExtensionTransferredModel");
gap> KOAHSS_ExtensionTransferredModel:=function(backend,k) local wrapped; wrapped:=transferBarWrapper(backend,k); if wrapped.status="computed" then wrapped.classEquivalenceVerified:=true; fi; return wrapped; end;;
gap> transferIdentityFull:=koFull(transferFallbackFull.ahss,rec(extensionModel:="transfer"));;
gap> KOAHSS_ExtensionTransferredModel:=transferSavedFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionTransferredModel");
gap> Assert(0,transferIdentityFull.invariants=transferFallbackFull.invariants);
gap> Assert(0,transferIdentityFull.degreeResults[5].modelSelection.selected="transfer" and not transferIdentityFull.degreeResults[5].modelSelection.fallback);
gap> Assert(0,transferIdentityFull.degreeResults[5].certificateLevel="transfer-R");
gap> # Preserve reference-model cache reuse and the zero/free-layer behavior.
gap> transferSavedBarFactory:=KOAHSS_ExtensionBarModel;; transferBarSetupCount:=0;;
gap> MakeReadWriteGlobal("KOAHSS_ExtensionBarModel");
gap> KOAHSS_ExtensionBarModel:=function(backend,k) transferBarSetupCount:=transferBarSetupCount+1; return rec(status:="unresolved",reason:="controlled bar resource refusal"); end;;
gap> transferBarZeroFree:=koFull(CyclicGroup(3),0,0,4);;
gap> KOAHSS_ExtensionBarModel:=transferSavedBarFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionBarModel");
gap> Assert(0,transferBarSetupCount=1);
gap> Assert(0,transferBarZeroFree.degreeResults[5].status="computed" and transferBarZeroFree.invariants[5]=[0,3]);
gap> Assert(0,transferBarZeroFree.degreeResults[6].status="computed" and transferBarZeroFree.invariants[6]=[]);
