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
gap> # The ordinary public API preserves an explicitly supplied resolution.
gap> transferHookResolution:=ResolutionFiniteGroup(CyclicGroup(2),5);;
gap> transferHookFull:=koFull(transferHookResolution,0,0,2);;
gap> Assert(0,transferHookFull.invariants=[[0],[],[2,2],[2,2]]);
gap> Assert(0,IsIdenticalObj(transferHookFull.ahss._context.resolution,transferHookResolution));
gap> Assert(0,transferHookFull.gaugeCompletenessAssumed and not IsBound(transferHookFull.extensionModel));
gap> koFull(CyclicGroup(2),0,0,1,rec(extensionModel:="transfer"));
Error, usage: koFull(group or HAP resolution,s,omega,k) or koFull(detailedE6Result)
gap> koFull(transferHookFull.ahss,rec(extensionModel:="bar"));
Error, usage: koFull(group or HAP resolution,s,omega,k) or koFull(detailedE6Result)
gap> koFull(transferHookFull.ahss,rec());
Error, usage: koFull(group or HAP resolution,s,omega,k) or koFull(detailedE6Result)
gap> # C8 exceeds the complete-bar budget. Its native result requires no bar
> # construction, certification hook, or exceptional C2 basis isomorphism.
gap> transferSavedBarFactory:=KOAHSS_ExtensionBarModel;;
gap> MakeReadWriteGlobal("KOAHSS_ExtensionBarModel");
gap> KOAHSS_ExtensionBarModel:=function(backend,k) Error("koFull attempted complete-bar construction"); end;;
gap> transferSavedFactory:=KOAHSS_ExtensionTransferredModel;; transferNativeCloses:=0;;
gap> MakeReadWriteGlobal("KOAHSS_ExtensionTransferredModel");
gap> KOAHSS_ExtensionTransferredModel:=function(backend,k)
> local native,close;
> native:=transferSavedFactory(backend,k);
> if native.status="computed" then
>     Assert(0,not IsBound(native.certifyPresentation));
>     close:=native.close;
>     native.close:=function() transferNativeCloses:=transferNativeCloses+1; close(); end;
> fi;
> return native;
> end;;
gap> transferLargeR:=ResolutionFiniteGroup(CyclicGroup(8),6);;
gap> transferLargeFull:=koFull(transferLargeR,0,0,3);;
gap> Assert(0,transferLargeFull.status="computed" and transferLargeFull.invariants[5]=[0,2,16]);
gap> Assert(0,IsIdenticalObj(transferLargeFull.ahss._context.resolution,transferLargeR));
gap> Assert(0,transferLargeFull.degreeResults[5].certificateLevel="transfer-R" and transferLargeFull.degreeResults[5].gaugeCompletenessAssumed);
gap> Assert(0,transferLargeFull.degreeResults[5].algebraAudit.status="computed" and transferNativeCloses=1);
gap> Assert(0,not IsBound(transferLargeFull.degreeResults[5].modelSelection) and not IsBound(transferLargeFull.degreeResults[5].barCertification));
gap> # Setup refusal remains unresolved, without a complete-bar retry.
gap> transferRefusalCloses:=0;;
gap> KOAHSS_ExtensionTransferredModel:=function(backend,k) return rec(status:="unresolved",reason:="controlled native setup refusal",close:=function() transferRefusalCloses:=transferRefusalCloses+1; end); end;;
gap> transferRefusedFull:=koFull(transferLargeFull.ahss);;
gap> Assert(0,transferRefusedFull.degreeResults[5].status="unresolved" and transferRefusedFull.degreeResults[5].pendingLayer="model-setup");
gap> Assert(0,transferRefusedFull.degreeResults[5].reason="controlled native setup refusal" and transferRefusalCloses=1);
gap> Assert(0,not IsBound(transferRefusedFull.degreeResults[5].transferAttempt));
gap> # Degree six is unresolved even for zero layers, without a model request.
gap> transferZeroGroup:=AbelianGroup(IsPcpGroup,[]);;
gap> transferZeroCell:=rec(group:=transferZeroGroup,lift:=x->x);;
gap> transferZeroAHSS:=rec(kind:="koAHSSResult",computedThrough:=6,maxDegree:=6,pages:=rec(pageNumbers:=[6]),_context:=rec(getCell:=function(page,p,q) return transferZeroCell; end,backend:=rec(cohomologyData:=function(p,q) return rec(represent:=x->[]); end)));;
gap> transferSixFull:=koFull(transferZeroAHSS);;
gap> Assert(0,transferRefusalCloses=4);
gap> Assert(0,transferSixFull.degreeResults[8].status="unresolved" and transferSixFull.degreeResults[8].pendingLayer="model-setup");
gap> Assert(0,transferSixFull.degreeResults[8].reason="native extension degree six is not implemented");
gap> KOAHSS_ExtensionTransferredModel:=transferSavedFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionTransferredModel");
gap> KOAHSS_ExtensionBarModel:=transferSavedBarFactory;;
gap> MakeReadOnlyGlobal("KOAHSS_ExtensionBarModel");
