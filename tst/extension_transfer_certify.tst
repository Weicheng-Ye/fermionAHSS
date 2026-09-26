# A native candidate is accepted by the same literal bar checks as the old engine.
gap> certificateR:=ResolutionFiniteGroup(CyclicGroup(2),6);;
gap> certificateBackend:=koAHSSHAPSpace(certificateR,koAHSSNaturalOperations()).koAHSS(0,0,6);;
gap> certificateModel:=rec(modelId:="identity-resolution-test",classEquivalenceVerified:=false,barState:=function(k,state,bar) return state; end);;
gap> certificateLayers:=rec();;
gap> for certificateName in ["D","C","B","A"] do
> certificateLayers.(certificateName):=rec(orders:=[],fullLifts:=[]);
> od;
gap> certificateLayers.D.orders:=[2];;
gap> certificateLayers.D.p:=4;; certificateLayers.D.cochains:=[[1]];;
gap> certificateLayers.D.fullLifts:=[rec(status:="computed",state:=rec(A:=[0],B:=[0],C:=[0],D:=[1]))];;
gap> certificatePresentation:=koAHSSExtensionFromLayers(rec(D:=[2]),fail);;
gap> certificateProof:=KOAHSS_ExtensionTransferCertify(certificateBackend,certificateModel,3,certificateLayers,certificatePresentation);;
gap> Assert(0,certificateProof.status="computed" and certificateProof.certificateLevel="complete-bar");
gap> Assert(0,not certificateProof.gaugeCompletenessAssumed);
gap> Assert(0,Length(certificateProof.relationWitnesses)=1);
gap> Assert(0,certificateProof.relationWitnesses[1].canonicalComparison.equalityVerified);
gap> Assert(0,certificateProof.relationWitnesses[1].stackedState.D=[2]);
gap> Assert(0,certificateProof.algebraAudit.associativityVerified and certificateProof.algebraAudit.commutativityVerified);
gap> Assert(0,certificateProof.barLayers.D.fullLifts[1].state.D=[1]);
gap> Assert(0,KOAHSS_ExtensionTransferCertify(certificateBackend,rec(),3,certificateLayers,certificatePresentation).status="unresolved");
gap> certificateLargeR:=ResolutionFiniteGroup(CyclicGroup(6),7);;
gap> certificateLargeBackend:=koAHSSHAPSpace(certificateLargeR,koAHSSNaturalOperations()).koAHSS(0,0,7);;
gap> certificateUnavailable:=KOAHSS_ExtensionTransferCertify(certificateLargeBackend,certificateModel,4,certificateLayers,certificatePresentation);;
gap> Assert(0,certificateUnavailable.status="unresolved" and certificateUnavailable.certificateLevel="complete-bar");
gap> certificateBadModel:=ShallowCopy(certificateModel);;
gap> certificateBadModel.barState:=function(k,state,bar) return bar.zero(k); end;;
gap> KOAHSS_ExtensionTransferCertify(certificateBackend,certificateBadModel,3,certificateLayers,certificatePresentation);
Error, koFull: transfer export changed its marked E6 leading representative
Error, koFull: independent bar certification failed; the original error is reported above
