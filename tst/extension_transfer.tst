# Audit the sparse comparison and the opt-in nonlinear transferred model.
gap> for transferOrder in [2,4] do
> transferR:=ResolutionFiniteGroup(CyclicGroup(transferOrder),8);
> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS([1],0,7);
> transferAudit:=KOAHSS_ExtensionTransferPreflight(transferBackend,5);
> Assert(0,transferAudit.status="checked" and transferAudit.chainRetractionVerified);
> Assert(0,transferAudit.checkedThrough=7 and Length(transferAudit.degrees)=8);
> Assert(0,not transferAudit.transferReady and not transferAudit.sideConditionsVerified);
> Assert(0,not transferAudit.nonlinearIdentitiesVerified);
> transferTr:=transferBackend.naturalBar();
> for transferDegree in [0..7] do
>     Assert(0,transferAudit.degrees[transferDegree+1].basesChecked=transferR!.dimension(transferDegree));
>     for transferSign in [0,[1]] do
>         for transferBasis in [1..transferR!.dimension(transferDegree)] do
>             transferVector:=List([1..transferR!.dimension(transferDegree)],i->0);
>             transferVector[transferBasis]:=1;
>             Assert(0,transferTr.pullback(transferDegree,
>                 transferTr.lift(transferDegree,transferVector,transferSign),transferSign)=transferVector);
>         od;
>     od;
> od;
> od;
gap> Assert(0,KOAHSS_ExtensionTransferPreflight(transferBackend,3,rec(maxTerms:=1)).code="term-budget");
gap> Assert(0,KOAHSS_ExtensionTransferPreflight(transferBackend,3,rec(maxSupport:=1)).code="support-budget");
gap> Assert(0,KOAHSS_ExtensionTransferPreflight(rec(),3).code="missing-bar-transport");
gap> transferShortR:=ResolutionFiniteGroup(CyclicGroup(2),3);;
gap> transferShortBackend:=rec(naturalBar:=function() return KOAHSS_NaturalBarTransport(transferShortR); end);;
gap> Assert(0,KOAHSS_ExtensionTransferPreflight(transferShortBackend,3).code="resolution-length");
gap> # Add a contractible free ZG pair u1,u2 to a valid C2 resolution. R0 is
gap> # unchanged, d(u2)=u1, h(u1)=u2, and g(u1)=0. Thus rank-one R0 is insufficient.
gap> transferBase:=ResolutionFiniteGroup(CyclicGroup(2),6);;
gap> transferExtra:=Objectify(TypeObj(transferBase),rec(
> group:=transferBase!.group,elts:=transferBase!.elts,properties:=transferBase!.properties,
> dimension:=function(n)
>     if n in [1,2] then return transferBase!.dimension(n)+1; fi;
>     return transferBase!.dimension(n);
> end,
> boundary:=function(n,j)
>     if n=1 and j=transferBase!.dimension(1)+1 then return []; fi;
>     if n=2 and j=transferBase!.dimension(2)+1 then
>         return [[transferBase!.dimension(1)+1,Position(transferBase!.elts,One(transferBase!.group))]];
>     fi;
>     return transferBase!.boundary(n,j);
> end,
> homotopy:=function(n,t)
>     if n=1 and AbsInt(t[1])=transferBase!.dimension(1)+1 then
>         return [[SignInt(t[1])*(transferBase!.dimension(2)+1),t[2]]];
>     fi;
>     if n=2 and AbsInt(t[1])=transferBase!.dimension(2)+1 then return []; fi;
>     return transferBase!.homotopy(n,t);
> end));;
gap> transferExtraTr:=KOAHSS_NaturalBarTransport(transferExtra);;
gap> transferExtraBackend:=rec(naturalBar:=function() return transferExtraTr; end);;
gap> transferFailure:=KOAHSS_ExtensionTransferPreflight(transferExtraBackend,3);;
gap> Assert(0,transferExtra!.dimension(0)=1 and transferFailure.status="unresolved");
gap> Assert(0,transferFailure.code="chain-retraction-failed");
gap> Assert(0,transferFailure.failure.degree=1 and transferFailure.failure.basis=2);
gap> Assert(0,transferFailure.failure.image=[] and transferFailure.failure.expected=[[2,1,1]]);
gap> Assert(0,transferExtraTr.pullback(1,transferExtraTr.lift(1,[0,1],0),0)=[0,0]);
gap> Assert(0,not transferFailure.chainRetractionVerified and not transferFailure.transferReady);
gap> # The normalized homotopy satisfies the integral group-ring identities,
gap> # including nonzero homotopies on a proper C4 retract.
gap> transferR:=ResolutionFiniteGroup(CyclicGroup(4),7);;
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS([1],0,7);;
gap> transferEngine:=KOAHSS_ExtensionNormalizedTransport(transferBackend,3);;
gap> transferTr:=transferBackend.naturalBar();;
gap> transferReduceBar:=function(chain)
> local keys;
> keys:=Set(List(Filtered(chain,t->not ForAny([2..Length(t[2])],i->t[2][i]=t[2][i-1])),t->t[2]));
> return Filtered(List(keys,sigma->[Sum(Filtered(chain,t->t[2]=sigma),t->t[1]),sigma]),t->t[1]<>0);
> end;;
gap> transferBoundary:=function(chain)
> local result,t,i;
> result:=[];
> for t in chain do
>     if Length(t[2])>1 then
>         for i in [1..Length(t[2])] do
>             Add(result,[t[1]*(-1)^(i-1),t[2]{Filtered([1..Length(t[2])],j->j<>i)}]);
>         od;
>     fi;
> od;
> return transferReduceBar(result);
> end;;
gap> transferApplyH:=function(chain)
> local result,t,h;
> result:=[];
> for t in chain do
>     for h in transferEngine.homotopy(t[2]) do Add(result,[t[1]*h[1],h[2]]); od;
> od;
> return transferReduceBar(result);
> end;;
gap> transferApplyGF:=function(chain)
> local result,t,f,g;
> result:=[];
> for t in chain do
>     for f in transferTr.f(t[2]) do
>         for g in transferTr.g(Length(t[2])-1,f[1]) do
>             Add(result,[t[1]*f[3]*g[1],List(g[2],v->f[2]*v)]);
>         od;
>     od;
> od;
> return transferReduceBar(result);
> end;;
gap> transferApplyF:=function(chain)
> local result,t,f,keys;
> result:=[];
> for t in chain do
>     for f in transferTr.f(t[2]) do Add(result,[t[1]*f[3],[f[1],f[2]]]); od;
> od;
> keys:=Set(List(result,t->t[2]));
> return Filtered(List(keys,key->[Sum(Filtered(result,t->t[2]=key),t->t[1]),key]),t->t[1]<>0);
> end;;
gap> transferUnit:=One(transferR!.group);; transferNonunit:=Filtered(AsList(transferR!.group),g->g<>transferUnit);;
gap> transferNonzeroH:=false;;
gap> for transferDegree in [0..3] do
> for transferIncrements in Tuples(transferNonunit,transferDegree) do
>     transferSigma:=[transferUnit];
>     for transferIncrement in transferIncrements do Add(transferSigma,Last(transferSigma)*transferIncrement); od;
>     transferChain:=[[1,transferSigma]];
>     transferH:=transferApplyH(transferChain);
>     if not IsEmpty(transferH) then transferNonzeroH:=true; fi;
>     transferLhs:=transferReduceBar(Concatenation(transferBoundary(transferH),transferApplyH(transferBoundary(transferChain))));
>     transferRhs:=transferReduceBar(Concatenation(transferChain,List(transferApplyGF(transferChain),t->[-t[1],t[2]])));
>     Assert(0,transferLhs=transferRhs);
>     Assert(0,transferApplyF(transferH)=[] and transferApplyH(transferH)=[]);
> od;
> Assert(0,ForAll([1..transferR!.dimension(transferDegree)],j->transferApplyH(transferTr.g(transferDegree,j))=[]));
> od;
gap> Assert(0,transferNonzeroH);
gap> # C2 is an exact bar-basis isomorphism. Check its native nonlinear
gap> # product against the reference and retain the nonzero integral carry.
gap> transferR:=ResolutionFiniteGroup(CyclicGroup(2),7);;
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS(0,0,7);;
gap> transferModel:=KOAHSS_ExtensionTransferredModel(transferBackend,3);;
gap> transferBar:=KOAHSS_ExtensionBarModel(transferBackend,3);;
gap> Assert(0,transferModel.classEquivalenceVerified and transferModel.dimension(5)=1);
gap> transferState:=transferModel.zero(3);; transferState.C:=[1];;
gap> Assert(0,transferModel.d(3,transferState)=transferModel.zero(4));
gap> transferProduct:=transferModel.xtimes(3,transferState,transferState);;
gap> Assert(0,transferProduct=transferBar.xtimes(3,transferState,transferState));
gap> Assert(0,transferProduct.D=[1] and transferProduct.C=[0]);
gap> Assert(0,transferModel.barState(3,transferState,transferBar)=transferState);
gap> Assert(0,transferModel.divideLeft(3,transferState,transferProduct)=transferState);
gap> transferGauge:=transferModel.zero(2);; transferGauge.C:=[1];;
gap> Assert(0,transferModel.act(3,transferGauge,transferState)=transferBar.xtimes(3,transferBar.d(2,transferGauge),transferState));
gap> transferModel.close();; transferBar.close();;
gap> # Nonclosed input reports its first obstruction with correct-sized
gap> # zero padding, rather than trying to evaluate an illegal upper formula.
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS([1],0,7);;
gap> transferModel:=KOAHSS_ExtensionTransferredModel(transferBackend,3);;
gap> transferState:=transferModel.zero(3);; transferState.A:=[1];;
gap> transferCurvature:=transferModel.d(3,transferState);;
gap> Assert(0,transferCurvature.A=[-2] and transferCurvature.B=[0] and transferCurvature.C=[0] and transferCurvature.D=[0]);
gap> transferModel.close();;
gap> # A proper retract: compare a native C3 product with the full bar
gap> # through the actual embedding and its explicit reflection gauge.
gap> transferR:=ResolutionFiniteGroup(CyclicGroup(3),6);;
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS(0,0,6);;
gap> transferModel:=KOAHSS_ExtensionTransferredModel(transferBackend,3);;
gap> transferBar:=KOAHSS_ExtensionBarModel(transferBackend,3);;
gap> Assert(0,not transferModel.classEquivalenceVerified);
gap> transferState:=transferModel.zero(3);; transferState.C:=[1];;
gap> Assert(0,KOAHSS_ExtensionStateIsZero(transferModel.d(3,transferState)));
gap> transferExport:=transferModel.barState(3,transferState,transferBar);;
gap> Assert(0,KOAHSS_ExtensionStateIsZero(transferBar.d(3,transferExport)));
gap> transferProduct:=transferModel.xtimes(3,transferState,transferState);;
gap> Assert(0,transferProduct.C=[0] and transferProduct.D=[1]);
gap> transferProductExport:=transferModel.barState(3,transferProduct,transferBar);;
gap> transferWitness:=transferModel.debugRequest(rec(operation:="gauge_values",degree:=3,state:=transferState,other:=transferState));;
gap> Assert(0,transferWitness.state=transferProduct);
gap> Assert(0,transferBar.xtimes(3,transferExport,transferExport)=transferBar.xtimes(3,transferBar.d(2,transferWitness.gauge),transferProductExport));
gap> Assert(0,transferModel.divideLeft(3,transferState,transferProduct)=transferState);
gap> transferModel.close();; transferBar.close();;
gap> # The valid contractible-summand resolution is still accepted by the
gap> # public API; only the strict native retraction is refused, with a
gap> # recorded reference-model fallback in the original supplied basis.
gap> transferExtraFull:=koFull(transferExtra,0,0,3,rec(extensionModel:="transfer"));;
gap> Assert(0,IsIdenticalObj(transferExtraFull.ahss._context.resolution,transferExtra));
gap> Assert(0,transferExtraFull.invariants=[[0],[],[2,2],[2,2],[0,8]]);
gap> Assert(0,transferExtraFull.degreeResults[5].modelSelection.selected="bar");
gap> Assert(0,transferExtraFull.degreeResults[5].modelSelection.fallback);
gap> Assert(0,transferExtraFull.degreeResults[5].transferAttempt.reason=transferFailure.reason);
gap> Assert(0,transferExtraFull.degreeResults[5].transferAttempt.pendingLayer="model-setup");
gap> # A nonzero reflection gauge on a proper C4 retract checks both the
gap> # binary C correction and the integral D correction on the whole bar.
gap> transferR:=ResolutionFiniteGroup(CyclicGroup(4),6);;
gap> transferBackend:=koAHSSHAPSpace(transferR,koAHSSNaturalOperations()).koAHSS(0,0,6);;
gap> transferModel:=KOAHSS_ExtensionTransferredModel(transferBackend,3);;
gap> transferBar:=KOAHSS_ExtensionBarModel(transferBackend,3);;
gap> Assert(0,not transferModel.classEquivalenceVerified);
gap> transferState:=transferModel.zero(3);; transferState.B:=[1];;
gap> Assert(0,KOAHSS_ExtensionStateIsZero(transferModel.d(3,transferState)));
gap> transferProduct:=transferModel.xtimes(3,transferState,transferState);;
gap> Assert(0,transferProduct=rec(A:=[0],B:=[0],C:=[0],D:=[6]));
gap> transferWitness:=transferModel.debugRequest(rec(operation:="gauge_values",degree:=3,state:=transferState,other:=transferState));;
gap> Assert(0,transferWitness.state=transferProduct);
gap> Assert(0,transferWitness.gauge.C=[0,1,1] and ForAny(transferWitness.gauge.D,x->x<>0));
gap> Assert(0,not KOAHSS_ExtensionStateIsZero(transferWitness.gauge));
gap> transferExport:=transferModel.barState(3,transferState,transferBar);;
gap> transferProductExport:=transferModel.barState(3,transferProduct,transferBar);;
gap> Assert(0,KOAHSS_ExtensionStateIsZero(transferBar.d(3,transferExport)));
gap> Assert(0,KOAHSS_ExtensionStateIsZero(transferBar.d(3,transferProductExport)));
gap> Assert(0,transferBar.xtimes(3,transferExport,transferExport)=transferBar.xtimes(3,transferBar.d(2,transferWitness.gauge),transferProductExport));
gap> transferModel.close();; transferBar.close();;
