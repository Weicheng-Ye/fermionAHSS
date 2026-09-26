# This audits the sparse comparison; it does not enable a nonlinear model.
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
