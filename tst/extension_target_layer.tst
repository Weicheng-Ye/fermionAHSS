# The target layer of a relation: the lowest layer of the recorded lower
# presentation with a generator outside n*G_low. Generators below it cannot
# change the extension class, and the certificate exhibits their multiples.
gap> # Pin-: lower group <b,c,d | 2b=c, 2c=d, 2d=0> = Z/8 (layers D, C, B added in this order).
gap> targetLowerPin := rec(generatorCount:=3,generatorIds:=["D:1","C:1","B:1"],relationMatrix:=[[2,0,0],[-1,2,0],[0,-1,2]],layers:=[rec(name:="D",startColumn:=1,orders:=[2]),rec(name:="C",startColumn:=2,orders:=[2]),rec(name:="B",startColumn:=3,orders:=[2])]);;
gap> targetPin := KOAHSS_ExtensionTargetLayer(targetLowerPin,2);;
gap> Assert(0,targetPin.layer="B" and targetPin.index=1 and targetPin.freeGenerators=["B:1"]);
gap> Assert(0,List(targetPin.certificate,c->c.generatorId)=["D:1","C:1"]);
gap> Assert(0,ForAll(targetPin.certificate,c->2*c.combination+c.relationCoefficients*targetLowerPin.relationMatrix=List([1..3],function(i) if i=c.column then return 1; fi; return 0; end)));
gap> # Order four on the same lower group: 4*G_low = <d>, so c is free as well and the target is C.
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerPin,4).layer="C");
gap> # Z/2 + Z/4 with 2b=0, 2c=d: b and c are free, d=2c is not; the target is C, not B.
gap> targetLowerMixed := rec(generatorCount:=3,generatorIds:=["D:1","C:1","B:1"],relationMatrix:=[[2,0,0],[-1,2,0],[0,0,2]],layers:=targetLowerPin.layers);;
gap> targetMixed := KOAHSS_ExtensionTargetLayer(targetLowerMixed,2);;
gap> Assert(0,targetMixed.layer="C" and targetMixed.index=2 and targetMixed.freeGenerators=["C:1","B:1"]);
gap> Assert(0,List(targetMixed.certificate,c->c.generatorId)=["D:1"]);
gap> # A split D layer (Z/2)^3: every generator is free and the complete measurement is needed.
gap> targetLowerSplit := rec(generatorCount:=3,generatorIds:=["D:1","D:2","D:3"],relationMatrix:=[[2,0,0],[0,2,0],[0,0,2]],layers:=[rec(name:="D",startColumn:=1,orders:=[2,2,2])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerSplit,2).layer="D");
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerSplit,2).index=3);
gap> # An odd lower group and an even order: nothing is free, the class vanishes.
gap> targetLowerOdd := rec(generatorCount:=1,generatorIds:=["D:1"],relationMatrix:=[[3]],layers:=[rec(name:="D",startColumn:=1,orders:=[3])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerOdd,2).layer=fail);
gap> Assert(0,Length(KOAHSS_ExtensionTargetLayer(targetLowerOdd,2).certificate)=1);
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerOdd,3).layer="D");
gap> # A free lower generator of infinite order is free for every order.
gap> targetLowerFree := rec(generatorCount:=1,generatorIds:=["B:1"],relationMatrix:=[],layers:=[rec(name:="B",startColumn:=1,orders:=[0])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerFree,2).layer="B");
gap> # No lower generators: no target.
gap> Assert(0,KOAHSS_ExtensionTargetLayer(rec(generatorCount:=0,generatorIds:=[],relationMatrix:=[],layers:=[]),2).layer=fail);
gap> Assert(0,KOAHSS_ExtensionStateIsZeroThrough(rec(A:=[0],B:=[0],C:=[1],D:=[3]),1));
gap> Assert(0,not KOAHSS_ExtensionStateIsZeroThrough(rec(A:=[0],B:=[0],C:=[1],D:=[3]),2));
gap> Assert(0,KOAHSS_ExtensionStatesAgreeThrough(rec(A:=[1],B:=[1],C:=[0],D:=[5]),rec(A:=[1],B:=[1],C:=[1],D:=[0]),1));
gap> Assert(0,IsBool(KOAHSS_LayeredRelationsEnabled()));
gap> # A row measured only through its target layer holds for a shifted lift of
gap> # its generator. All layers Z/2 with the true relations 2c=d, 2b=d, 2a=b of
gap> # the stored lifts: the B row has target C and drops d=2c, so it holds for
gap> # b-c. Measured against the stored lift of b, the A row 2a=b (target C) then
gap> # gives Z/4+Z/4; recording the rows measures the B row through D, keeps the
gap> # A row, restates its certificate, and the group is Z/2+Z/8.
gap> shiftLayers := rec(D:=rec(name:="D",orders:=[2]),C:=rec(name:="C",orders:=[2]),B:=rec(name:="B",orders:=[2]),A:=rec(name:="A",orders:=[2]));;
gap> shiftTruth := rec(C:=[1],B:=[1,0],A:=[0,0,1]);;
gap> shiftCalls := [];;
gap> shiftEngine := rec(answer:=function(arg)
>     local layer,order,lower,fields,target,coordinates,stored,j,below;
>     layer:=arg[1]; order:=arg[3]; lower:=arg[4]; fields:=["A","B","C","D"];
>     coordinates:=ShallowCopy(shiftTruth.(layer.name));
>     if Length(arg)>=5 and arg[5].complete then
>         Add(shiftCalls,[layer.name,"complete"]);
>         return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>             lowerCoordinates:=coordinates,witness:=rec(model:="complete",truncatedBelow:=fail,
>             measuredLayers:=fields{[Position(fields,layer.name)+1..4]}));
>     fi;
>     Add(shiftCalls,[layer.name,"target"]);
>     target:=KOAHSS_ExtensionTargetLayer(lower,order);
>     if target.layer=fail then
>         return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>             lowerCoordinates:=0*coordinates,witness:=rec(model:="complete",measuredLayers:=[],
>             truncatedBelow:=fields[Position(fields,layer.name)+1],sufficiency:=target));
>     fi;
>     for stored in lower.layers do
>         if Position(fields,stored.name)>Position(fields,target.layer) then
>             for j in [1..Length(stored.orders)] do coordinates[stored.startColumn+j-1]:=0; od;
>         fi;
>     od;
>     below:=fail; if target.layer<>"D" then below:=fields[Position(fields,target.layer)+1]; fi;
>     return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>         lowerCoordinates:=coordinates,witness:=rec(model:="complete",truncatedBelow:=below,
>         measuredLayers:=fields{[Position(fields,layer.name)+1..Position(fields,target.layer)]},
>         sufficiency:=target));
> end);;
gap> Assert(0,koAHSSExtensionFromLayers(shiftLayers,shiftEngine.answer).invariants=[4,4]);
gap> shiftCalls := [];;
gap> shiftRows := KOAHSS_ExtensionPrimeRows(shiftEngine,shiftLayers,KOAHSS_ExtensionPrimeParts(shiftLayers),4);;
gap> Assert(0,shiftCalls=[["C","target"],["B","target"],["A","target"],["B","complete"]]);
gap> Assert(0,shiftRows.B[1].lowerCoordinates=[1,0] and IsBound(shiftRows.B[1].witness.measuredThroughD));
gap> Assert(0,shiftRows.A[1].lowerCoordinates=[0,0,1]);
gap> Assert(0,koAHSSExtensionFromLayers(shiftLayers,KOAHSS_ExtensionReplayOracle(shiftRows)).invariants=[2,8]);
gap> shiftLower := KOAHSS_ExtensionLowerPresentation(shiftLayers,["D","C","B"],KOAHSS_ExtensionReplayOracle(shiftRows));;
gap> shiftCertificate := shiftRows.A[1].witness.sufficiency.certificate;;
gap> Assert(0,List(shiftCertificate,c->c.generatorId)=["D:1"] and ForAll(shiftCertificate,c->2*c.combination
>     +c.relationCoefficients*shiftLower.relationMatrix=List([1..3],function(i) if i=c.column then return 1; fi; return 0; end)));
gap> # A later relation whose target is at or above the layer of the shifted
gap> # generator keeps its row: with 2c=d, 2b=c+d, 2a=b the B row drops d and
gap> # holds for b-c, the lower group of a is Z/8 generated by b, the target
gap> # of the A relation is B, and C and D lie in 2*Z/8. The group is Z/16.
gap> shiftTruth := rec(C:=[1],B:=[1,1],A:=[0,0,1]);; shiftCalls := [];;
gap> shiftRows := KOAHSS_ExtensionPrimeRows(shiftEngine,shiftLayers,KOAHSS_ExtensionPrimeParts(shiftLayers),4);;
gap> Assert(0,shiftCalls=[["C","target"],["B","target"],["A","target"]]);
gap> Assert(0,shiftRows.B[1].lowerCoordinates=[0,1] and shiftRows.A[1].witness.truncatedBelow="C");
gap> Assert(0,koAHSSExtensionFromLayers(shiftLayers,KOAHSS_ExtensionReplayOracle(shiftRows)).invariants=[16]);
gap> # When the measurement through D is unresolved, so is the later relation,
gap> # and the recorded row is kept.
gap> shiftTruth := rec(C:=[1],B:=[1,0],A:=[0,0,1]);;
gap> shiftFailing := rec(answer:=function(arg)
>     if Length(arg)>=5 then return rec(status:="unresolved",reason:="no complete lift"); fi;
>     return CallFuncList(shiftEngine.answer,arg);
> end);;
gap> shiftRows := KOAHSS_ExtensionPrimeRows(shiftFailing,shiftLayers,KOAHSS_ExtensionPrimeParts(shiftLayers),4);;
gap> Assert(0,shiftRows.A[1].status="unresolved" and shiftRows.A[1].pendingRelation.reason="no complete lift");
gap> Assert(0,shiftRows.B[1].lowerCoordinates=[0,0] and shiftRows.B[1].witness.truncatedBelow="D");
gap> Assert(0,koAHSSExtensionFromLayers(shiftLayers,KOAHSS_ExtensionReplayOracle(shiftRows)).status="unresolved");
gap> # A primary-operation row of a C generator over D is determined modulo 2H,
gap> # the relation of a lift shifted within D. A later relation measured
gap> # through D with a coefficient on that generator measures it through D.
gap> primaryLayers := rec(D:=rec(name:="D",orders:=[2,2]),C:=rec(name:="C",orders:=[2]),
>     B:=rec(name:="B",orders:=[2]),A:=rec(name:="A",orders:=[]));;
gap> primaryCalls := [];;
gap> primaryEngine := rec(answer:=function(arg)
>     local layer,lower,complete,witness;
>     layer:=arg[1]; lower:=arg[4]; complete:=Length(arg)>=5 and arg[5].complete;
>     Add(primaryCalls,[layer.name,complete]);
>     witness:=rec(model:="complete",truncatedBelow:=fail,measuredLayers:=["D"]);
>     if layer.name="C" then
>         if not complete then
>             witness:=rec(model:="primary-R",truncatedBelow:=fail,shiftedLift:="D",measuredLayers:=["D"]);
>         fi;
>         return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>             lowerCoordinates:=[1,0],witness:=witness);
>     fi;
>     return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>         lowerCoordinates:=[0,0,1],witness:=witness);
> end);;
gap> primaryRows := KOAHSS_ExtensionPrimeRows(primaryEngine,primaryLayers,KOAHSS_ExtensionPrimeParts(primaryLayers),4);;
gap> Assert(0,primaryCalls=[["C",false],["B",false],["C",true]]);
gap> Assert(0,IsBound(primaryRows.C[1].witness.measuredThroughD) and not IsBound(primaryRows.C[1].witness.shiftedLift));
