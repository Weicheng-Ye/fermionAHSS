# Localization of the relation measurement at the primes: the prime of a
# generator decides its model (complete, three-local or split), and the
# odd-primary relations of degree five are measured without the complete
# four-layer formulas.
gap> Assert(0,KOAHSS_RelationPrime(9)=3 and KOAHSS_RelationPrime(8)=2 and KOAHSS_RelationPrime(2)=2);
gap> Assert(0,KOAHSS_RelationPrime(6)=fail and KOAHSS_RelationPrime(0)=fail and KOAHSS_RelationPrime(1)=fail);
gap> Assert(0,KOAHSS_ExtensionRelationModel(5,"A",3).model="three-local");
gap> Assert(0,KOAHSS_ExtensionRelationModel(5,"A",9).model="three-local");
gap> Assert(0,KOAHSS_ExtensionRelationModel(6,"A",3).model="three-local");
gap> Assert(0,KOAHSS_ExtensionRelationModel(4,"A",3).model="split");
gap> Assert(0,KOAHSS_ExtensionRelationModel(5,"A",5).model="split" and KOAHSS_ExtensionRelationModel(6,"A",25).model="split");
gap> Assert(0,KOAHSS_ExtensionRelationModel(5,"A",4).model="complete" and KOAHSS_ExtensionRelationModel(5,"A",0).model="complete");
gap> Assert(0,KOAHSS_ExtensionRelationModel(5,"D",3).model="complete");
gap> Assert(0,IsBool(KOAHSS_PrimeLocalizationEnabled()));
gap> # The generators are grouped by prime right after the layers are read;
gap> # free generators belong to no prime and are lower generators of all.
gap> splitParts := KOAHSS_ExtensionPrimeParts(rec(D:=rec(orders:=[5,2,0]),C:=rec(orders:=[2]),B:=rec(orders:=[]),A:=rec(orders:=[5,2,4,3,0])));;
gap> Assert(0,List(splitParts.parts,part->part.prime)=[2,3,5] and splitParts.names=["D","C","B","A"]);
gap> Assert(0,splitParts.parts[1].generators=rec(D:=[2],C:=[1],B:=[],A:=[2,3]));
gap> Assert(0,splitParts.parts[3].generators=rec(D:=[1],C:=[],B:=[],A:=[1]));
gap> Assert(0,splitParts.free=rec(D:=[3],C:=[],B:=[],A:=[5]));
gap> Assert(0,KOAHSS_ExtensionPrimeParts(rec(D:=rec(orders:=[2]),C:=rec(status:="unresolved"),B:=rec(orders:=[2]),A:=rec(orders:=[3]))).names=["D"]);
gap> # A relation that splits at its prime needs no stacking model.
gap> splitLayers := rec(A:=[5],D:=[5]);;
gap> splitResult := koAHSSExtensionFromLayers(splitLayers,function(layer,index,order,lower)
>     return KOAHSS_ExtensionSplitResponse(order,5,lower); end);;
gap> Assert(0,splitResult.invariants=[5,5] and Last(splitResult.extensionVectors).result.witness.model="split");
gap> Assert(0,KOAHSS_ExtensionPrimeSummary(splitResult)=[rec(prime:=5,relations:=1,models:=["split"])]);
gap> # Z/3, untwisted, degree five: A=Z/3 at (2,0) and D=Z/3 at (6,-4); the
gap> # three-local relation 3a = 2d (mod 3) gives Z/9, as ko_(3)^2(BZ/3) requires.
gap> z3 := koFull(CyclicGroup(3),0,0,5);;
gap> Assert(0,z3.invariants=[9]);
gap> z3A := First(z3.degreeResult.extensionVectors,v->v.layer="A");;
gap> Assert(0,z3A.result.witness.model="three-local" and z3A.result.witness.prime=3);
gap> Assert(0,z3A.result.witness.measuredLayers=["D"] and z3A.result.lowerCoordinates[1] mod 3<>0);
gap> Assert(0,z3.degreeResult.primes=[rec(prime:=3,relations:=1,models:=["three-local"])]);
gap> Assert(0,z3.degreeResult.layers.A.fullLifts[1].model="three-local");
gap> # Z/5: the rows q=0 and q=-4 lie in different Adams summands of ko_(5); the
gap> # degree splits without a stacking model.
gap> z5 := koFull(CyclicGroup(5),0,0,5);;
gap> Assert(0,z5.invariants=[5,5] and z5.degreeResult.certificateLevel="prime-split");
gap> measured := result -> Filtered(result.degreeResult.extensionVectors,v->not IsEmpty(v.lowerGeneratorIds));;
gap> Assert(0,Length(measured(z5))=1 and ForAll(measured(z5),v->v.result.witness.model="split"));
gap> Assert(0,z5.degreeResult.primes=[rec(prime:=5,relations:=1,models:=["split"])]);
gap> # Z/9: 9a = 6d modulo 9, so Z/3 + Z/27, the truncation of ko^2(BZ/9)_(3)
gap> # computed independently from the representation ring of Z/9.
gap> z9 := koFull(CyclicGroup(9),0,0,5);;
gap> Assert(0,z9.invariants=[3,27]);
gap> # Z/6 = Z/2 x Z/3: the two-primary part of its E6 line is that of Z/2
gap> # (here zero), the three-primary relation is measured three-locally.
gap> z6 := koFull(CyclicGroup(6),0,0,5);;
gap> c2 := koFull(CyclicGroup(2),0,0,5);;
gap> Assert(0,SortedList(z6.invariants)=SortedList(Concatenation(c2.invariants,[9])));
gap> Assert(0,ForAll(measured(z6),v->v.result.witness.model in ["complete","three-local"]));
gap> Assert(0,First(measured(z6),v->v.order=3).result.witness.model="three-local");
gap> Assert(0,List(z6.degreeResult.primeParts,part->part.prime)=[3]);
gap> Assert(0,z6.degreeResult.primeParts[1].invariants=[9]);
gap> # A row keeps only its coordinates on the generators of its prime and the
gap> # free generators; here every measurement returns all ones.
gap> rowLayers := rec(D:=rec(orders:=[2,3]),C:=rec(orders:=[2]),B:=rec(orders:=[]),A:=rec(orders:=[3,0]));;
gap> rowEngine := rec(answer:=function(layer,index,order,lower)
>     return rec(status:="computed",lowerPresentationId:=lower.presentationId,
>         lowerCoordinates:=List([1..lower.generatorCount],j->1),witness:=rec(model:="stub"));
> end);;
gap> rows := KOAHSS_ExtensionPrimeRows(rowEngine,rowLayers,KOAHSS_ExtensionPrimeParts(rowLayers),5);;
gap> Assert(0,rows.C[1].lowerCoordinates=[1,0] and rows.C[1].witness.fullLowerCoordinates=[1,1]);
gap> Assert(0,rows.A[1].lowerCoordinates=[0,1,0] and not IsBound(rows.A[2]));
gap> Assert(0,koAHSSExtensionFromLayers(rowLayers,KOAHSS_ExtensionReplayOracle(rows)).relationMatrix
>     =[[2,0,0,0,0],[0,3,0,0,0],[-1,0,2,0,0],[0,-1,0,3,0]]);
gap> # Over a lower group with no generator of its prime and none free, a
gap> # relation needs no measurement.
gap> rowLayers := rec(D:=rec(orders:=[3]),C:=rec(orders:=[2]),B:=rec(orders:=[]),A:=rec(orders:=[]));;
gap> rows := KOAHSS_ExtensionPrimeRows(rec(answer:=function(arg) Error("measured"); end),rowLayers,KOAHSS_ExtensionPrimeParts(rowLayers),5);;
gap> Assert(0,rows.C[1].witness.kind="zero-local-lower-group" and rows.C[1].lowerCoordinates=[0]);
