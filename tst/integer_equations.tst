# Reusing a Smith preparation preserves the original exact affine family.
gap> uncachedIntegerSolve := function(matrix,rhs)
> local snf,transformed,y,i;
> snf:=SmithNormalFormIntegerMatTransforms(matrix);
> transformed:=rhs*snf.coltrans; y:=List(matrix,row->0);
> for i in [1..snf.rank] do
>     if transformed[i] mod snf.normal[i][i]<>0 then return fail; fi;
>     y[i]:=QuoInt(transformed[i],snf.normal[i][i]);
> od;
> if ForAny([snf.rank+1..Length(rhs)],i->transformed[i]<>0) then return fail; fi;
> return rec(particular:=y*snf.rowtrans,
>     homogeneousGenerators:=List([snf.rank+1..Length(matrix)],i->ShallowCopy(snf.rowtrans[i])),
>     rank:=snf.rank,coefficientRing:=Integers);
> end;;
gap> integerMatrices := [[[2,4,0],[0,6,3],[4,8,0]],[[2,0,0],[0,4,0]],[[0,0],[0,0]],[[10^40,-2*10^40],[3,-6]],[[2,-4],[6,8],[0,10]]];;
gap> for integerMatrix in integerMatrices do
> integerPrepared:=KOAHSS_PrepareIntegerSystem(integerMatrix,Length(integerMatrix[1]));
> Assert(0,not IsMutable(integerPrepared) and not IsMutable(integerPrepared.matrix));
> Assert(0,IsIdenticalObj(integerPrepared,KOAHSS_PrepareIntegerSystem(StructuralCopy(integerMatrix),Length(integerMatrix[1]))));
> for integerCoefficients in Cartesian(List(integerMatrix,row->[-2..2])) do
>     integerRhs:=integerCoefficients*integerMatrix;
>     integerSolution:=koAHSSSolveIntegerSystem(integerMatrix,integerRhs);
>     Assert(0,integerSolution=uncachedIntegerSolve(integerMatrix,integerRhs));
>     Assert(0,integerSolution.particular*integerMatrix=integerRhs);
>     Assert(0,ForAll(integerSolution.homogeneousGenerators,v->ForAll(v*integerMatrix,x->x=0)));
> od;
> od;
gap> Assert(0,koAHSSSolveIntegerSystem([[2,0,0],[0,4,0]],[1,0,0])=fail);
gap> Assert(0,koAHSSSolveIntegerSystem([[2,0,0],[0,4,0]],[0,0,1])=fail);
gap> Assert(0,koAHSSSolveIntegerSystem([[0,0],[0,0]],[0,1])=fail);
gap> # Mutating caller inputs and returned families cannot alter a cached matrix.
gap> integerMatrix := [[2,4],[4,8]];; integerOriginal:=StructuralCopy(integerMatrix);;
gap> integerPrepared := KOAHSS_PrepareIntegerSystem(integerMatrix,2);;
gap> integerSolution := koAHSSSolveIntegerSystem(integerMatrix,[6,12]);;
gap> integerSolution.particular[1]:=999;; integerSolution.homogeneousGenerators[1][1]:=999;;
gap> Assert(0,koAHSSSolveIntegerSystem(integerMatrix,[6,12])=uncachedIntegerSolve(integerMatrix,[6,12]));
gap> integerMatrix[1][1]:=3;;
gap> Assert(0,integerPrepared.matrix=integerOriginal and IsMutable(integerMatrix) and IsMutable(integerMatrix[1]));
gap> Assert(0,not IsIdenticalObj(integerPrepared,KOAHSS_PrepareIntegerSystem(integerMatrix,2)));
gap> Assert(0,koAHSSSolveIntegerSystem(integerMatrix,[7,12])=uncachedIntegerSolve(integerMatrix,[7,12]));
gap> integerMatrix[1][1]:=2;;
gap> Assert(0,IsIdenticalObj(integerPrepared,KOAHSS_PrepareIntegerSystem(integerMatrix,2)));
gap> # Empty source/target modules retain their distinct exact dimensions.
gap> Assert(0,koAHSSSolveIntegerSystem([],[]).particular=[]);
gap> Assert(0,koAHSSSolveIntegerSystem([],[0,0]).homogeneousGenerators=[]);
gap> Assert(0,koAHSSSolveIntegerSystem([],[0,1])=fail);
gap> Assert(0,koAHSSSolveIntegerSystem([[],[],[]],[]).homogeneousGenerators=IdentityMat(3));
gap> Assert(0,koAHSSSolveIntegerSystem([[],[],[]],[]).particular=[0,0,0]);
gap> Assert(0,KOAHSS_PrepareIntegerSystem([],0).cols=0 and KOAHSS_PrepareIntegerSystem([],2).cols=2);
gap> Assert(0,not IsIdenticalObj(KOAHSS_PrepareIntegerSystem([],0),KOAHSS_PrepareIntegerSystem([],2)));
gap> # Preparations are bounded: old entries are evicted without changing answers.
gap> integerPrepared:=KOAHSS_PrepareIntegerSystem([[101]],1);;
gap> for integerValue in [102..109] do KOAHSS_PrepareIntegerSystem([[integerValue]],1); od;
gap> Assert(0,not IsIdenticalObj(integerPrepared,KOAHSS_PrepareIntegerSystem([[101]],1)));
gap> Assert(0,koAHSSSolveIntegerSystem([[101]],[303]).particular=[3]);
