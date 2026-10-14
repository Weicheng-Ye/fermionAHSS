# Cochain equations and mod-two Adem primitives. Matrices act on row vectors.
# Descending pivots make the particular solution lexicographically least.
InstallGlobalFunction(koAHSSSolveMod2System, function(matrix, rhs)
    local variables, equations, zero, one, pivots, column, row, j, k,
          pivot, particular, generators, free, vector, value, solve;
    if not IsList(matrix) or not IsList(rhs) or not ForAll(rhs, IsInt) then
        Error("koAHSS: mod-two matrix and rhs must contain integers");
    fi;
    variables := Length(matrix); equations := Length(rhs);
    if not ForAll(matrix, r -> IsList(r) and Length(r) = equations
                                      and ForAll(r, IsInt)) then
        Error("koAHSS: mod-two matrix rows must have Length(rhs) integer entries");
    fi;
    zero := Zero(GF(2)); one := One(GF(2)); pivots := [];
    for column in [1..equations] do
        row := List([1..variables], i -> (matrix[i][column] mod 2)*one);
        Add(row, (rhs[column] mod 2)*one);
        ConvertToVectorRep(row, 2);
        pivot := 0;
        for j in Reversed([1..variables]) do
            if row[j] <> zero then
                if IsBound(pivots[j]) then AddRowVector(row, pivots[j]);
                else pivot := j; pivots[j] := row; break; fi;
            fi;
        od;
        if pivot = 0 and row[variables+1] <> zero then return fail; fi;
    od;
    solve := function(freeIndex, inhomogeneous)
        local solution, i, h, v;
        solution := List([1..variables], i -> 0);
        if freeIndex > 0 then solution[freeIndex] := 1; fi;
        for i in [1..variables] do
            if IsBound(pivots[i]) then
                v := 0;
                if inhomogeneous and pivots[i][variables+1] <> zero then v := 1; fi;
                for h in [1..i-1] do
                    if solution[h] = 1 and pivots[i][h] <> zero then v := 1-v; fi;
                od;
                solution[i] := v;
            fi;
        od;
        return solution;
    end;
    particular := solve(0, true);
    free := Filtered([1..variables], i -> not IsBound(pivots[i]));
    generators := List(free, i -> solve(i, false));
    return rec(particular := particular, homogeneousGenerators := generators,
               homogeneousOrders := List(free, i -> 2), modulus := 2,
               rank := variables-Length(free));
end);

BindGlobal("KOAHSS_SolveF2System", koAHSSSolveMod2System);

BindGlobal("KOAHSS_EquationBackend", function(backend, degree)
    if not IsRecord(backend) or not IsBound(backend.dimension)
       or not IsFunction(backend.dimension) or not IsBound(backend.coboundary)
       or not IsFunction(backend.coboundary) then
        Error("koAHSS: equation solver needs backend.dimension and backend.coboundary");
    fi;
    if not IsInt(degree) or degree < 0 then
        Error("koAHSS: cochain degree must be a nonnegative integer");
    fi;
end);

BindGlobal("KOAHSS_EquationDimension", function(backend, degree)
    local d;
    d := backend.dimension(degree);
    if not IsInt(d) or d < 0 then
        Error("koAHSS: cochain dimension must be a nonnegative integer");
    fi;
    return d;
end);

BindGlobal("KOAHSS_EquationDifferential", function(backend, degree, vector)
    local value;
    value := backend.coboundary(degree, vector, false);
    KOAHSS_CC_CheckVector(value, KOAHSS_EquationDimension(backend,degree+1),
                         "ordinary coboundary output");
    return value;
end);

InstallGlobalFunction(koAHSSSolveCochainEquation, function(arg)
    local backend, degree, rhs, dimension, width, matrix, basis, i,
          constraints, extra, result, image, vector;
    if not Length(arg) in [3,4] then
        Error("usage: koAHSSSolveCochainEquation(backend,degree,rhs[,constraints])");
    fi;
    backend := arg[1]; degree := arg[2]; rhs := arg[3];
    KOAHSS_EquationBackend(backend,degree);
    dimension := KOAHSS_EquationDimension(backend,degree);
    width := KOAHSS_EquationDimension(backend,degree+1);
    KOAHSS_CC_CheckVector(rhs,width,"cochain equation rhs");
    rhs := List(rhs,x -> x mod 2);
    matrix := [];
    for i in [1..dimension] do
        basis := List([1..dimension], j -> 0); basis[i] := 1;
        Add(matrix, KOAHSS_EquationDifferential(backend,degree,basis));
    od;
    extra := [];
    if Length(arg) = 4 then
        constraints := arg[4];
        if not IsRecord(constraints) or not IsBound(constraints.matrix)
           or not IsBound(constraints.rhs) or not IsList(constraints.matrix)
           or Length(constraints.matrix) <> dimension
           or not IsList(constraints.rhs) or not ForAll(constraints.rhs,IsInt)
           or not ForAll(constraints.matrix, row -> IsList(row)
                 and Length(row)=Length(constraints.rhs) and ForAll(row,IsInt)) then
            Error("koAHSS: constraints need an integer matrix with one row per unknown and an rhs");
        fi;
        extra := constraints.rhs;
        matrix := List([1..dimension],i -> Concatenation(matrix[i],constraints.matrix[i]));
    fi;
    result := koAHSSSolveMod2System(matrix,Concatenation(rhs,extra));
    if result = fail then return fail; fi;
    # Check the action too, so a malformed custom backend cannot pass merely
    # because its values on basis vectors happened to define a matrix.
    image := List(KOAHSS_EquationDifferential(backend,degree,result.particular),x -> x mod 2);
    if image <> rhs then Error("koAHSS: cochain equation verification failed"); fi;
    for vector in result.homogeneousGenerators do
        if ForAny(KOAHSS_EquationDifferential(backend,degree,vector),x -> x mod 2 <> 0) then
            Error("koAHSS: homogeneous cochain equation verification failed");
        fi;
    od;
    result.status := "solved";
    result.primitive := ShallowCopy(result.particular);
    result.primitiveDegree := degree;
    result.rhs := rhs;
    result.constraintsApplied := Length(extra);
    result.choice := "lexicographically least in the supplied cochain basis";
    return result;
end);


# Exact integral affine systems, x * matrix = rhs.  Unlike phase equations,
# witnesses for integral coboundaries must not be reduced modulo a modulus.
# The same differential occurs with many right-hand sides in lift and gauge
# searches. Cache only its exact Smith preparation, never a solved affine
# family. Keys are immutable value snapshots, so caller mutation cannot make
# a preparation apply to a different matrix. The column count distinguishes
# empty matrices with different target ranks.
BindGlobal("KOAHSS_PrepareIntegerSystem",CallFuncList(function()
    local cache,cells,maxEntries,maxCells;
    cache:=[]; cells:=0; maxEntries:=8; maxCells:=2000000;
    return function(matrix,cols)
        local rows,index,prepared,snapshot,snf,size,old;
        rows:=Length(matrix);
        index:=PositionProperty(cache,item->item.cols=cols and item.matrix=matrix);
        if index<>fail then
            prepared:=Remove(cache,index); Add(cache,prepared); return prepared;
        fi;
        snapshot:=StructuralCopy(matrix);
        if rows=0 or cols=0 then
            snf:=fail; size:=rows;
        else
            snf:=SmithNormalFormIntegerMatTransforms(snapshot);
            # Discard the implementation's auxiliary transforms. These are
            # the only fields used by the affine solver below.
            snf:=rec(rank:=snf.rank,normal:=snf.normal,
                rowtrans:=snf.rowtrans,coltrans:=snf.coltrans);
            size:=2*rows*cols+rows^2+cols^2;
        fi;
        prepared:=rec(matrix:=snapshot,rows:=rows,cols:=cols,snf:=snf,cells:=size);
        MakeImmutable(prepared);
        if size<=maxCells then
            while Length(cache)>=maxEntries or cells+size>maxCells do
                old:=Remove(cache,1); cells:=cells-old.cells;
            od;
            Add(cache,prepared); cells:=cells+size;
        fi;
        return prepared;
    end;
end,[]));

InstallGlobalFunction(koAHSSSolveIntegerSystem, function(matrix,rhs)
    local rows, cols, prepared, snf, transformed, y, i, generators, solution;
    if not IsList(matrix) or not IsList(rhs) or not ForAll(rhs,IsInt) then
        Error("koAHSS: integer equation matrix and rhs must contain integers");
    fi;
    rows:=Length(matrix); cols:=Length(rhs);
    if not ForAll(matrix,r -> IsList(r) and Length(r)=cols and ForAll(r,IsInt)) then
        Error("koAHSS: integer equation matrix has incorrect dimensions");
    fi;
    prepared:=KOAHSS_PrepareIntegerSystem(matrix,cols);
    if rows=0 then
        if ForAny(rhs,x -> x<>0) then return fail; fi;
        return rec(particular:=[],homogeneousGenerators:=[],rank:=0,coefficientRing:=Integers);
    fi;
    if cols=0 then
        return rec(particular:=List([1..rows],i->0),homogeneousGenerators:=IdentityMat(rows),
                   rank:=0,coefficientRing:=Integers);
    fi;
    snf:=prepared.snf;
    transformed:=rhs*snf.coltrans;
    y:=List([1..rows],i->0);
    for i in [1..snf.rank] do
        if transformed[i] mod snf.normal[i][i]<>0 then return fail; fi;
        y[i]:=QuoInt(transformed[i],snf.normal[i][i]);
    od;
    if ForAny([snf.rank+1..cols],i->transformed[i]<>0) then return fail; fi;
    solution:=y*snf.rowtrans;
    generators:=List([snf.rank+1..rows],i->ShallowCopy(snf.rowtrans[i]));
    if solution*matrix<>rhs or ForAny(generators,v->ForAny(v*matrix,x->x<>0)) then
        Error("koAHSS: integer equation verification failed");
    fi;
    return rec(particular:=solution,homogeneousGenerators:=generators,
               rank:=snf.rank,coefficientRing:=Integers);
end);
