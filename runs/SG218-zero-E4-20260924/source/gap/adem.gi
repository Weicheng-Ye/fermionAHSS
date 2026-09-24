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

BindGlobal("KOAHSS_AdemCocycle", function(backend, degree, a)
    KOAHSS_EquationBackend(backend,degree);
    KOAHSS_CC_CheckVector(a,KOAHSS_EquationDimension(backend,degree),"Adem input");
    a := List(a,x -> x mod 2);
    if ForAny(KOAHSS_EquationDifferential(backend,degree,a),x -> x mod 2 <> 0) then
        Error("koAHSS: Steenrod and Adem inputs must be mod-two cocycles");
    fi;
    return a;
end);

InstallGlobalFunction(koAHSSSq, function(backend, degree, a, k)
    local value;
    a := KOAHSS_AdemCocycle(backend,degree,a);
    if not IsInt(k) or k < 0 then Error("koAHSS: square index must be nonnegative"); fi;
    if k = 0 then return a; fi;
    if k > degree then
        return List([1..KOAHSS_EquationDimension(backend,degree+k)],i -> 0);
    fi;
    if not IsBound(backend.cupMod2) or not IsFunction(backend.cupMod2) then
        Error("koAHSS: Steenrod squares require backend.cupMod2");
    fi;
    value := backend.cupMod2(degree-k,degree,a,degree,a);
    KOAHSS_CC_CheckVector(value,KOAHSS_EquationDimension(backend,degree+k),"square output");
    value := List(value,x -> x mod 2);
    if ForAny(KOAHSS_EquationDifferential(backend,degree+k,value),x -> x mod 2 <> 0) then
        Error("koAHSS: the supplied cups did not produce a Steenrod cocycle");
    fi;
    return value;
end);

InstallGlobalFunction(koAHSSSolveAdemRelation, function(arg)
    local backend, degree, a, r, s, total, relation, terms, t, b, value,
          equationArgs, result;
    if not Length(arg) in [5,6] then
        Error("usage: koAHSSSolveAdemRelation(backend,degree,a,r,s[,constraints])");
    fi;
    backend := arg[1]; degree := arg[2]; a := arg[3]; r := arg[4]; s := arg[5];
    a := KOAHSS_AdemCocycle(backend,degree,a);
    if not IsInt(r) or not IsInt(s) or r < 1 or s < 1 or r >= 2*s then
        Error("koAHSS: an Adem pair requires positive r,s with r < 2*s");
    fi;
    total := degree+r+s;
    b := koAHSSSq(backend,degree,a,s);
    relation := koAHSSSq(backend,degree+s,b,r);
    terms := [];
    for t in [0..QuoInt(r,2)] do
        if Binomial(s-t-1,r-2*t) mod 2 = 1 then
            b := koAHSSSq(backend,degree,a,t);
            value := koAHSSSq(backend,degree+t,b,r+s-t);
            relation := List(relation+value,x -> x mod 2);
            Add(terms,[r+s-t,t]);
        fi;
    od;
    equationArgs := [backend,total-1,relation];
    if Length(arg)=6 then Add(equationArgs,arg[6]); fi;
    result := CallFuncList(koAHSSSolveCochainEquation,equationArgs);
    if result=fail then return fail; fi;
    result.relationCocycle := relation;
    result.inputDegree := degree;
    result.ademPair := [r,s];
    result.rightHandSquares := terms;
    result.scope := "this cochain complex and this input cocycle";
    return result;
end);

# The supplied notes use the reduction of the integral Bockstein for Sq1.
# It need not be cochain-equal to cup_(n-1)(a,a) on a chosen resolution.
InstallGlobalFunction(koAHSSAdem22Primitive, function(arg)
    local backend, degree, a, e, b, relation, equationArgs, result;
    if not Length(arg) in [3,4] then
        Error("usage: koAHSSAdem22Primitive(backend,degree,a[,constraints])");
    fi;
    backend := arg[1]; degree := arg[2];
    a := KOAHSS_AdemCocycle(backend,degree,arg[3]);
    e := List(KOAHSS_EquationDifferential(backend,degree,a),x -> (x/2) mod 2);
    b := koAHSSSq(backend,degree,a,2);
    relation := List(koAHSSSq(backend,degree+2,b,2)
                   +koAHSSSq(backend,degree+1,e,3),x -> x mod 2);
    equationArgs := [backend,degree+3,relation];
    if Length(arg)=4 then Add(equationArgs,arg[4]); fi;
    result := CallFuncList(koAHSSSolveCochainEquation,equationArgs);
    if result=fail then return fail; fi;
    result.relationCocycle := relation;
    result.inputDegree := degree;
    result.bocksteinCocycle := e;
    result.ademPair := [2,2];
    result.scope := "this cochain complex and this input cocycle";
    return result;
end);
