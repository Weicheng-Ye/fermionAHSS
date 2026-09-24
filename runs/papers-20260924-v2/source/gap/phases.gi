# Exact affine solution set for finite-denominator phase cochains.
# Convention: unknown row vector x, equation x * matrix = rhs (mod modulus).
# Extra boundary/coherence equations can be appended as columns.
InstallGlobalFunction(koAHSSSolvePhaseSystem, function(matrix, rhs, modulus)
    local rows, cols, snf, transformed, y, generators, orders, i, d, g,
          reducedModulus, inverse, vec, particular;
    if not IsInt(modulus) or modulus < 1 then
        Error("phase modulus must be a positive integer");
    fi;
    if not IsList(matrix) or not IsList(rhs) or not ForAll(rhs, IsInt) then
        Error("phase matrix and rhs must contain integers");
    fi;
    rows := Length(matrix); cols := Length(rhs);
    if not ForAll(matrix, r -> IsList(r) and Length(r) = cols and ForAll(r, IsInt)) then
        Error("phase matrix rows must have Length(rhs) integer entries");
    fi;
    if rows = 0 then
        if ForAny(rhs, x -> x mod modulus <> 0) then return fail; fi;
        return rec(particular := [], homogeneousGenerators := [],
                   homogeneousOrders := [], modulus := modulus);
    fi;
    if cols = 0 then
        generators := []; orders := [];
        if modulus > 1 then
            for i in [1..rows] do
                vec := List([1..rows], j -> 0); vec[i] := 1;
                Add(generators, vec); Add(orders, modulus);
            od;
        fi;
        return rec(particular := List([1..rows], j -> 0),
                   homogeneousGenerators := generators,
                   homogeneousOrders := orders, modulus := modulus);
    fi;
    snf := SmithNormalFormIntegerMatTransforms(matrix);
    # snf.normal = snf.rowtrans * matrix * snf.coltrans.
    transformed := List(rhs * snf.coltrans, x -> x mod modulus);
    y := List([1..rows], j -> 0);
    generators := []; orders := [];
    for i in [1..snf.rank] do
        d := snf.normal[i][i]; g := Gcd(d, modulus);
        if transformed[i] mod g <> 0 then return fail; fi;
        reducedModulus := QuoInt(modulus, g);
        if reducedModulus > 1 then
            inverse := PowerModInt(QuoInt(d,g), -1, reducedModulus);
            y[i] := (QuoInt(transformed[i],g) * inverse) mod reducedModulus;
        fi;
        if g > 1 then
            vec := List([1..rows], j -> 0); vec[i] := reducedModulus;
            Add(generators, List(vec * snf.rowtrans, x -> x mod modulus));
            Add(orders, g);
        fi;
    od;
    for i in [snf.rank+1..cols] do
        if transformed[i] <> 0 then return fail; fi;
    od;
    if modulus > 1 then
        for i in [snf.rank+1..rows] do
            vec := List([1..rows], j -> 0); vec[i] := 1;
            Add(generators, List(vec * snf.rowtrans, x -> x mod modulus));
            Add(orders, modulus);
        od;
    fi;
    particular := List(y * snf.rowtrans, x -> x mod modulus);
    return rec(particular := particular, homogeneousGenerators := generators,
               homogeneousOrders := orders, modulus := modulus);
end);
