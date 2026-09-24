# Primary cochain operations on an augmented integral HAP resolution.
# Higher diagonals are available both over F_2 and over Z. Their cup products
# do not by themselves specify secondary or tertiary coherence choices.

BindGlobal("KOAHSS_HAP_Reduce", function(terms)
    return List(Filtered(Collected(terms), x -> IsOddInt(x[2])), x -> x[1]);
end);

# Signed tensor terms have their coefficient in position seven.  Collect on
# the six tensor coordinates, retaining all integer multiplicities.
BindGlobal("KOAHSS_HAP_ReduceIntegral", function(terms)
    local sorted, result, term, previous;
    sorted := ShallowCopy(terms);
    Sort(sorted, function(a, b) return a{[1..6]} < b{[1..6]}; end);
    result := [];
    for term in sorted do
        if Length(result) > 0 then
            previous := result[Length(result)];
            if previous{[1..6]} = term{[1..6]} then
                previous[7] := previous[7] + term[7];
            else
                Add(result, ShallowCopy(term));
            fi;
        else
            Add(result, ShallowCopy(term));
        fi;
    od;
    return Filtered(result, t -> t[7] <> 0);
end);

InstallGlobalFunction(koAHSSHAPSpace, function(arg)
    local R, operations, length, dimension, identity, index, contract, property,
          sign, differential, diagonals, diagonal, action, tensorH, cup, data,
          integralDiagonals, integralDiagonal, integralH, integralAction,
          integralCup, naturalBar, naturalBarFactory, native, nativeFactory;
    if Length(arg) < 1 or Length(arg) > 2 then
        Error("koAHSSHAPSpace(resolution[, operations])");
    fi;
    R := arg[1];
    operations := rec();
    if Length(arg) >= 2 then operations := arg[2]; fi;
    if not IsRecord(operations) then Error("operations must be a record"); fi;
    if not IsBoundGlobal("IsHapResolution") then
        Error("load HAP before calling koAHSSHAPSpace");
    fi;
    if not CallFuncList(ValueGlobal("IsHapResolution"), [R]) or not IsBound(R!.homotopy) then
        Error("expected a HAP resolution with a contracting homotopy");
    fi;
    property := ValueGlobal("EvaluateProperty");
    if property(R, "characteristic") <> 0 then
        Error("an integral resolution is required, not a mod-two resolution");
    fi;
    if R!.dimension(0) < 1 then
        Error("HAP adapter requires at least one augmented degree-zero generator");
    fi;
    length := property(R, "length");
    if not IsInt(length) or length < 1 then
        Error("resolution must have a known positive length");
    fi;
    identity := One(R!.group);
    dimension := function(n)
        if n < 0 then return 0; fi;
        if n > length then
            Error("resolution too short: need chain degree ", n,
                  "; available through ", length);
        fi;
        return R!.dimension(n);
    end;
    index := function(g)
        local k;
        k := Position(R!.elts, g);
        if k = fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(g);
            else Add(R!.elts, g); fi;
            k := Position(R!.elts, g);
            if k = fail then Error("could not index resolution group element"); fi;
        fi;
        return k;
    end;
    # Every degree-zero generator has augmentation 1, and the homotopy
    # contracts to the first generator at the identity. Multiple vertex
    # orbits, as in SGC resolutions, are allowed under this convention.
    # Always use positive basis indices: not every HAP homotopy handles signs.
    contract := function(n, j, g)
        return R!.homotopy(n, [j, index(g)]);
    end;
    sign := function(s, g)
        local word;
        if IsInt(s) and s = 0 then return 1; fi;
        word := contract(0, 1, g);
        return (-1) ^ (Sum(word, x -> SignInt(x[1]) * s[AbsInt(x[1])]) mod 2);
    end;
    differential := function(n, s)
        local mat, i, term;
        mat := List([1..dimension(n)], j -> List([1..dimension(n+1)], k -> 0));
        for i in [1..dimension(n+1)] do
            for term in R!.boundary(n+1, i) do
                mat[AbsInt(term[1])][i] := mat[AbsInt(term[1])][i]
                    + SignInt(term[1]) * sign(s, R!.elts[term[2]]);
            od;
        od;
        return mat;
    end;
    # Tensor terms [p,j,g,q,k,h] represent g e_j^p tensor h e_k^q.
    # Cropping either factor beyond the resolution length cannot affect any
    # retained component: tensor H only raises degrees, and swaps preserve max.
    tensorH := function(terms)
        local result, term, u;
        result := [];
        for term in terms do
            if term[1] < length then
                for u in contract(term[1], term[2], term[3]) do
                    Add(result, [term[1]+1, AbsInt(u[1]), R!.elts[u[2]],
                                 term[4], term[5], term[6]]);
                od;
            fi;
            if term[1] = 0 and term[4] < length then
                for u in contract(term[4], term[5], term[6]) do
                    Add(result, [0, 1, identity,
                                 term[4]+1, AbsInt(u[1]), R!.elts[u[2]]]);
                od;
            fi;
        od;
        return KOAHSS_HAP_Reduce(result);
    end;
    action := function(terms, g)
        return List(terms, t -> [t[1],t[2],g*t[3],t[4],t[5],g*t[6]]);
    end;
    diagonals := [];
    diagonal := function(i, n, j)
        local rhs, t, prev;
        if i < 0 or n < 0 then return []; fi;
        if not IsBound(diagonals[i+1]) then diagonals[i+1] := []; fi;
        if not IsBound(diagonals[i+1][n+1]) then diagonals[i+1][n+1] := []; fi;
        if IsBound(diagonals[i+1][n+1][j]) then
            return diagonals[i+1][n+1][j];
        fi;
        if i = 0 and n = 0 then
            rhs := [[0,j,identity,0,j,identity]];
        else
            rhs := [];
            if n > 0 then
                for t in R!.boundary(n,j) do
                    Append(rhs, action(diagonal(i,n-1,AbsInt(t[1])), R!.elts[t[2]]));
                od;
            fi;
            if i > 0 then
                prev := diagonal(i-1,n,j);
                Append(rhs, prev);
                Append(rhs, List(prev, t -> [t[4],t[5],t[6],t[1],t[2],t[3]]));
            fi;
            rhs := tensorH(KOAHSS_HAP_Reduce(rhs));
        fi;
        diagonals[i+1][n+1][j] := rhs;
        return rhs;
    end;
    cup := function(i, p, a, q, b)
        local n, result, j, t;
        n := p+q-i;
        if n < 0 then return []; fi;
        result := List([1..dimension(n)], j -> 0);
        if i < 0 then return result; fi;
        if Length(a) <> dimension(p) or Length(b) <> dimension(q) then
            Error("cup input has wrong cochain dimension");
        fi;
        for j in [1..Length(result)] do
            for t in diagonal(i,n,j) do
                if t[1] = p and t[4] = q then
                    result[j] := (result[j]+a[t[2]]*b[t[5]]) mod 2;
                fi;
            od;
        od;
        return result;
    end;

    # The tensor contraction H=h tensor 1 + eta epsilon tensor h respects
    # integral boundary/homotopy signs.  Its second summand is present only
    # when the first factor has degree zero, so no extra Koszul sign occurs.
    integralH := function(terms)
        local result, term, u;
        result := [];
        for term in terms do
            if term[1] < length then
                for u in contract(term[1], term[2], term[3]) do
                    Add(result, [term[1]+1, AbsInt(u[1]), R!.elts[u[2]],
                        term[4], term[5], term[6], term[7] * SignInt(u[1])]);
                od;
            fi;
            if term[1] = 0 and term[4] < length then
                for u in contract(term[4], term[5], term[6]) do
                    Add(result, [0, 1, identity, term[4]+1, AbsInt(u[1]),
                        R!.elts[u[2]], term[7] * SignInt(u[1])]);
                od;
            fi;
        od;
        return KOAHSS_HAP_ReduceIntegral(result);
    end;
    integralAction := function(terms, g, coefficient)
        return List(terms, t -> [t[1], t[2], g*t[3], t[4], t[5], g*t[6],
            coefficient*t[7]]);
    end;
    integralDiagonals := [];
    integralDiagonal := function(i, n, j)
        local rhs, t, previous;
        if i < 0 or n < 0 then return []; fi;
        if not IsBound(integralDiagonals[i+1]) then integralDiagonals[i+1] := []; fi;
        if not IsBound(integralDiagonals[i+1][n+1]) then
            integralDiagonals[i+1][n+1] := [];
        fi;
        if IsBound(integralDiagonals[i+1][n+1][j]) then
            return integralDiagonals[i+1][n+1][j];
        fi;
        if i = 0 and n = 0 then
            rhs := [[0, j, identity, 0, j, identity, 1]];
        else
            rhs := [];
            if n > 0 then
                for t in R!.boundary(n, j) do
                    Append(rhs, integralAction(integralDiagonal(i, n-1, AbsInt(t[1])),
                        R!.elts[t[2]], SignInt(t[1])));
                od;
            fi;
            if i > 0 then
                previous := integralDiagonal(i-1, n, j);
                Append(rhs, List(previous, t -> [t[1], t[2], t[3], t[4], t[5], t[6],
                    (-1)^n * t[7]]));
                Append(rhs, List(previous, t -> [t[4], t[5], t[6], t[1], t[2], t[3],
                    (-1)^(n+i+t[1]*t[4]) * t[7]]));
            fi;
            rhs := integralH(KOAHSS_HAP_ReduceIntegral(rhs));
        fi;
        integralDiagonals[i+1][n+1][j] := rhs;
        return rhs;
    end;
    integralCup := function(i, p, a, q, b, sFirst, sSecond)
        local n, result, j, t;
        n := p+q-i;
        if n < 0 then return []; fi;
        result := List([1..dimension(n)], j -> 0);
        if i < 0 then return result; fi;
        KOAHSS_CC_CheckVector(a, dimension(p), "integral cup first input");
        KOAHSS_CC_CheckVector(b, dimension(q), "integral cup second input");
        for j in [1..Length(result)] do
            for t in integralDiagonal(i, n, j) do
                if t[1] = p and t[4] = q then
                    result[j] := result[j] + t[7] * sign(sFirst, t[3]) * a[t[2]]
                        * sign(sSecond, t[6]) * b[t[5]];
                fi;
            od;
        od;
        return result;
    end;
    # The simplicial comparison is constructed only if a natural cochain
    # operation requests it; ordinary primary computations keep their costs.
    naturalBar := fail;
    naturalBarFactory := function()
        if naturalBar=fail then
            if not IsBoundGlobal("KOAHSS_NaturalBarTransport") then
                Error("koAHSS: the natural bar transport module is unavailable");
            fi;
            naturalBar := ValueGlobal("KOAHSS_NaturalBarTransport")(R);
        fi;
        return naturalBar;
    end;
    # Higher tensor cells are built directly from R's contraction, on demand.
    # This diagnostic engine is separate from the calibrated simplicial model.
    native := fail;
    nativeFactory := function()
        if native=fail then native:=koAHSSNativeCoherence(R); fi;
        return native;
    end;
    data := rec(dimension := dimension, differential := differential,
                cupMod2 := cup, cupIntegral := integralCup, operations := operations,
                untwistedConstantsLift := true,naturalBar:=naturalBarFactory,
                nativeCoherence:=nativeFactory,naturalTransport:=naturalBarFactory);
    return koAHSSCochainSpace(data);
end);
