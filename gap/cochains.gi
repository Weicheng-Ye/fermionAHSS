# Finite-rank cochains. Matrices act on row vectors on the right.
# Higher operations require a supplied, coherently normalized implementation.

BindGlobal("KOAHSS_CC_Multiply", function(vector, matrix, width)
    if Length(vector) = 0 then return List([1..width], i -> 0); fi;
    if width = 0 then return []; fi;
    return vector * matrix;
end);

BindGlobal("KOAHSS_CC_CheckVector", function(vector, dimension, label)
    if not IsList(vector) or Length(vector) <> dimension
       or not ForAll(vector, IsInt) then
        Error(Concatenation("koAHSS: ", label,
              " must be an integer cochain of the required dimension"));
    fi;
end);

# Present a quotient of a free cocycle lattice by integral relation rows.
# Each nontrivial Smith coordinate becomes one Pcp generator.
BindGlobal("KOAHSS_CC_Quotient", function(cycles, relations, dimension, modTwo)
    local rank, coordinates, row, solution, smith, transform, inverse,
          invariants, active, i, diagonal, group, generators, result;
    rank := Length(cycles);
    coordinates := [];
    for row in relations do
        if rank = 0 then
            if not ForAll(row, x -> x = 0) then
                Error("koAHSS: a boundary is not in the cocycle lattice");
            fi;
        else
            solution := SolutionMat(cycles, row);
            if solution = fail or not ForAll(solution, IsInt) then
                Error("koAHSS: a boundary is not in the integral cocycle lattice");
            fi;
            Add(coordinates, solution);
        fi;
    od;
    if rank = 0 then
        transform := [];
        inverse := [];
        diagonal := [];
    elif Length(coordinates) = 0 then
        transform := IdentityMat(rank);
        inverse := transform;
        diagonal := List([1..rank], i -> 0);
    else
        smith := SmithNormalFormIntegerMatTransforms(coordinates);
        transform := smith.coltrans;
        inverse := Inverse(transform);
        diagonal := List([1..rank], i -> 0);
        for i in [1..Minimum(rank, Length(smith.normal))] do
            diagonal[i] := AbsInt(smith.normal[i][i]);
        od;
    fi;
    active := Filtered([1..rank], i -> diagonal[i] <> 1);
    invariants := List(active, i -> diagonal[i]);
    group := AbelianGroup(IsPcpGroup, invariants);
    generators := GeneratorsOfGroup(group);
    result := rec(group := group, invariants := invariants,
                  dimension := dimension, modTwo := modTwo);
    result.represent := function(element)
        local values, exponents, j, cochain;
        if not element in group then Error("koAHSS: wrong cohomology group"); fi;
        values := List([1..rank], i -> 0);
        exponents := Exponents(element);
        for j in [1..Length(active)] do values[active[j]] := exponents[j]; od;
        values := KOAHSS_CC_Multiply(values, inverse, rank);
        cochain := KOAHSS_CC_Multiply(values, cycles, dimension);
        if modTwo then cochain := List(cochain, x -> x mod 2); fi;
        return cochain;
    end;
    result.class := function(cochain)
        local values, element, j;
        KOAHSS_CC_CheckVector(cochain, dimension, "cocycle");
        if rank = 0 then
            if not ForAll(cochain, x -> x = 0) then
                Error("koAHSS: cochain is not a cocycle");
            fi;
            return One(group);
        fi;
        values := SolutionMat(cycles, cochain);
        if values = fail or not ForAll(values, IsInt) then
            Error("koAHSS: cochain is not a cocycle in the indicated coefficients");
        fi;
        values := values * transform;
        element := One(group);
        for j in [1..Length(active)] do
            element := element * generators[j]^values[active[j]];
        od;
        return element;
    end;
    return result;
end);

InstallGlobalFunction(koAHSSCochainSpace, function(data)
    local space;
    if not IsRecord(data) or not IsBound(data.dimension)
       or not IsFunction(data.dimension) or not IsBound(data.differential)
       or not IsFunction(data.differential) then
        Error("koAHSS: cochain data require dimension(n) and differential(n,s)");
    fi;
    if IsBound(data.cupMod2) and not IsFunction(data.cupMod2) then
        Error("koAHSS: cupMod2 must be a function");
    fi;
    if IsBound(data.operations) and not IsRecord(data.operations) then
        Error("koAHSS: operations must be a record of functions");
    fi;
    if IsBound(data.cupIntegral) and not IsFunction(data.cupIntegral) then
        Error("koAHSS: cupIntegral must be a function");
    fi;
    if IsBound(data.nativeCoherence) and not IsFunction(data.nativeCoherence) then
        Error("koAHSS: nativeCoherence must be a function");
    fi;
    if IsBound(data.naturalTransport) and not IsFunction(data.naturalTransport) then
        Error("koAHSS: naturalTransport must be a function");
    fi;
    space := rec(cochainData := data);
    space.koAHSS := function(s, omega, maxDegree)
        local backend, dimensions, ordinary, twisted, cache, zeroS,
              dimension, getMatrix, validateTwist, cohomologyData,
              cup, cupIntegral, reduce, differentialOnCochain, operation, operations;
        dimensions := [];
        ordinary := [];
        twisted := [];
        cache := [];
        dimension := function(n)
            local value;
            if n < 0 then return 0; fi;
            if not IsBound(dimensions[n+1]) then
                value := data.dimension(n);
                if not IsInt(value) or value < 0 then
                    Error("koAHSS: dimension(n) must be a nonnegative integer");
                fi;
                dimensions[n+1] := value;
            fi;
            return dimensions[n+1];
        end;
        if s = 0 then s := List([1..dimension(1)],i -> 0); fi;
        if omega = 0 then omega := List([1..dimension(2)],i -> 0); fi;
        KOAHSS_CC_CheckVector(s, dimension(1), "s");
        KOAHSS_CC_CheckVector(omega, dimension(2), "omega");
        if not ForAll(s, x -> x in [0,1])
           or not ForAll(omega, x -> x in [0,1]) then
            Error("koAHSS: twists must use binary cocycle representatives");
        fi;
        s := ShallowCopy(s);
        omega := ShallowCopy(omega);
        zeroS := List(s, x -> 0);
        getMatrix := function(n, sign)
            local matrices, matrix, d, e, previous, i;
            if n < 0 then return []; fi;
            if sign then matrices := twisted; else matrices := ordinary; fi;
            if not IsBound(matrices[n+1]) then
                d := dimension(n);
                e := dimension(n+1);
                if sign then matrix := data.differential(n, s);
                else matrix := data.differential(n, zeroS); fi;
                if not IsList(matrix) or Length(matrix) <> d
                   or not ForAll(matrix, row -> IsList(row)
                       and Length(row) = e and ForAll(row, IsInt)) then
                    Error("koAHSS: differential has incorrect dimensions or noninteger entries");
                fi;
                matrix := List(matrix, ShallowCopy);
                matrices[n+1] := matrix;
                if n > 0 then
                    previous := getMatrix(n-1, sign);
                    for i in [1..Length(previous)] do
                        if not ForAll(KOAHSS_CC_Multiply(previous[i], matrix, e),
                                      x -> x = 0) then
                            Error("koAHSS: the supplied coboundary does not square to zero");
                        fi;
                    od;
                fi;
                if sign then
                    previous := getMatrix(n, false);
                    for i in [1..d] do
                        if not ForAll(matrix[i]-previous[i], x -> x mod 2 = 0) then
                            Error("koAHSS: sign and ordinary coboundaries disagree modulo two");
                        fi;
                    od;
                fi;
            fi;
            return matrices[n+1];
        end;
        validateTwist := function(cochain, n)
            if not ForAll(KOAHSS_CC_Multiply(cochain, getMatrix(n,false),
                           dimension(n+1)), x -> x mod 2 = 0) then
                Error("koAHSS: s and omega must be mod-two cocycles");
            fi;
        end;
        validateTwist(s,1);
        validateTwist(omega,2);
        reduce := vector -> List(vector, x -> x mod 2);
        cohomologyData := function(n, q)
            local key, d, e, matrix, cycles, relations, augmented,
                  nullspace, info, i, modTwo;
            if not q in [0,-1,-2,-3,-4] then
                Error("koAHSS: this cochain backend supports rows -4 through 0");
            fi;
            if n < 0 or q = -3 then
                if not IsBound(backend.zero) then
                    backend.zero := KOAHSS_CC_Quotient([],[],0,false);
                fi;
                return backend.zero;
            fi;
            modTwo := q in [-1,-2];
            if modTwo then key := 2*n+2; else key := 2*n+1; fi;
            if IsBound(cache[key]) then return cache[key]; fi;
            d := dimension(n);
            e := dimension(n+1);
            matrix := getMatrix(n, not modTwo);
            if n = 0 then relations := [];
            else relations := List(getMatrix(n-1,not modTwo),ShallowCopy); fi;
            if d = 0 then cycles := [];
            elif e = 0 then cycles := IdentityMat(d);
            elif modTwo then
                augmented := Concatenation(matrix, -2*IdentityMat(e));
                nullspace := NullspaceIntMat(augmented);
                cycles := List(nullspace, row -> row{[1..d]});
            else cycles := NullspaceIntMat(matrix); fi;
            if modTwo and d > 0 then
                Append(relations, 2*IdentityMat(d));
            fi;
            info := KOAHSS_CC_Quotient(cycles,relations,d,modTwo);
            info.degree := n;
            info.row := q;
            cache[key] := info;
            return info;
        end;
        cup := function(index, degreeA, a, degreeB, b)
            local outDegree, value;
            outDegree := degreeA + degreeB - index;
            if outDegree < 0 then return []; fi;
            if index < 0 or ForAll(a,x -> x = 0) or ForAll(b,x -> x = 0) then
                return List([1..dimension(outDegree)], i -> 0);
            fi;
            if not IsBound(data.cupMod2) then
                Error("koAHSS: primary operations require data.cupMod2(i,degreeA,a,degreeB,b)");
            fi;
            value := data.cupMod2(index,degreeA,a,degreeB,b);
            KOAHSS_CC_CheckVector(value,dimension(outDegree),"higher cup output");
            return reduce(value);
        end;
        differentialOnCochain := function(n, a, sign)
            return KOAHSS_CC_Multiply(a,getMatrix(n,sign),dimension(n+1));
        end;
        cupIntegral := function(index, degreeA, a, degreeB, b, firstSign, secondSign)
            local outDegree, value;
            outDegree := degreeA + degreeB - index;
            if outDegree < 0 then return []; fi;
            KOAHSS_CC_CheckVector(a,dimension(degreeA),"integral cup input");
            KOAHSS_CC_CheckVector(b,dimension(degreeB),"integral cup input");
            if index < 0 or ForAll(a,x -> x=0) or ForAll(b,x -> x=0) then
                return List([1..dimension(outDegree)], i -> 0);
            fi;
            if not IsBound(data.cupIntegral) then return fail; fi;
            if firstSign = true then firstSign := s;
            elif firstSign = false then firstSign := zeroS; fi;
            if secondSign = true then secondSign := s;
            elif secondSign = false then secondSign := zeroS; fi;
            value := data.cupIntegral(index,degreeA,a,degreeB,b,firstSign,secondSign);
            if value = fail then return fail; fi;
            KOAHSS_CC_CheckVector(value,dimension(outDegree),"integral cup output");
            return value;
        end;
        operations := rec();
        if IsBound(data.operations) then operations := data.operations; fi;
        backend := rec();
        backend.dimension := dimension;
        backend.cohomologyData := cohomologyData;
        backend.data := cohomologyData;
        backend.cohomology := function(n,q) return cohomologyData(n,q).group; end;
        backend.coboundary := differentialOnCochain;
        backend.cupMod2 := cup;
        backend.cupIntegral := cupIntegral;
        backend.hasIntegralCups := IsBound(data.cupIntegral);
        backend.twists := rec(s := s, omega := omega);
        backend.tertiaryCoefficient := 1;
        if IsBound(operations.usesNaturalT) and operations.usesNaturalT=true then
            if not IsBound(operations.T) or not IsFunction(operations.T)
               or not IsBound(operations.tertiaryReference) then
                Error("koAHSS: natural tertiary metadata requires its final T callback");
            fi;
            backend.tertiaryCoefficient:=0;
            backend.tertiaryReference:=operations.tertiaryReference;
        fi;
        if IsBound(data.naturalTransport) then backend.naturalTransport:=data.naturalTransport; fi;
        if IsBound(data.naturalBar) then backend.naturalBar:=data.naturalBar; fi;
        if IsBound(data.nativeCoherence) then
            backend.nativeCoherence:=data.nativeCoherence;
        fi;
        backend.usesNaturalPrimary:=IsBound(operations.useNaturalPrimary)
            and operations.useNaturalPrimary=true;
        backend.primary := function(name, n, input)
            local a, sq1, sq2, value;
            KOAHSS_CC_CheckVector(input,dimension(n),"primary input");
            if name = "Dbar" and not ForAll(
                differentialOnCochain(n,input,true),x -> x = 0) then
                Error("koAHSS: Dbar input must be a sign-integral cocycle");
            fi;
            a := reduce(input);
            value := differentialOnCochain(n,a,false);
            if not ForAll(value,x -> x mod 2 = 0) then
                Error("koAHSS: primary input is not a mod-two cocycle");
            fi;
            if backend.usesNaturalPrimary then
                return CallFuncList(ValueGlobal("KOAHSS_NaturalPrimary"),[backend,name,n,input]);
            fi;
            sq1 := reduce(List(value,x -> x/2));
            sq2 := cup(n-2,n,a,n,a);
            value := sq2 + cup(0,2,omega,n,a);
            if name = "D" or name = "Dbar" then
                return reduce(value+cup(0,1,s,n+1,sq1));
            elif name = "Dtilde" then
                value := differentialOnCochain(n+2,reduce(value),true);
                if not ForAll(value,x -> x mod 2 = 0) then
                    Error("koAHSS: Dtilde Bockstein numerator is not even");
                fi;
                return List(value,x -> x/2);
            fi;
            Error("koAHSS: unknown primary operation");
        end;
        operation := function(r,n,q,cochain,source,target)
            local name, value, context, correction;
            if r = 2 and q = 0 then name := "Dbar";
            elif r = 2 and q = -1 then name := "D";
            elif r = 3 and q = -2 then name := "Dtilde";
            elif r = 3 and q = 0 then name := "Tau";
            elif r = 4 and q = -1 then name := "Psi";
            elif r = 5 and q = 0 then name := "T";
            else return List([1..dimension(n+r)],i -> 0); fi;
            if name in ["Dbar","D","Dtilde"] then
                return backend.primary(name,n,cochain);
            fi;
            # The untwisted coefficient classes in column zero pull back from
            # a point and already lift to ko. HAP group cochains have this unit.
            if n = 0 and IsBound(data.untwistedConstantsLift)
               and data.untwistedConstantsLift = true
               and ForAll(s,x->x=0) and ForAll(omega,x->x=0)
               and not IsBound(operations.(name))
               and (name <> "T" or not IsBound(operations.TReference)) then
                return List([1..dimension(n+r)],i->0);
            fi;
            context := rec(degree:=n, cochain:=cochain, s:=s, omega:=omega,
                backend:=backend, source:=source, target:=target,
                tertiaryCoefficient:=backend.tertiaryCoefficient);
            if IsBound(backend.tertiaryReference) then
                context.tertiaryReference:=backend.tertiaryReference;
            fi;
            if name = "T" and not IsBound(operations.T) then
                if not IsBound(operations.TReference) then
                    return rec(status:="unresolved",operations:=["T"],
                        reasons:=["T reference operation is unavailable; its binary correction coefficient is fixed to 1"]);
                fi;
                if not IsFunction(operations.TReference) then
                    Error("koAHSS: TReference must be a function");
                fi;
                value := operations.TReference(context);
                if value = fail or KOAHSS_IsUnresolved(value) then return value; fi;
                KOAHSS_CC_CheckVector(value,dimension(n+5),"T reference output");
                correction := koAHSSTertiaryCorrection(backend,n,cochain);
                if correction = fail then
                    return rec(status:="unresolved",operations:=["T"],
                        reasons:=["T with coefficient 1 needs signed integral cup products"]);
                fi;
                value := value + correction;
            else
                if not IsBound(operations.(name)) then return fail; fi;
                if not IsFunction(operations.(name)) then
                    Error("koAHSS: ", name, " must be a function");
                fi;
                # An explicit T callback supplies the final operation. Do not
                # add the convention's correction twice to that legacy API.
                value := operations.(name)(context);
                if value = fail or KOAHSS_IsUnresolved(value) then return value; fi;
            fi;
            KOAHSS_CC_CheckVector(value,dimension(n+r),"higher operation output");
            return value;
        end;
        backend.differential := function(r,n,q,source,target)
            local generators, images, sourceData, targetData, element,
                  baseElement, cochain, value, hom;
            generators := GeneratorsOfGroup(source.group);
            images := [];
            if IsTrivial(source.group) or IsTrivial(target.group) then
                return GroupHomomorphismByImages(source.group,target.group,
                    generators,List(generators,g -> One(target.group)));
            fi;
            sourceData := cohomologyData(n,q);
            targetData := cohomologyData(n+r,q-r+1);
            for element in generators do
                baseElement := source.lift(element);
                cochain := sourceData.represent(baseElement);
                value := operation(r,n,q,cochain,source,target);
                if value = fail or KOAHSS_IsUnresolved(value) then return value; fi;
                baseElement := targetData.class(value);
                Add(images,target.project(baseElement));
            od;
            hom := GroupHomomorphismByImages(source.group,target.group,generators,images);
            if hom = fail then Error("koAHSS: operation images do not define a group homomorphism"); fi;
            return hom;
        end;
        return backend;
    end;
    return space;
end);
