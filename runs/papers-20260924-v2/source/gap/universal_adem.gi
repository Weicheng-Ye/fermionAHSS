# A normalized universal primitive on the simplicial cocycle model K(F_2,n).
# This construction uses the note's interval-cut cups, not a resolution's
# independently chosen higher diagonals.
InstallGlobalFunction(koAHSSUniversalAdem22, function(arg)
    local n, options, limit, sourceDegree, targetDegree, targetBits, required,
          sourceCount, targetCount, degrees, degreeData, mask, coordinates,
          checkCode, fullValues, encodeCoordinates, encode, faceValue,
          restrictCode, degeneracy, isDegenerate, cupCache, cupTerms,
          sourceData, targetData, sourceDegenerates, sourceCodes, sourceIndex,
          code, lowerCode, j, i, k, face, faces, data, indices, count,
          boundaryForms, outerSquare, innerSquare, outerBockstein,
          squareExpressions, bocksteinExpressions, mapFace, squareExpression,
          bocksteinExpression, evalSquare, evalBockstein, obstructionValue,
          obstructionValues, basis, basisRhs, values, equation, rhs,
          unknownValues, normalizedValues, boundaryCodes, coboundaryValue,
          evaluate, pullback, result, rank, normalization, addEquation;
    if not Length(arg) in [1,2] then
        Error("koAHSSUniversalAdem22(n[, options])");
    fi;
    n := arg[1];
    if not IsInt(n) or n < 0 then
        Error("koAHSS: universal Adem degree must be a nonnegative integer");
    fi;
    options := rec();
    if Length(arg) = 2 then options := arg[2]; fi;
    if not IsRecord(options) then Error("koAHSS: Adem options must be a record"); fi;
    limit := 65536;
    if IsBound(options.maxSimplices) then limit := options.maxSimplices; fi;
    if not IsInt(limit) or limit < 1 then
        Error("koAHSS: maxSimplices must be a positive integer");
    fi;
    sourceDegree := n+3;
    targetDegree := n+4;
    normalization:=rec(kind:="unconstrained lexicographic minimum",
        constraintApplied:=false);
    if n=9 then
        normalization:=rec(kind:="chi9_C1",sourceDegree:=12,
            simplexCode:=87112285934295547847080502420714144399368,
            coordinateIndicesZeroBased:=[83,118,162,216],value:=1,
            constraintApplied:=false,
            requirement:="lexicographically least normalized solution with v9(C)=1");
    fi;
    targetBits := Binomial(targetDegree,n);
    if targetBits > LogInt(limit,2) then
        required := fail;
        # Avoid constructing enormous integers merely to explain the guard.
        if targetBits <= 1024 then required := 2^targetBits; fi;
        return rec(status := "resource-limit", inputDegree := n,
            requiredTargetSimplices := required, targetSimplexExponent := targetBits,
            maxSimplices := limit,normalization:=normalization);
    fi;
    sourceCount := 2^Binomial(sourceDegree,n);
    targetCount := 2^targetBits;
    degrees := [];
    mask := vertices -> Sum(vertices, v -> 2^v);
    degreeData := function(r)
        local independent, allFaces, positions, forms, f, form, omitted, z, entry;
        if not IsInt(r) or r < 0 or r > targetDegree then
            Error("koAHSS: simplex degree must lie between zero and n+4");
        fi;
        if IsBound(degrees[r+1]) then return degrees[r+1]; fi;
        independent := List(Combinations([1..r],n), f -> Concatenation([0],f));
        allFaces := Combinations([0..r],n+1);
        positions := [];
        for z in [1..Length(allFaces)] do positions[mask(allFaces[z])+1] := z; od;
        forms := [];
        for f in allFaces do
            if f[1] = 0 then
                Add(forms,[Position(independent,f)]);
            else
                form := [];
                for omitted in [1..Length(f)] do
                    entry := Concatenation([0],f{Filtered([1..Length(f)], z -> z<>omitted)});
                    Add(form,Position(independent,entry));
                od;
                Add(forms,form);
            fi;
        od;
        entry := rec(degree:=r,independent:=independent,faces:=allFaces,
            positions:=positions,forms:=forms,bits:=Length(independent),
            count:=2^Length(independent));
        degrees[r+1] := entry;
        return entry;
    end;
    checkCode := function(r, simplex)
        local d;
        d := degreeData(r);
        if not IsInt(simplex) or simplex < 0 or simplex >= d.count then
            Error("koAHSS: simplex code is outside the cocycle-model range");
        fi;
        return d;
    end;
    coordinates := function(r, simplex)
        local d;
        d := checkCode(r,simplex);
        return List([1..d.bits], t -> QuoInt(simplex,2^(d.bits-t)) mod 2);
    end;
    fullValues := function(r, simplex)
        local d, bits;
        d := checkCode(r,simplex);
        bits := coordinates(r,simplex);
        return List(d.forms, form -> Sum(form,t -> bits[t]) mod 2);
    end;
    encodeCoordinates := function(r, bits)
        local d, simplex, bit;
        d := degreeData(r);
        if not IsList(bits) or Length(bits) <> d.bits
           or not ForAll(bits, x -> x in [0,1]) then
            Error("koAHSS: independent simplex coordinates must be binary");
        fi;
        simplex := 0;
        for bit in bits do simplex := 2*simplex+bit; od;
        return simplex;
    end;
    encode := function(r, allValues)
        local d, simplex;
        d := degreeData(r);
        if not IsList(allValues) or Length(allValues) <> Length(d.faces)
           or not ForAll(allValues, x -> x in [0,1]) then
            Error("koAHSS: simplex face values have incorrect dimension or are not binary");
        fi;
        # Containing-0 faces occur first in the full lexicographic face order.
        simplex := encodeCoordinates(r,allValues{[1..d.bits]});
        if fullValues(r,simplex) <> allValues then
            Error("koAHSS: simplex face values do not define an n-cocycle");
        fi;
        return simplex;
    end;
    faceValue := function(r, simplex, vertices)
        local d, position;
        d := checkCode(r,simplex);
        if not IsList(vertices) or Length(vertices) <> n+1
           or not ForAll(vertices, IsInt) or Set(vertices) <> vertices
           or ForAny(vertices,v -> v<0 or v>r) then
            Error("koAHSS: n-face vertices must be ordered, distinct, and in the simplex");
        fi;
        position := d.positions[mask(vertices)+1];
        return fullValues(r,simplex)[position];
    end;
    restrictCode := function(r, simplex, vertices)
        local d, smaller, allValues, bits;
        d := checkCode(r,simplex);
        if not IsList(vertices) or Length(vertices)=0
           or not ForAll(vertices,IsInt) or Set(vertices)<>vertices
           or ForAny(vertices,v -> v<0 or v>r) then
            Error("koAHSS: restriction vertices must be nonempty, ordered, and in the simplex");
        fi;
        smaller := degreeData(Length(vertices)-1);
        allValues := fullValues(r,simplex);
        bits := List(smaller.independent,
            f -> allValues[d.positions[mask(List(f,v -> vertices[v+1]))+1]]);
        return encodeCoordinates(Length(vertices)-1,bits);
    end;
    degeneracy := function(r, simplex, repeatVertex)
        local d, larger, allValues, bits, f, image;
        d := checkCode(r,simplex);
        if r >= targetDegree or not IsInt(repeatVertex)
           or repeatVertex<0 or repeatVertex>r then
            Error("koAHSS: invalid simplex degeneracy");
        fi;
        larger := degreeData(r+1);
        allValues := fullValues(r,simplex);
        bits := [];
        for f in larger.independent do
            image := List(f,function(v)
                if v<=repeatVertex then return v; else return v-1; fi;
            end);
            if Length(Set(image))<>Length(image) then Add(bits,0);
            else Add(bits,allValues[d.positions[mask(image)+1]]); fi;
        od;
        return encodeCoordinates(r+1,bits);
    end;
    isDegenerate := function(r, simplex)
        local j, smaller;
        checkCode(r,simplex);
        for j in [0..r-1] do
            smaller := restrictCode(r,simplex,Filtered([0..r],v -> v<>j));
            if degeneracy(r-1,smaller,j)=simplex then return true; fi;
        od;
        return false;
    end;

    # Enumerate the finite interval-cut formula once per pair of degrees.
    cupCache := [];
    cupTerms := function(p,q,index)
        local cached, degree, length, cuts, terms, visit;
        if index<0 then return []; fi;
        degree := p+q-index;
        if degree<Maximum(p,q) then return []; fi;
        for cached in cupCache do
            if cached.key=[p,q,index] then return cached.terms; fi;
        od;
        length := index+2;
        cuts := [0];
        terms := [];
        visit := function(position,previous)
            local next, left, right, segment, term, l;
            if position=length then
                cuts[length+1]:=degree;
                left:=[]; right:=[];
                for l in [1..length] do
                    segment:=[cuts[l]..cuts[l+1]];
                    if IsOddInt(l) then UniteSet(left,segment);
                    else UniteSet(right,segment); fi;
                od;
                if Length(left)=p+1 and Length(right)=q+1 then
                    term:=[ShallowCopy(left),ShallowCopy(right)];
                    if term in terms then RemoveSet(terms,term);
                    else AddSet(terms,term); fi;
                fi;
            else
                for next in [previous..degree] do
                    cuts[position+1]:=next;
                    visit(position+1,next);
                od;
            fi;
        end;
        visit(1,0);
        Add(cupCache,rec(key:=[p,q,index],terms:=terms));
        return terms;
    end;

    sourceData := degreeData(sourceDegree);
    targetData := degreeData(targetDegree);
    sourceDegenerates := [];
    for lowerCode in [0..degreeData(sourceDegree-1).count-1] do
        for j in [0..sourceDegree-1] do
            AddSet(sourceDegenerates,degeneracy(sourceDegree-1,lowerCode,j));
        od;
    od;
    sourceCodes := Difference([0..sourceCount-1],sourceDegenerates);
    sourceIndex := List([1..sourceCount],j -> 0);
    for j in [1..Length(sourceCodes)] do sourceIndex[sourceCodes[j]+1]:=j; od;

    # Express all facet restrictions directly in the target's full face values.
    boundaryForms := [];
    for j in [0..targetDegree] do
        face := Filtered([0..targetDegree],v -> v<>j);
        Add(boundaryForms,List(sourceData.independent,
            f -> targetData.positions[mask(List(f,v -> face[v+1]))+1]));
    od;
    boundaryCodes := function(allValues)
        local codes, form, simplex, t;
        codes := [];
        for form in boundaryForms do
            simplex:=0;
            for t in form do simplex:=2*simplex+allValues[t]; od;
            Add(codes,simplex);
        od;
        return codes;
    end;

    innerSquare := cupTerms(n,n,n-2);
    outerSquare := cupTerms(n+2,n+2,n);
    outerBockstein := cupTerms(n+1,n+1,n-2);
    mapFace := function(f,vertices)
        return targetData.positions[mask(List(f,v -> vertices[v+1]))+1];
    end;
    squareExpression := vertices -> List(innerSquare,
        pair -> [mapFace(pair[1],vertices),mapFace(pair[2],vertices)]);
    bocksteinExpression := vertices -> List([1..Length(vertices)],
        omitted -> targetData.positions[mask(vertices{
            Filtered([1..Length(vertices)],t -> t<>omitted)})+1]);
    squareExpressions := List(outerSquare,
        pair -> [squareExpression(pair[1]),squareExpression(pair[2])]);
    bocksteinExpressions := List(outerBockstein,
        pair -> [bocksteinExpression(pair[1]),bocksteinExpression(pair[2])]);
    evalSquare := function(expression,allValues)
        return Sum(expression,pair -> allValues[pair[1]]*allValues[pair[2]]) mod 2;
    end;
    evalBockstein := function(expression,allValues)
        local numerator;
        numerator:=Sum([1..Length(expression)],
            t -> (-1)^(t-1)*allValues[expression[t]]);
        if numerator mod 2<>0 then
            Error("koAHSS: universal cocycle has a nonintegral Bockstein");
        fi;
        return QuoInt(numerator,2) mod 2;
    end;
    obstructionValue := function(allValues)
        local value, pair;
        value:=0;
        for pair in squareExpressions do
            value:=value+evalSquare(pair[1],allValues)*evalSquare(pair[2],allValues);
        od;
        for pair in bocksteinExpressions do
            value:=value+evalBockstein(pair[1],allValues)*evalBockstein(pair[2],allValues);
        od;
        return value mod 2;
    end;

    # Equations have at most n+5 nonzero entries. Highest-index pivots
    # express each pivot variable in earlier variables; zero free variables
    # therefore produce the lexicographically least entire value vector.
    basis:=[]; basisRhs:=[]; obstructionValues:=[]; rank:=0;
    addEquation:=function(equation,rhs)
        local pivot;
        while Length(equation)>0 do
            pivot:=Last(equation);
            if IsBound(basis[pivot]) then
                equation:=Union(Difference(equation,basis[pivot]),
                                Difference(basis[pivot],equation));
                rhs:=(rhs+basisRhs[pivot]) mod 2;
            else
                basis[pivot]:=equation; basisRhs[pivot]:=rhs;
                rank:=rank+1;
                return true;
            fi;
        od;
        return rhs=0;
    end;
    for code in [0..targetCount-1] do
        values:=fullValues(targetDegree,code);
        rhs:=obstructionValue(values);
        obstructionValues[code+1]:=rhs;
        equation:=[];
        for lowerCode in boundaryCodes(values) do
            i:=sourceIndex[lowerCode+1];
            if i<>0 then
                if i in equation then RemoveSet(equation,i); else AddSet(equation,i); fi;
            fi;
        od;
        if not addEquation(equation,rhs) then return fail; fi;
    od;
    if n=9 then
        i:=sourceIndex[normalization.simplexCode+1];
        if i=0 then Error("koAHSS: the note's chi9 calibration simplex is degenerate"); fi;
        if not addEquation([i],1) then return fail; fi;
        normalization.constraintApplied:=true;
    fi;
    unknownValues:=List(sourceCodes,j -> 0);
    for i in [1..Length(sourceCodes)] do
        if IsBound(basis[i]) then
            unknownValues[i]:=(basisRhs[i]+Sum(basis[i],j -> unknownValues[j])) mod 2;
        fi;
    od;
    normalizedValues:=List([1..sourceCount],j -> 0);
    for i in [1..Length(sourceCodes)] do
        normalizedValues[sourceCodes[i]+1]:=unknownValues[i];
    od;
    evaluate := function(simplex)
        checkCode(sourceDegree,simplex);
        return normalizedValues[simplex+1];
    end;
    pullback := function(faceEvaluator)
        local allValues;
        if not IsFunction(faceEvaluator) then
            Error("koAHSS: pullback requires a function on ordered n-face vertex lists");
        fi;
        allValues:=List(sourceData.faces,f -> faceEvaluator(ShallowCopy(f)));
        return evaluate(encode(sourceDegree,allValues));
    end;
    coboundaryValue := function(simplex)
        local allValues;
        allValues:=fullValues(targetDegree,simplex);
        return Sum(boundaryCodes(allValues),evaluate) mod 2;
    end;
    result:=rec(status:="ok",inputDegree:=n,cochainDegree:=sourceDegree,
        equationDegree:=targetDegree,sourceSimplexCount:=sourceCount,
        targetSimplexCount:=targetCount,nondegenerateSourceCodes:=sourceCodes,
        values:=unknownValues,equationRank:=rank,
        solutionDimension:=Length(sourceCodes)-rank,
        normalization:=normalization,
        isUniversalSimplicialHelper:=true,isNormalizedCochain:=true,
        isNormalizedOperation:=false,isStableFamily:=false,
        coordinates:=coordinates,encodeCoordinates:=encodeCoordinates,
        encode:=encode,faceValue:=faceValue,restrict:=restrictCode,
        degeneracy:=degeneracy,isDegenerate:=isDegenerate,
        evaluate:=evaluate,pullback:=pullback,coboundary:=coboundaryValue);
    result.obstruction:=function(simplex)
        checkCode(targetDegree,simplex);
        return obstructionValues[simplex+1];
    end;
    result.faces:=r -> List(degreeData(r).faces,ShallowCopy);
    return result;
end);
