# Sparse comparison with the normalized homogeneous bar resolution.
# Chain maps are integral and equivariant. Cochain-facing methods use LOCAL
# values in the fiber at the first vertex, so binary lifts stay binary.
BindGlobal("KOAHSS_NaturalBarTransport", function(R)
    local identity, length, dimension, index, checkSimplex, anchor, degenerate,
          reduceR, reduceBar, actR, actBar, barCone, fCache, gCache, hCache,
          key, f, g, homotopy, character, characterCache, evaluateR, pullback,
          result;
    if not IsBoundGlobal("IsHapResolution")
       or not CallFuncList(ValueGlobal("IsHapResolution"),[R])
       or not IsBound(R!.homotopy) then
        Error("koAHSS: natural bar transport requires a HAP resolution with homotopy");
    fi;
    if ValueGlobal("EvaluateProperty")(R,"characteristic") <> 0 then
        Error("koAHSS: natural bar transport requires an integral resolution");
    fi;
    identity := One(R!.group);
    length := ValueGlobal("EvaluateProperty")(R,"length");
    if not IsInt(length) or length < 1 or R!.dimension(0) < 1 then
        Error("koAHSS: natural bar transport requires a known positive resolution length");
    fi;
    dimension := function(n)
        if not IsInt(n) or n < 0 or n > length then
            Error("koAHSS: natural bar transport requested an unavailable resolution degree");
        fi;
        return R!.dimension(n);
    end;
    index := function(element)
        local position;
        position := Position(R!.elts,element);
        if position = fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(element);
            else Add(R!.elts,element); fi;
            position := Position(R!.elts,element);
            if position = fail then Error("koAHSS: could not index a bar vertex"); fi;
        fi;
        return position;
    end;
    checkSimplex := function(simplex)
        if not IsList(simplex) or IsEmpty(simplex)
           or not ForAll(simplex,element -> IsMultiplicativeElementWithInverse(element)
               and FamilyObj(element)=FamilyObj(identity)) then
            Error("koAHSS: bar simplex vertices must have the resolution group element type");
        fi;
        # Vertices come from R.elts and their products. Do not ask GAP to test
        # membership in R.group: some infinite affine crystallographic groups
        # implement that query via an unavailable group Enumerator.
    end;
    anchor := function(simplex)
        local inverse;
        if simplex[1]=identity then return simplex; fi;
        inverse:=simplex[1]^-1;
        return List(simplex,element -> inverse*element);
    end;
    degenerate := simplex -> ForAny([2..Length(simplex)],
        j -> simplex[j]=simplex[j-1]);
    # R terms are [positive basis index, group element, integer coefficient].
    reduceR := function(terms)
        local sorted, answer, term, previous;
        sorted := ShallowCopy(terms);
        Sort(sorted,function(a,b) return a{[1,2]} < b{[1,2]}; end);
        answer := [];
        for term in sorted do
            if not IsEmpty(answer) then
                previous := Last(answer);
                if previous{[1,2]}=term{[1,2]} then
                    previous[3] := previous[3]+term[3];
                else Add(answer,ShallowCopy(term)); fi;
            else Add(answer,ShallowCopy(term)); fi;
        od;
        return Filtered(answer,term -> term[3]<>0);
    end;
    # Bar terms are [integer coefficient, homogeneous vertex list].
    reduceBar := function(terms)
        local sorted, answer, term, previous;
        sorted := Filtered(terms,term -> term[1]<>0 and not degenerate(term[2]));
        Sort(sorted,function(a,b) return a[2] < b[2]; end);
        answer := [];
        for term in sorted do
            if not IsEmpty(answer) then
                previous := Last(answer);
                if previous[2]=term[2] then previous[1] := previous[1]+term[1];
                else Add(answer,[term[1],ShallowCopy(term[2])]); fi;
            else Add(answer,[term[1],ShallowCopy(term[2])]); fi;
        od;
        return Filtered(answer,term -> term[1]<>0);
    end;
    actR := function(terms,element,coefficient)
        return List(terms,t -> [t[1],element*t[2],coefficient*t[3]]);
    end;
    actBar := function(terms,element,coefficient)
        return List(terms,t -> [coefficient*t[1],List(t[2],v -> element*v)]);
    end;
    barCone := terms -> reduceBar(List(terms,
        t -> [t[1],Concatenation([identity],t[2])]));
    fCache := rec(); gCache := []; hCache := rec(); characterCache := rec();
    key := simplex -> JoinStringsWithSeparator(List(simplex,
        element -> String(index(element))),"_");
    f := function(simplex)
        local normalized, cacheKey, n, rhs, j, face, term, contraction, answer;
        checkSimplex(simplex); n := Length(simplex)-1; dimension(n);
        if degenerate(simplex) then return []; fi;
        normalized := anchor(simplex); cacheKey := key(normalized);
        if not IsBound(fCache.(cacheKey)) then
            if n=0 then answer := [[1,identity,1]];
            else
                rhs := [];
                for j in [1..n+1] do
                    face := normalized{Filtered([1..n+1],k -> k<>j)};
                    Append(rhs,actR(f(face),identity,(-1)^(j-1)));
                od;
                answer := [];
                for term in reduceR(rhs) do
                    contraction := R!.homotopy(n-1,[term[1],index(term[2])]);
                    for j in contraction do
                        Add(answer,[AbsInt(j[1]),R!.elts[j[2]],term[3]*SignInt(j[1])]);
                    od;
                od;
                answer := reduceR(answer);
            fi;
            MakeImmutable(answer); fCache.(cacheKey) := answer;
        fi;
        if simplex[1]=identity then return fCache.(cacheKey); fi;
        return actR(fCache.(cacheKey),simplex[1],1);
    end;
    g := function(n,j)
        local rhs, term, answer;
        if not IsInt(j) or j < 1 or j > dimension(n) then
            Error("koAHSS: invalid resolution basis index for bar comparison");
        fi;
        if not IsBound(gCache[n+1]) then gCache[n+1] := []; fi;
        if not IsBound(gCache[n+1][j]) then
            if n=0 then answer := [[1,[identity]]];
            else
                rhs := [];
                for term in R!.boundary(n,j) do
                    Append(rhs,actBar(g(n-1,AbsInt(term[1])),
                        R!.elts[term[2]],SignInt(term[1])));
                od;
                answer := barCone(reduceBar(rhs));
            fi;
            MakeImmutable(answer); gCache[n+1][j] := answer;
        fi;
        return gCache[n+1][j];
    end;
    homotopy := function(simplex)
        local normalized, cacheKey, n, rhs, term, j, face, answer;
        checkSimplex(simplex); n := Length(simplex)-1; dimension(n);
        if degenerate(simplex) then return []; fi;
        normalized := anchor(simplex); cacheKey := key(normalized);
        if not IsBound(hCache.(cacheKey)) then
            rhs := [[1,normalized]];
            for term in f(normalized) do
                Append(rhs,actBar(g(n,term[1]),term[2],-term[3]));
            od;
            if n>0 then
                for j in [1..n+1] do
                    face := normalized{Filtered([1..n+1],k -> k<>j)};
                    Append(rhs,actBar(homotopy(face),identity,(-1)^j));
                od;
            fi;
            answer := barCone(reduceBar(rhs));
            MakeImmutable(answer); hCache.(cacheKey) := answer;
        fi;
        if simplex[1]=identity then return hCache.(cacheKey); fi;
        return actBar(hCache.(cacheKey),simplex[1],1);
    end;
    character := function(sign,element)
        local cacheKey, contraction;
        if IsInt(sign) and sign=0 then return 1; fi;
        if not IsList(sign) or Length(sign)<>dimension(1)
           or not ForAll(sign,x -> x in [0,1]) then
            Error("koAHSS: a bar sign must be zero or a binary degree-one cochain");
        fi;
        cacheKey := Concatenation(JoinStringsWithSeparator(List(sign,String),""),
            "_",String(index(element)));
        if not IsBound(characterCache.(cacheKey)) then
            contraction := R!.homotopy(0,[1,index(element)]);
            characterCache.(cacheKey) := (-1)^(
                Sum(contraction,t -> SignInt(t[1])*sign[AbsInt(t[1])]) mod 2);
        fi;
        return characterCache.(cacheKey);
    end;
    evaluateR := function(n,vector,sign,simplex)
        if not IsList(vector) or Length(vector)<>dimension(n)
           or not ForAll(vector,IsInt) then
            Error("koAHSS: invalid resolution cochain for bar evaluation");
        fi;
        checkSimplex(simplex);
        if Length(simplex)<>n+1 then Error("koAHSS: wrong simplex degree for bar evaluation"); fi;
        return Sum(f(anchor(simplex)),term -> term[3]*character(sign,term[2])*vector[term[1]]);
    end;
    pullback := function(n,cochain,sign)
        local vector, j, term, value;
        if not IsFunction(cochain) then Error("koAHSS: bar cochain must be a function"); fi;
        vector := List([1..dimension(n)],j -> 0);
        for j in [1..Length(vector)] do
            for term in g(n,j) do
                value := cochain(term[2]);
                if not IsInt(value) then Error("koAHSS: bar cochains must return integers"); fi;
                vector[j] := vector[j]+term[1]*character(sign,term[2][1])*value;
            od;
        od;
        return vector;
    end;
    result := rec(f:=f,g:=g,homotopy:=homotopy,character:=character,
        evaluateR:=evaluateR,pullback:=pullback,project:=pullback,
        normalizeSimplex:=anchor,
        convention:="normalized-homogeneous-bar-local-cochains",
        resolution:=R,dimension:=dimension);
    result.lift := function(n,vector,sign)
        return simplex -> evaluateR(n,vector,sign,simplex);
    end;
    result.homotopyValue := function(n,cochain,simplex)
        local term, value, answer;
        if not IsInt(n) or n<1 or not IsFunction(cochain) then
            Error("koAHSS: homotopyValue requires a positive cochain degree and function");
        fi;
        checkSimplex(simplex);
        if Length(simplex)<>n then Error("koAHSS: wrong simplex degree for bar homotopy"); fi;
        answer := 0;
        for term in homotopy(anchor(simplex)) do
            value := cochain(term[2]);
            if not IsInt(value) then Error("koAHSS: bar cochains must return integers"); fi;
            answer := answer+term[1]*value;
        od;
        return answer;
    end;
    return result;
end);
