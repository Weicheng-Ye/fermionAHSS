# Interval-cut operations on an explicitly supplied simplicial cochain model.
# This module does not turn an arbitrary HAP diagonal into a surjection action.
# Its chi values are the calibrated chi7_tail ANF in input degrees 0-7.
CallFuncList(function()
    local name, slashes, directory;
    name := INPUT_FILENAME();
    if not IsEmpty(name) and name[1] <> '/' then
        name := Filename(DirectoryCurrent(),name);
    fi;
    slashes := Positions(name, '/');
    if IsEmpty(slashes) then directory := Directory(".");
    else directory := Directory(name{[1..Last(slashes)]}); fi;
    BindGlobal("KOAHSS_NATURAL_WORD_DATA_DIRECTORY",
        Directory(Concatenation(Filename(directory, ""), "../data/")));
end, []);

BindGlobal("KOAHSS_NATURAL_CHI_CONVENTION", "chi-tail-suspension-degree7");

# Cache only small interval-cut shapes, independently of all cochain values.
# Raw large Adem tables retain the streaming evaluator and cannot fill this.
BindGlobal("KOAHSS_NATURAL_WORD_PATTERN_CACHE",rec(
    entries:=rec(),count:=0,patternCount:=0));

# The interval cuts of a word depend only on the word, the input degrees and
# the mode. An evaluator checks these and collects the cuts once; the function
# it returns only multiplies input values on each simplex.
BindGlobal("KOAHSS_NaturalWordEvaluator", function(arg)
    local word, degrees, inputs, integral, cupIndex, arity, degree, labels,
          last, occurrences, fixedSign, checked, cuts, evaluate, cacheable,
          cache, key, patterns;
    if not Length(arg) in [3,4] then
        Error("KOAHSS_NaturalWordEvaluator(word,degrees,inputs[,integralCupIndex])");
    fi;
    word := arg[1]; degrees := arg[2]; inputs := arg[3];
    if IsString(word) then
        if IsEmpty(word) or not ForAll(word, c -> c in "123456789") then
            Error("koAHSS: a natural word must contain positive integer labels");
        fi;
        labels := List(word, c -> Int([c]));
    elif IsList(word) and not IsEmpty(word)
         and ForAll(word, x -> IsInt(x) and x > 0) then
        labels := ShallowCopy(word);
    else Error("koAHSS: a natural word must contain positive integer labels"); fi;
    if not IsList(degrees) or IsEmpty(degrees)
       or not ForAll(degrees, x -> IsInt(x) and x >= 0) then
        Error("koAHSS: word input degrees must be nonnegative integers");
    fi;
    arity := Length(degrees);
    if Set(labels) <> [1..arity]
       or ForAny([2..Length(labels)], j -> labels[j] = labels[j-1]) then
        Error("koAHSS: a natural word must be a nondegenerate surjection");
    fi;
    if not IsList(inputs) or Length(inputs) <> arity
       or not ForAll(inputs, IsFunction) then
        Error("koAHSS: word inputs must be one cochain function per degree");
    fi;
    integral := Length(arg) = 4; cupIndex := 0;
    if integral then
        cupIndex := arg[4];
        if not IsInt(cupIndex) or cupIndex < 0 or arity <> 2
           or labels <> List([0..cupIndex+1], j -> 1+j mod 2) then
            Error("koAHSS: integral mode requires the matching alternating binary cup word");
        fi;
    fi;
    degree := Sum(degrees) - Length(labels) + arity;
    if degree < 0 then
        return function(simplex)
            if not IsList(simplex) then Error("koAHSS: simplex vertices must be a list"); fi;
            return 0;
        end;
    fi;
    checked := function(simplex)
        if not IsList(simplex) then Error("koAHSS: simplex vertices must be a list"); fi;
        if Length(simplex) <> degree+1 then
            Error("koAHSS: simplex has the wrong degree for the word operation");
        fi;
        return simplex;
    end;
    occurrences := List([1..arity], j -> Number(labels, x -> x=j));
    if degree < Maximum(degrees)
       or ForAny([1..arity], j -> occurrences[j] > degrees[j]+1) then
        return function(simplex) checked(simplex); return 0; end;
    fi;
    last := List([1..arity], j -> Last(Positions(labels,j)));
    fixedSign := cupIndex*Sum(degrees) + Binomial(cupIndex,2);
    # Enumerate the cuts: collect them as [faces,sign] when simplex=fail,
    # otherwise evaluate the word on the simplex while enumerating.
    cuts := function(simplex)
        local collectOnly, result, found, visit;
        collectOnly := simplex = fail; result := 0; found := [];
        # Positions here are one-based. The note's cut i_l is endpoint-1.
        # Repeated vertices within one input cannot occur in any retained cut:
        # their multiplicity would make the total face size too small. Rejecting
        # them early is therefore equivalent to the note's increasing-union rule.
        visit := function(position, previous, faces, masses, parity, product)
            local label, ends, endpoint, face, updated, mass, sign, k, value;
            label := labels[position];
            if position = Length(labels) then ends := [degree+1];
            else ends := [previous..degree+1]; fi;
            for endpoint in ends do
                if not IsEmpty(faces[label]) and Last(faces[label]) >= previous then
                    continue;
                fi;
                face := Concatenation(faces[label], [previous..endpoint]);
                if Length(face) > degrees[label]+1 then continue; fi;
                if position = last[label] and Length(face) <> degrees[label]+1 then
                    continue;
                fi;
                updated := ShallowCopy(faces); updated[label] := face;
                value := product;
                if position = last[label] and not collectOnly then
                    k := inputs[label](simplex{face});
                    if not IsInt(k) then Error("koAHSS: word cochains must return integers"); fi;
                    if not integral then k := k mod 2; fi;
                    value := value*k;
                    if value = 0 then continue; fi;
                fi;
                sign := parity; mass := endpoint-previous;
                if position < last[label] then
                    mass := mass+1; sign := sign+endpoint-1;
                fi;
                if integral then
                    for k in [1..position-1] do
                        if labels[k] > label then sign := sign + masses[k]*mass; fi;
                    od;
                fi;
                if position = Length(labels) then
                    if collectOnly then
                        if integral then Add(found,[updated,(-1)^(fixedSign+sign)]);
                        else Add(found,[updated,1]); fi;
                    elif integral then result := result + (-1)^(fixedSign+sign)*value;
                    else result := (result+value) mod 2; fi;
                else
                    visit(position+1, endpoint, updated,
                        Concatenation(masses,[mass]), sign, value);
                fi;
            od;
        end;
        visit(1, 1, List([1..arity], j -> []), [], 0, 1);
        if collectOnly then return found; fi;
        return result;
    end;
    evaluate := function(patterns, simplex)
        local result, pattern, value, j, entry;
        result := 0;
        for pattern in patterns do
            value:=pattern[2];
            for j in [1..arity] do
                entry:=inputs[j](simplex{pattern[1][j]});
                if not IsInt(entry) then Error("koAHSS: word cochains must return integers"); fi;
                if not integral then entry:=entry mod 2; fi;
                value:=value*entry;
                if value=0 then break; fi;
            od;
            result:=result+value;
        od;
        if not integral then result:=result mod 2; fi;
        return result;
    end;
    cacheable:=arity<=4 and degree<=12 and Length(labels)<=12;
    if not cacheable then return simplex -> cuts(checked(simplex)); fi;
    cache:=KOAHSS_NATURAL_WORD_PATTERN_CACHE;
    key:=Concatenation(JoinStringsWithSeparator(List(labels,String),","),"/",
        JoinStringsWithSeparator(List(degrees,String),","),"/",String(integral));
    if integral then key:=Concatenation(key,"/",String(cupIndex)); fi;
    if IsBound(cache.entries.(key)) then patterns:=cache.entries.(key);
    else
        patterns:=cuts(fail);
        # Bound both the number of shapes and their retained interval cuts.
        if Length(patterns)<=4096 then
            if cache.count>=512 or cache.patternCount+Length(patterns)>65536 then
                cache.entries:=rec();cache.count:=0;cache.patternCount:=0;
            fi;
            MakeImmutable(patterns);cache.entries.(key):=patterns;
            cache.count:=cache.count+1;
            cache.patternCount:=cache.patternCount+Length(patterns);
        fi;
    fi;
    # Larger shapes are not retained; their cuts are collected on each call.
    if Length(patterns)>4096 then
        return simplex -> evaluate(cuts(fail), checked(simplex));
    fi;
    return simplex -> evaluate(patterns, checked(simplex));
end);

InstallGlobalFunction(koAHSSNaturalWordValue, function(arg)
    if not Length(arg) in [4,5] then
        Error("koAHSSNaturalWordValue(word,degrees,inputs,simplex[,integralCupIndex])");
    fi;
    return CallFuncList(KOAHSS_NaturalWordEvaluator,
        arg{Difference([1..Length(arg)],[4])})(arg[4]);
end);

BindGlobal("KOAHSS_NATURAL_CHI_ANF_COMPILED",[]);

# Compile the finite ANF once into a prefix tree. At its last level a bit-list
# dot product replaces a loop over fourth factors. No cochain identities or
# cocycle assumptions are used in this optimization: a(face)^2=a(face) in F2.
BindGlobal("KOAHSS_NaturalChiANFValue",function(n,a,simplex)
    local allData,data,compiled,terms,code,term,radix,build,tree,evaluate,
          values,value,face,answer,cacheRows;
    if not IsBoundGlobal("KOAHSS_NATURAL_CHI_CALIBRATED_ANF") then
        Read(Filename(KOAHSS_NATURAL_WORD_DATA_DIRECTORY,"chi-calibrated-degree7-anf.g"));
    fi;
    allData:=ValueGlobal("KOAHSS_NATURAL_CHI_CALIBRATED_ANF");
    data:=allData.degrees[n+1];
    if IsEmpty(data.monomials) then return 0; fi;
    cacheRows:=KOAHSS_NATURAL_CHI_ANF_COMPILED;
    if not IsBound(cacheRows[n+1]) then
        radix:=2^data.bitsPerFace; terms:=[];
        for code in data.monomials do
            term:=[];
            while code>0 do
                Add(term,code mod radix); code:=QuoInt(code,radix);
            od;
            Add(terms,term);
        od;
        Sort(terms);
        build:=function(rows,depth)
            local constant,position,label,indices,children,group;
            constant:=0; position:=1;
            if not IsEmpty(rows) and IsEmpty(rows[1]) then
                constant:=1; position:=2;
            fi;
            if depth=3 then
                indices:=List(rows{[position..Length(rows)]},row->row[1]);
                return [constant,BlistList([1..Length(data.faces)],indices)];
            fi;
            indices:=[]; children:=[];
            while position<=Length(rows) do
                label:=rows[position][1]; group:=[];
                while position<=Length(rows) and rows[position][1]=label do
                    Add(group,rows[position]{[2..Length(rows[position])]});
                    position:=position+1;
                od;
                Add(indices,label); Add(children,build(group,depth+1));
            od;
            return [constant,indices,children];
        end;
        tree:=build(terms,0); MakeImmutable(tree);
        compiled:=rec(tree:=tree,cache:=NewDictionary([],true),cacheCount:=0);
        cacheRows[n+1]:=compiled;
    fi;
    compiled:=cacheRows[n+1]; values:=[];
    for face in data.faces do
        value:=a(simplex{face});
        if not IsInt(value) then Error("koAHSS: chi cochains must return integers"); fi;
        Add(values,value mod 2=1);
    od;
    answer:=LookupDictionary(compiled.cache,values);
    if answer<>fail then return answer; fi;
    evaluate:=function(node)
        local result,j;
        result:=node[1];
        if Length(node)=2 then
            return result+SizeBlist(IntersectionBlist(node[2],values));
        fi;
        for j in [1..Length(node[2])] do
            if values[node[2][j]] then result:=result+evaluate(node[3][j]); fi;
        od;
        return result;
    end;
    answer:=evaluate(compiled.tree) mod 2;
    # This cache depends only on the finite input-face bit vector, and is
    # consequently reusable across simplices and cochain callbacks. Bound its
    # size so long group batches do not retain every past universal input.
    if compiled.cacheCount>=8192 then
        compiled.cache:=NewDictionary([],true); compiled.cacheCount:=0;
    fi;
    MakeImmutable(values); AddDictionary(compiled.cache,values,answer);
    compiled.cacheCount:=compiled.cacheCount+1;
    return answer;
end);

InstallGlobalFunction(koAHSSNaturalChiValue, function(arg)
    local n,a,simplex;
    if not Length(arg) in [3,4] then Error("koAHSSNaturalChiValue(n,a,simplex[,convention])"); fi;
    n:=arg[1]; a:=arg[2]; simplex:=arg[3];
    if Length(arg)=4 and not arg[4] in ["tail","tail7"] then
        Error("koAHSS: the chi convention is the calibrated chi7_tail (tail)");
    fi;
    if not IsInt(n) or n < 0 or n > 7 then
        Error("koAHSS: chi input degree is outside the selected finite family");
    fi;
    if not IsFunction(a) or not IsList(simplex) or Length(simplex) <> n+4 then
        Error("koAHSS: chi needs a cochain function and an (n+3)-simplex");
    fi;
    return KOAHSS_NaturalChiANFValue(n,a,simplex);
end);

# The Cartan word of zeta_kind in input degree n>0, as a word evaluator.
BindGlobal("KOAHSS_NaturalZetaEvaluator", function(kind, n, x, y)
    local word;
    if kind = 1 then
        word := [1,2,3,2];
        Append(word, List([0..n-1], j -> 4-j mod 2));
    else
        word := [1,2,3,1];
        Append(word, List([0..n], j -> 3+j mod 2));
    fi;
    return KOAHSS_NaturalWordEvaluator(word, [kind,kind,n,n], [x,x,y,y]);
end);

InstallGlobalFunction(koAHSSNaturalZetaValue, function(kind, n, x, y, simplex)
    if not kind in [1,2] or not IsInt(n) or n < 0 then
        Error("koAHSS: Cartan word helpers require kind 1 or 2 and nonnegative input degree");
    fi;
    if n=0 then
        if not IsFunction(x) or not IsFunction(y) or not IsList(simplex)
           or Length(simplex)<>kind+2 then Error("invalid degree-zero Cartan input"); fi;
        return 0;
    fi;
    return KOAHSS_NaturalZetaEvaluator(kind, n, x, y)(simplex);
end);
