# Interval-cut operations on an explicitly supplied simplicial cochain model.
# This module does not turn an arbitrary HAP diagonal into a surjection action.
# Its head-suspension chi family is an Adem primitive, not a ko normalization.
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

InstallGlobalFunction(koAHSSNaturalWordValue, function(arg)
    local word, degrees, inputs, simplex, integral, cupIndex, arity, degree,
          labels, j, last, occurrences, fixedSign, visit, result,
          collectOnly,patterns,cache,key,cacheable,pattern,value,entry;
    if not Length(arg) in [4,5] then
        Error("koAHSSNaturalWordValue(word,degrees,inputs,simplex[,integralCupIndex])");
    fi;
    word := arg[1]; degrees := arg[2]; inputs := arg[3]; simplex := arg[4];
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
    if not IsList(simplex) then Error("koAHSS: simplex vertices must be a list"); fi;
    integral := Length(arg) = 5; cupIndex := 0;
    if integral then
        cupIndex := arg[5];
        if not IsInt(cupIndex) or cupIndex < 0 or arity <> 2
           or labels <> List([0..cupIndex+1], j -> 1+j mod 2) then
            Error("koAHSS: integral mode requires the matching alternating binary cup word");
        fi;
    fi;
    degree := Sum(degrees) - Length(labels) + arity;
    if degree < 0 then return 0; fi;
    if Length(simplex) <> degree+1 then
        Error("koAHSS: simplex has the wrong degree for the word operation");
    fi;
    if degree < Maximum(degrees) then return 0; fi;
    occurrences := List([1..arity], j -> Number(labels, x -> x=j));
    if ForAny([1..arity], j -> occurrences[j] > degrees[j]+1) then return 0; fi;
    last := List([1..arity], j -> Last(Positions(labels,j)));
    fixedSign := cupIndex*Sum(degrees) + Binomial(cupIndex,2);
    result := 0; collectOnly:=false;
    cacheable:=arity<=4 and degree<=12 and Length(labels)<=12;
    patterns:=fail;
    if cacheable then
        cache:=KOAHSS_NATURAL_WORD_PATTERN_CACHE;
        key:=Concatenation(JoinStringsWithSeparator(List(labels,String),","),"/",
            JoinStringsWithSeparator(List(degrees,String),","),"/",String(integral));
        if integral then key:=Concatenation(key,"/",String(cupIndex)); fi;
        if IsBound(cache.entries.(key)) then patterns:=cache.entries.(key); fi;
    fi;
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
                    if integral then Add(patterns,[updated,(-1)^(fixedSign+sign)]);
                    else Add(patterns,[updated,1]); fi;
                elif integral then result := result + (-1)^(fixedSign+sign)*value;
                else result := (result+value) mod 2; fi;
            else
                visit(position+1, endpoint, updated,
                    Concatenation(masses,[mass]), sign, value);
            fi;
        od;
    end;
    if not cacheable then
        visit(1, 1, List([1..arity], j -> []), [], 0, 1);
        return result;
    fi;
    if patterns=fail then
        collectOnly:=true; patterns:=[];
        visit(1, 1, List([1..arity], j -> []), [], 0, 1);
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
end);

BindGlobal("KOAHSS_NATURAL_CHI_ANF_COMPILED",rec(head:=[],tail:=[]));

# Compile the finite ANF once into a prefix tree. At its last level a bit-list
# dot product replaces a loop over fourth factors. No cochain identities or
# cocycle assumptions are used in this optimization: a(face)^2=a(face) in F2.
BindGlobal("KOAHSS_NaturalChiANFValue",function(n,a,simplex,family)
    local allData,data,compiled,terms,code,term,radix,build,tree,evaluate,
          values,value,face,answer,globalName,filename,cacheRows;
    if family="head" then
        globalName:="KOAHSS_NATURAL_CHI_HEAD4_ANF";
        filename:="chi-head-degree4-anf.g";
    else
        globalName:="KOAHSS_NATURAL_CHI_CALIBRATED_ANF";
        filename:="chi-calibrated-degree7-anf.g";
    fi;
    if not IsBoundGlobal(globalName) then
        Read(Filename(KOAHSS_NATURAL_WORD_DATA_DIRECTORY,filename));
    fi;
    allData:=ValueGlobal(globalName);
    data:=allData.degrees[n+1];
    if IsEmpty(data.monomials) then return 0; fi;
    cacheRows:=KOAHSS_NATURAL_CHI_ANF_COMPILED.(family);
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
    local result, word, words,n,a,simplex,convention,maximum;
    if not Length(arg) in [3,4] then Error("koAHSSNaturalChiValue(n,a,simplex[,convention])"); fi;
    n:=arg[1]; a:=arg[2]; simplex:=arg[3]; convention:="tail";
    if Length(arg)=4 then convention:=arg[4]; fi;
    if not convention in ["head","tail","tail6","tail7"] then
        Error("chi convention must be head, tail, tail6, or tail7");
    fi;
    maximum:=6;
    if convention in ["tail","tail7"] then maximum:=7; fi;
    if not IsInt(n) or n < 0 or n > maximum then
        Error("koAHSS: chi input degree is outside the selected finite family");
    fi;
    if not IsFunction(a) or not IsList(simplex) or Length(simplex) <> n+4 then
        Error("koAHSS: chi needs a cochain function and an (n+3)-simplex");
    fi;
    if convention in ["tail","tail7"] then
        return KOAHSS_NaturalChiANFValue(n,a,simplex,"tail");
    fi;
    if convention="head" and n<=4 then
        return KOAHSS_NaturalChiANFValue(n,a,simplex,"head");
    fi;
    if not IsBoundGlobal("KOAHSS_NATURAL_CHI_HEAD_WORDS") then
        Read(Filename(KOAHSS_NATURAL_WORD_DATA_DIRECTORY, "chi-head-degree6.g"));
    fi;
    words := ValueGlobal("KOAHSS_NATURAL_CHI_HEAD_WORDS");
    if convention="tail6" and n<6 then
        if not IsBoundGlobal("KOAHSS_NATURAL_CHI_TAIL_WORDS") then
            Read(Filename(KOAHSS_NATURAL_WORD_DATA_DIRECTORY,"chi-tail-degree6.g"));
        fi;
        words:=ValueGlobal("KOAHSS_NATURAL_CHI_TAIL_WORDS");
    fi;
    result := 0;
    for word in words[n+1] do
        result := (result + koAHSSNaturalWordValue(word, [n,n,n,n],
            [a,a,a,a], simplex)) mod 2;
    od;
    return result;
end);

InstallGlobalFunction(koAHSSNaturalZetaValue, function(kind, n, x, y, simplex)
    local word, j;
    if not kind in [1,2] or not IsInt(n) or n < 0 then
        Error("koAHSS: Cartan word helpers require kind 1 or 2 and nonnegative input degree");
    fi;
    if n=0 then
        if not IsFunction(x) or not IsFunction(y) or not IsList(simplex)
           or Length(simplex)<>kind+2 then Error("invalid degree-zero Cartan input"); fi;
        return 0;
    fi;
    if kind = 1 then
        word := [1,2,3,2];
        Append(word, List([0..n-1], j -> 4-j mod 2));
    else
        word := [1,2,3,1];
        Append(word, List([0..n], j -> 3+j mod 2));
    fi;
    return koAHSSNaturalWordValue(word, [kind,kind,n,n], [x,x,y,y], simplex);
end);
