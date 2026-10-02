# Exact, bounded memoization of the supplied resolution's contraction.
# R!.elts is append-only and may be extended by any of its consumers.
InstallGlobalFunction(KOAHSS_ResolutionMemo,function(R)
    local indices,synced,key,index,identity,words,order,retained,contract,
        products,productCount,product;
    indices:=NewDictionary([],true); synced:=0;
    key:=function(g)
        if IsPcpElementRep(g) then return g!.exponents; fi;
        return [g];
    end;
    index:=function(g)
        local k,position;
        # HAP lazy lists can be infinite and provide their own inverse
        # indexing map. Enumerating them would never terminate.
        if IsPseudoList(R!.elts) then
            position:=Position(R!.elts,g);
            if position=fail then Error("could not index resolution group element"); fi;
            return position;
        fi;
        while synced<Length(R!.elts) do
            synced:=synced+1; k:=key(R!.elts[synced]);
            if LookupDictionary(indices,k)=fail then AddDictionary(indices,k,synced); fi;
        od;
        position:=LookupDictionary(indices,key(g));
        if position=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(g); else Add(R!.elts,g); fi;
            while synced<Length(R!.elts) do
                synced:=synced+1; k:=key(R!.elts[synced]);
                if LookupDictionary(indices,k)=fail then AddDictionary(indices,k,synced); fi;
            od;
            position:=LookupDictionary(indices,key(g));
            if position=fail then Error("could not index resolution group element"); fi;
        fi;
        return position;
    end;
    identity:=index(One(R!.group));
    words:=NewDictionary([0,1,1],true); order:=[]; retained:=0;
    contract:=function(n,j,g)
        local k,word,old;
        k:=[n,j,g]; word:=LookupDictionary(words,k);
        if word<>fail then return word; fi;
        # Copy before making immutable: a HAP implementation may own this
        # list. Preserve its signs, order, and multiplicities exactly.
        word:=Immutable(R!.homotopy(n,[j,g]));
        if Length(word)>200000 then return word; fi;
        while not IsEmpty(order) and
            (Length(order)>=2048 or retained+Length(word)>200000) do
            old:=Remove(order,1);
            retained:=retained-Length(LookupDictionary(words,old));
            RemoveDictionary(words,old);
        od;
        AddDictionary(words,k,word); Add(order,k); retained:=retained+Length(word);
        return word;
    end;
    products:=NewDictionary([1,1],true); productCount:=0;
    product:=function(g,h)
        local value;
        if g=identity then return h; fi;
        if h=identity then return g; fi;
        value:=LookupDictionary(products,[g,h]);
        if value<>fail then return value; fi;
        value:=index(R!.elts[g]*R!.elts[h]);
        if productCount>=65536 then products:=NewDictionary([1,1],true); productCount:=0; fi;
        AddDictionary(products,[g,h],value); productCount:=productCount+1;
        return value;
    end;
    return rec(index:=index,identity:=identity,contract:=contract,product:=product);
end);
