# The arity-two native contraction and higher diagonals over F_2.
# Terms [p,j,g,q,k,h] use indices in R!.elts for both group elements.
# This is the reduction of the integral tensor recurrence, with its same
# augmentation and contraction. No cohomology representative is changed.
InstallGlobalFunction(KOAHSS_BinaryTensorEngine,function(arg)
    local R,memo,length,identity,index,product,reduce,newCache,remember,
        cells,cell,boundaries,cellBoundary,act,contract,boundary,
        diagonals,diagonal,cup,stats;
    R:=arg[1];
    if Length(arg)=2 then memo:=arg[2]; else memo:=KOAHSS_ResolutionMemo(R); fi;
    length:=EvaluateProperty(R,"length");
    index:=memo.index; identity:=memo.identity; product:=memo.product;
    stats:=rec(contractions:=0,contractionHits:=0,
        diagonalRequests:=0,diagonalHits:=0,largestTensor:=0);
    reduce:=function(terms)
        return List(Filtered(Collected(terms),t->IsOddInt(t[2])),t->t[1]);
    end;
    newCache:=function()
        return rec(values:=NewDictionary([],true),keys:=[],terms:=0);
    end;
    remember:=function(cache,k,value)
        local old;
        MakeImmutable(value);
        if Length(value)>200000 then return value; fi;
        while not IsEmpty(cache.keys) and
            (Length(cache.keys)>=2048 or cache.terms+Length(value)>200000) do
            old:=Remove(cache.keys,1);
            cache.terms:=cache.terms-Length(LookupDictionary(cache.values,old));
            RemoveDictionary(cache.values,old);
        od;
        AddDictionary(cache.values,k,value); Add(cache.keys,k);
        cache.terms:=cache.terms+Length(value);
        return value;
    end;
    cells:=newCache(); boundaries:=[]; diagonals:=newCache();
    cell:=function(n,j,g)
        local value;
        stats.contractions:=stats.contractions+1;
        value:=LookupDictionary(cells.values,[n,j,g]);
        if value<>fail then stats.contractionHits:=stats.contractionHits+1; return value; fi;
        # Positive basis indices, as in the integral contraction. Reduce
        # signed word multiplicities only after the complete word is read.
        value:=reduce(List(memo.contract(n,j,g),u->[AbsInt(u[1]),u[2]]));
        return remember(cells,[n,j,g],value);
    end;
    cellBoundary:=function(n,j)
        if not IsBound(boundaries[n+1]) then boundaries[n+1]:=[]; fi;
        if not IsBound(boundaries[n+1][j]) then
            boundaries[n+1][j]:=Immutable(reduce(List(R!.boundary(n,j),
                u->[AbsInt(u[1]),u[2]])));
        fi;
        return boundaries[n+1][j];
    end;
    act:=function(terms,g)
        if g=identity then return terms; fi;
        return List(terms,t->[t[1],t[2],product(g,t[3]),t[4],t[5],product(g,t[6])]);
    end;
    contract:=function(terms,cap)
        local result,t,u;
        result:=[];
        for t in terms do
            if t[1]<cap and t[4]<=cap then
                for u in cell(t[1],t[2],t[3]) do
                    Add(result,[t[1]+1,u[1],u[2],t[4],t[5],t[6]]);
                od;
            fi;
            if t[1]=0 and t[4]<cap then
                for u in cell(t[4],t[5],t[6]) do
                    Add(result,[0,1,identity,t[4]+1,u[1],u[2]]);
                od;
            fi;
        od;
        result:=reduce(result);
        stats.largestTensor:=Maximum(stats.largestTensor,Length(result));
        return result;
    end;
    boundary:=function(terms)
        local result,t,u;
        result:=[];
        for t in terms do
            if t[1]>0 then
                for u in cellBoundary(t[1],t[2]) do
                    Add(result,[t[1]-1,u[1],product(t[3],u[2]),t[4],t[5],t[6]]);
                od;
            fi;
            if t[4]>0 then
                for u in cellBoundary(t[4],t[5]) do
                    Add(result,[t[1],t[2],t[3],t[4]-1,u[1],product(t[6],u[2])]);
                od;
            fi;
        od;
        return reduce(result);
    end;
    diagonal:=function(i,n,j,cap)
        local value,rhs,t,previous;
        if i<0 or n<0 then return []; fi;
        if not IsInt(cap) or cap<0 or cap>length then Error("invalid binary tensor degree cap"); fi;
        stats.diagonalRequests:=stats.diagonalRequests+1;
        value:=LookupDictionary(diagonals.values,[cap,i,n,j]);
        if value<>fail then stats.diagonalHits:=stats.diagonalHits+1; return value; fi;
        if i=0 and n=0 then value:=[[0,j,identity,0,j,identity]];
        else
            rhs:=[];
            if n>0 then
                for t in cellBoundary(n,j) do
                    Append(rhs,act(diagonal(i,n-1,t[1],cap),t[2]));
                od;
            fi;
            if i>0 then
                previous:=diagonal(i-1,n,j,cap);
                Append(rhs,previous);
                Append(rhs,List(previous,t->[t[4],t[5],t[6],t[1],t[2],t[3]]));
            fi;
            value:=contract(reduce(rhs),cap);
        fi;
        return remember(diagonals,[cap,i,n,j],value);
    end;
    cup:=function(i,p,a,q,b)
        local n,result,j,t;
        n:=p+q-i;
        if n<0 then return []; fi;
        result:=List([1..R!.dimension(n)],j->0);
        if i<0 then return result; fi;
        if Length(a)<>R!.dimension(p) or Length(b)<>R!.dimension(q) then
            Error("cup input has wrong cochain dimension");
        fi;
        for j in [1..Length(result)] do
            for t in diagonal(i,n,j,Maximum(p,q)) do
                if t[1]=p and t[4]=q then result[j]:=result[j]+a[t[2]]*b[t[5]]; fi;
            od;
            result[j]:=result[j] mod 2;
        od;
        return result;
    end;
    return rec(index:=index,identity:=identity,reduce:=reduce,act:=act,
        contract:=contract,boundary:=boundary,diagonal:=diagonal,cup:=cup,stats:=stats);
end);
