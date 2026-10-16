# The arity-two native contraction and higher diagonals over F_2.
# A term [p,j,g,q,k,h] (degrees p,q, basis indices j,k, and the indices g,h
# in R!.elts of both group elements) is kept as one integer in a mixed radix
# (2^8 for a degree, 2^20 for a basis index, 2^40 for an element index, as
# gap/cell_transport.gi keys its chains): comparing and sorting terms then
# runs in the kernel, a tensor is reduced over F_2 by one sort and a parity
# scan, and the caches are keyed by single integers.
# The cache budgets (entries, retained terms, largest stored value) are the
# fields cellCache and diagonalCache of KOAHSS_EXTENSION_TRANSPORT_LIMITS.
# This is the reduction of the integral tensor recurrence, with its same
# augmentation and contraction. No cohomology representative is changed.
InstallGlobalFunction(KOAHSS_BinaryTensorEngine,function(arg)
    local R,memo,length,identity,index,product,BQ,BJ,BG,BL,BGJ,limits,element,basis,
        pack,unpack,swap,reduce,newCache,remember,cells,cell,boundaries,cellBoundary,
        act,contract,boundary,diagonals,diagonal,cup,stats;
    R:=arg[1];
    if Length(arg)=2 then memo:=arg[2]; else memo:=KOAHSS_ResolutionMemo(R); fi;
    length:=EvaluateProperty(R,"length");
    index:=memo.index; identity:=memo.identity; product:=memo.product;
    BQ:=2^8; BJ:=2^20; BG:=2^40; BGJ:=BG*BJ; BL:=BQ*BGJ;
    limits:=ValueGlobal("KOAHSS_EXTENSION_TRANSPORT_LIMITS");
    stats:=rec(contractions:=0,contractionHits:=0,
        diagonalRequests:=0,diagonalHits:=0,largestTensor:=0,
        refusedValues:=0,evictions:=0);
    # Element indices grow with R!.elts; every digit must stay below its base.
    element:=function(x)
        if x>=BG then Error("binary tensor engine: an element index exceeds 2^40"); fi;
        return x;
    end;
    basis:=function(x)
        if x>=BJ then Error("binary tensor engine: a basis index exceeds 2^20"); fi;
        return x;
    end;
    # term = ((p*BJ+j)*BG+g)*BL + (q*BJ+k)*BG+h
    pack:=function(p,j,g,q,k,h)
        if p>=BQ or q>=BQ then Error("binary tensor engine: a degree exceeds 2^8"); fi;
        return ((p*BJ+basis(j))*BG+element(g))*BL+(q*BJ+basis(k))*BG+element(h);
    end;
    unpack:=function(t)
        local low,hi;
        low:=t mod BL; hi:=QuoInt(t,BL);
        return [QuoInt(hi,BGJ),QuoInt(hi,BG) mod BJ,hi mod BG,
                QuoInt(low,BGJ),QuoInt(low,BG) mod BJ,low mod BG];
    end;
    # [p,j,g,q,k,h] -> [q,k,h,p,j,g]
    swap:=t->(t mod BL)*BL+QuoInt(t,BL);
    # The terms of odd multiplicity, sorted: one kernel sort and a parity scan.
    reduce:=function(terms)
        local out,i,n,current,count;
        if not IsMutable(terms) then terms:=ShallowCopy(terms); fi;
        Sort(terms); out:=[]; i:=1; n:=Length(terms);
        while i<=n do
            current:=terms[i]; count:=1; i:=i+1;
            while i<=n and terms[i]=current do count:=count+1; i:=i+1; od;
            if IsOddInt(count) then Add(out,current); fi;
        od;
        return out;
    end;
    newCache:=function(budget)
        return rec(values:=NewDictionary(1,true),keys:=[],terms:=0,
            entries:=budget[1],retained:=budget[2],largest:=budget[3]);
    end;
    remember:=function(cache,k,value)
        local old;
        MakeImmutable(value);
        if Length(value)>cache.largest then
            stats.refusedValues:=stats.refusedValues+1; return value;
        fi;
        while not IsEmpty(cache.keys) and
            (Length(cache.keys)>=cache.entries or cache.terms+Length(value)>cache.retained) do
            old:=Remove(cache.keys,1);
            cache.terms:=cache.terms-Length(LookupDictionary(cache.values,old));
            RemoveDictionary(cache.values,old); stats.evictions:=stats.evictions+1;
        od;
        AddDictionary(cache.values,k,value); Add(cache.keys,k);
        cache.terms:=cache.terms+Length(value);
        return value;
    end;
    cells:=newCache(limits.cellCache); boundaries:=[]; diagonals:=newCache(limits.diagonalCache);
    # A contraction or boundary word is a list of j'*BG+g' (basis index, element index).
    cell:=function(n,j,g)
        local value,key;
        stats.contractions:=stats.contractions+1;
        key:=(n*BJ+j)*BG+g;
        value:=LookupDictionary(cells.values,key);
        if value<>fail then stats.contractionHits:=stats.contractionHits+1; return value; fi;
        # Positive basis indices, as in the integral contraction. Reduce
        # signed word multiplicities only after the complete word is read.
        value:=reduce(List(memo.contract(n,j,g),u->basis(AbsInt(u[1]))*BG+element(u[2])));
        return remember(cells,key,value);
    end;
    cellBoundary:=function(n,j)
        if not IsBound(boundaries[n+1]) then boundaries[n+1]:=[]; fi;
        if not IsBound(boundaries[n+1][j]) then
            boundaries[n+1][j]:=Immutable(reduce(List(R!.boundary(n,j),
                u->basis(AbsInt(u[1]))*BG+element(u[2]))));
        fi;
        return boundaries[n+1][j];
    end;
    act:=function(terms,g)
        local out,t,low,hi;
        if g=identity then return terms; fi;
        out:=[];
        for t in terms do
            low:=t mod BL; hi:=QuoInt(t,BL);
            Add(out,(QuoInt(hi,BG)*BG+element(product(g,hi mod BG)))*BL
                +QuoInt(low,BG)*BG+element(product(g,low mod BG)));
        od;
        return out;
    end;
    contract:=function(terms,cap)
        local result,t,u,low,hi,h,k,q,g,j,p;
        result:=[];
        for t in terms do
            low:=t mod BL; hi:=QuoInt(t,BL);
            h:=low mod BG; k:=QuoInt(low,BG) mod BJ; q:=QuoInt(low,BGJ);
            g:=hi mod BG; j:=QuoInt(hi,BG) mod BJ; p:=QuoInt(hi,BGJ);
            if p<cap and q<=cap then
                for u in cell(p,j,g) do
                    # [p+1, u1, u2, q, k, h]
                    Add(result,(((p+1)*BJ+QuoInt(u,BG))*BG+(u mod BG))*BL+low);
                od;
            fi;
            if p=0 and q<cap then
                for u in cell(q,k,h) do
                    # [0, 1, identity, q+1, u1, u2]
                    Add(result,(BG+identity)*BL+((q+1)*BJ+QuoInt(u,BG))*BG+(u mod BG));
                od;
            fi;
        od;
        result:=reduce(result);
        stats.largestTensor:=Maximum(stats.largestTensor,Length(result));
        return result;
    end;
    boundary:=function(terms)
        local result,t,u,low,hi,h,k,q,g,j,p;
        result:=[];
        for t in terms do
            low:=t mod BL; hi:=QuoInt(t,BL);
            h:=low mod BG; k:=QuoInt(low,BG) mod BJ; q:=QuoInt(low,BGJ);
            g:=hi mod BG; j:=QuoInt(hi,BG) mod BJ; p:=QuoInt(hi,BGJ);
            if p>0 then
                for u in cellBoundary(p,j) do
                    # [p-1, u1, g*u2, q, k, h]
                    Add(result,(((p-1)*BJ+QuoInt(u,BG))*BG+element(product(g,u mod BG)))*BL+low);
                od;
            fi;
            if q>0 then
                for u in cellBoundary(q,k) do
                    # [p, j, g, q-1, u1, h*u2]
                    Add(result,hi*BL+((q-1)*BJ+QuoInt(u,BG))*BG+element(product(h,u mod BG)));
                od;
            fi;
        od;
        return reduce(result);
    end;
    diagonal:=function(i,n,j,cap)
        local value,rhs,t,previous,key;
        if i<0 or n<0 then return []; fi;
        if not IsInt(cap) or cap<0 or cap>length then Error("invalid binary tensor degree cap"); fi;
        stats.diagonalRequests:=stats.diagonalRequests+1;
        key:=((cap*BQ+i)*BQ+n)*BJ+j;
        value:=LookupDictionary(diagonals.values,key);
        if value<>fail then stats.diagonalHits:=stats.diagonalHits+1; return value; fi;
        if i=0 and n=0 then value:=[pack(0,j,identity,0,j,identity)];
        else
            rhs:=[];
            if n>0 then
                for t in cellBoundary(n,j) do
                    Append(rhs,act(diagonal(i,n-1,QuoInt(t,BG),cap),t mod BG));
                od;
            fi;
            if i>0 then
                previous:=diagonal(i-1,n,j,cap);
                Append(rhs,previous);
                Append(rhs,List(previous,swap));
            fi;
            value:=contract(reduce(rhs),cap);
        fi;
        return remember(diagonals,key,value);
    end;
    cup:=function(i,p,a,q,b)
        local n,result,j,t,low,hi,tp,tq,tj,tk;
        n:=p+q-i;
        if n<0 then return []; fi;
        result:=List([1..R!.dimension(n)],j->0);
        if i<0 then return result; fi;
        if Length(a)<>R!.dimension(p) or Length(b)<>R!.dimension(q) then
            Error("cup input has wrong cochain dimension");
        fi;
        for j in [1..Length(result)] do
            for t in diagonal(i,n,j,Maximum(p,q)) do
                low:=t mod BL; hi:=QuoInt(t,BL);
                tq:=QuoInt(low,BGJ); tp:=QuoInt(hi,BGJ);
                if tp=p and tq=q then
                    tk:=QuoInt(low,BG) mod BJ; tj:=QuoInt(hi,BG) mod BJ;
                    result[j]:=result[j]+a[tj]*b[tk];
                fi;
            od;
            result[j]:=result[j] mod 2;
        od;
        return result;
    end;
    return rec(index:=index,identity:=identity,reduce:=reduce,act:=act,
        contract:=contract,boundary:=boundary,diagonal:=diagonal,cup:=cup,stats:=stats,
        pack:=pack,unpack:=unpack,swap:=swap,digits:=rec(Q:=BQ,J:=BJ,G:=BG,L:=BL));
end);
