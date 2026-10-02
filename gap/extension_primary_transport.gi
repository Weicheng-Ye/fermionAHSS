# Binary comparison of native and simplicial higher diagonals. This keeps
# the fixed f and uses only the tensor contraction of R, never g or bar H.
# Tensor coefficients and cancellations are binary; group actions are kept
# until the tensor has been contracted. See light-transport-reduction (12).
BindGlobal("KOAHSS_ExtensionPrimaryTransport",function(tr,native,maximum)
    local R,cache,order,retained,failure,cuts,reduce,add,index,vertexGroup,
        comparison,compute,tensorF,used;
    R:=tr.resolution; cache:=NewDictionary("",true); order:=[]; retained:=0;
    failure:=false; cuts:=[]; used:=0;
    vertexGroup:=v->v;
    if IsBound(tr.vertexGroup) then vertexGroup:=tr.vertexGroup; fi;
    reduce:=function(terms)
        local value;
        used:=used+Length(terms);
        if used>maximum then failure:=true; return []; fi;
        value:=native.reduce(terms);
        value:=List(Filtered(value,t->IsOddInt(t[1])),t->[1,t[2]]);
        if Length(value)>maximum then failure:=true; return []; fi;
        return value;
    end;
    add:=function(a,b) return reduce(Concatenation(a,b)); end;
    index:=function(g)
        local j;
        j:=Position(R!.elts,g);
        if j=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(g); else Add(R!.elts,g); fi;
            j:=Position(R!.elts,g);
        fi;
        return j;
    end;
    tensorF:=function(left,right)
        local out,a,b;
        out:=[];
        for a in tr.f(left) do
            if IsEvenInt(a[3]) then continue; fi;
            for b in tr.f(right) do
                if IsOddInt(b[3]) then
                    Add(out,[1,[[Length(left)-1,a[1],index(a[2])],
                                [Length(right)-1,b[1],index(b[2])]]]);
                fi;
            od;
        od;
        return reduce(out);
    end;
    comparison:=function(i,sigma)
        local anchor,key,value,old;
        if i<0 or ForAny([2..Length(sigma)],j->sigma[j]=sigma[j-1]) then return []; fi;
        anchor:=tr.normalizeSimplex(sigma);
        key:=Concatenation(String(i),":",String(anchor));
        value:=LookupDictionary(cache,key);
        if value=fail then
            value:=compute(i,anchor);
            if failure then return []; fi;
            MakeImmutable(value);
            if Length(value)<=200000 then
                while not IsEmpty(order) and (Length(order)>=256 or retained+Length(value)>200000) do
                    old:=Remove(order,1); retained:=retained-Length(LookupDictionary(cache,old));
                    RemoveDictionary(cache,old);
                od;
                AddDictionary(cache,key,value); Add(order,key); retained:=retained+Length(value);
            fi;
        fi;
        return reduce(native.act(value,vertexGroup(sigma[1])));
    end;
    compute:=function(i,sigma)
        local n,rhs,t,faces,lower,j,value;
        n:=Length(sigma)-1; rhs:=[];
        if n+i+1>native.maxTotalDegree then
            Error("primary comparison exceeds the supplied tensor degrees");
        fi;
        for t in tr.f(sigma) do
            if IsOddInt(t[3]) then
                rhs:=add(rhs,native.act(native.diagonal(i).evaluate(n,t[1]),t[2]));
            fi;
        od;
        for faces in cuts[i+1][n+1] do
            rhs:=add(rhs,tensorF(sigma{faces[1]},sigma{faces[2]}));
        od;
        if i>0 then
            lower:=comparison(i-1,sigma);
            rhs:=add(rhs,add(lower,native.permute(lower,[2,1])));
        fi;
        if n>0 then
            for j in [1..n+1] do
                rhs:=add(rhs,comparison(i,sigma{Filtered([1..n+1],q->q<>j)}));
            od;
        fi;
        if failure then return []; fi;
        value:=reduce(native.boundary(rhs));
        if failure then return []; fi;
        if value<>[] then Error("primary comparison source is not closed"); fi;
        value:=reduce(native.contract(rhs));
        if failure then return []; fi;
        lower:=reduce(native.boundary(value));
        if failure then return []; fi;
        if lower<>rhs then Error("primary comparison contraction failed"); fi;
        return value;
    end;
    return function(i,sigma,patterns)
        local value;
        # The cuts come from the one interval-cut engine, cochain_tools.py.
        cuts:=patterns; used:=0;
        value:=comparison(i,sigma);
        if failure then
            return rec(status:="unresolved",code:="primary-comparison-term-budget",
                reason:="native primary comparison exceeded its tensor term budget");
        fi;
        return rec(status:="computed",terms:=List(value,t->
            [1,t[2][1][1],t[2][1][2]-1,t[2][2][1],t[2][2][2]-1]));
    end;
end);
