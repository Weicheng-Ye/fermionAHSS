# Binary comparison of native and simplicial higher diagonals. This keeps
# the fixed f and uses only the tensor contraction of R, never g or bar H.
# Tensor coefficients and cancellations are binary; group actions are kept
# until the tensor has been contracted. See light-transport-reduction (12).
# The same recurrence in flat binary tensors. Keep the full tensors for the
# two exact boundary checks; collect only the final trivial-character pairing.
BindGlobal("KOAHSS_ExtensionPrimaryTransport",function(tr,native,maximum)
    local R,tensor,cache,order,retained,failure,cuts,used,extend,account,
        vertexGroup,labels,nextLabel,label,comparison,compute,tensorF,indexed,elementIndex,firstIndex;
    if not IsBound(native.binaryTensor) then
        Error("koFull: the shared binary tensor engine of the resolution is required");
    fi;
    R:=tr.resolution; tensor:=native.binaryTensor();
    cache:=NewDictionary([],true); order:=[]; retained:=0;
    failure:=false; cuts:=[]; used:=0; labels:=fail; nextLabel:=0;
    vertexGroup:=v->v;
    if IsBound(tr.vertexGroup) then vertexGroup:=tr.vertexGroup; fi;
    # The cell comparison carries element indices (the positions in R!.elts
    # that the tensor engine's index returns); the bar carries elements.
    indexed:=IsBound(tr.indexedElements) and tr.indexedElements=true;
    if indexed then
        elementIndex:=x->x; firstIndex:=sigma->tr.vertexIndex(sigma[1]);
    else
        elementIndex:=tensor.index; firstIndex:=sigma->tensor.index(vertexGroup(sigma[1]));
    fi;
    label:=function(v)
        local value;
        if labels=fail then labels:=NewDictionary(v,true); fi;
        value:=LookupDictionary(labels,v);
        if value=fail then
            nextLabel:=nextLabel+1; value:=nextLabel;
            AddDictionary(labels,Immutable(v),value);
        fi;
        return value;
    end;
    account:=function(value)
        used:=used+Length(value);
        if used>maximum then failure:=true; return []; fi;
        return value;
    end;
    extend:=function(target,source)
        Append(target,account(source));
    end;
    tensorF:=function(left,right)
        local out,a,b,x,y;
        out:=[];
        x:=Filtered(tr.f(left),a->IsOddInt(a[3]));
        y:=Filtered(tr.f(right),b->IsOddInt(b[3]));
        for a in x do
            for b in y do
                Add(out,[Length(left)-1,a[1],elementIndex(a[2]),
                         Length(right)-1,b[1],elementIndex(b[2])]);
            od;
        od;
        return tensor.reduce(out);
    end;
    comparison:=function(i,sigma)
        local anchor,key,value,old;
        if i<0 or ForAny([2..Length(sigma)],j->sigma[j]=sigma[j-1]) then return []; fi;
        anchor:=tr.normalizeSimplex(sigma);
        key:=Concatenation([i],List(anchor,label));
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
        return tensor.act(value,firstIndex(sigma));
    end;
    compute:=function(i,sigma)
        local n,rhs,t,faces,lower,j,value;
        n:=Length(sigma)-1; rhs:=[];
        if n+i+1>native.maxTotalDegree then
            Error("primary comparison exceeds the supplied tensor degrees");
        fi;
        for t in tr.f(sigma) do
            if IsOddInt(t[3]) then
                extend(rhs,tensor.act(tensor.diagonal(i,n,t[1],native.maxTotalDegree),
                    elementIndex(t[2])));
            fi;
        od;
        for faces in cuts[i+1][n+1] do
            extend(rhs,tensorF(sigma{faces[1]},sigma{faces[2]}));
        od;
        if i>0 then
            lower:=comparison(i-1,sigma);
            extend(rhs,lower);
            extend(rhs,List(lower,t->[t[4],t[5],t[6],t[1],t[2],t[3]]));
        fi;
        if n>0 then
            for j in [1..n+1] do
                extend(rhs,comparison(i,sigma{Filtered([1..n+1],q->q<>j)}));
            od;
        fi;
        if failure then return []; fi;
        rhs:=tensor.reduce(rhs);
        value:=account(tensor.boundary(rhs));
        if failure then return []; fi;
        if value<>[] then Error("primary comparison source is not closed"); fi;
        value:=account(tensor.contract(rhs,native.maxTotalDegree));
        if failure then return []; fi;
        lower:=account(tensor.boundary(value));
        if failure then return []; fi;
        if lower<>rhs then Error("primary comparison contraction failed"); fi;
        return value;
    end;
    return function(i,sigma,patterns)
        local value;
        cuts:=patterns; used:=0;
        value:=comparison(i,sigma);
        if failure then
            return rec(status:="unresolved",code:="primary-comparison-term-budget",
                reason:="native primary comparison exceeded its tensor term budget");
        fi;
        # Only here do the group actions disappear: both coefficients are
        # trivial over F_2. Equal native pairs cancel before serialization.
        return rec(status:="computed",terms:=tensor.reduce(List(value,t->
            [1,t[1],t[2]-1,t[4],t[5]-1])));
    end;
end);
