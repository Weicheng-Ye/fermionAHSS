# Binary comparison of native and simplicial higher diagonals. This keeps
# the fixed f and uses only the tensor contraction of R, never g or bar H.
# Tensor coefficients and cancellations are binary; group actions are kept
# until the tensor has been contracted. See light-transport-reduction (12).
# The same recurrence in flat binary tensors, with the integer terms of
# KOAHSS_BinaryTensorEngine. The comparisons are cached under the budget
# comparisonCache of KOAHSS_EXTENSION_TRANSPORT_LIMITS (entries, retained
# terms, largest stored value); the recursion of a degree-n simplex needs
# the comparisons of its faces and of the lower index, so an evicted value
# recomputes the tree below it. The two exact checks (closedness of the
# source, exactness of the contraction) guard real invariants and cost
# about as much as the contraction; in production they run on the first
# computed simplex of every (i, n) and on one in every primaryAuditEvery
# thereafter (primaryAuditEvery = 1 checks every simplex). A failed check
# is an error, as before. Only the final trivial-character pairing is
# collected, restricted to the requested (p, q) split when one is given.
BindGlobal("KOAHSS_ExtensionPrimaryTransport",function(tr,native,maximum)
    local R,tensor,LB,digits,limits,cache,order,retained,failure,cuts,used,extend,account,
        vertexGroup,labels,nextLabel,label,keyOf,computed,comparison,compute,tensorF,
        indexed,elementIndex,firstIndex;
    if not IsBound(native.binaryTensor) then
        Error("koFull: the shared binary tensor engine of the resolution is required");
    fi;
    R:=tr.resolution; tensor:=native.binaryTensor();
    LB:=2^20; digits:=tensor.digits;
    limits:=ValueGlobal("KOAHSS_EXTENSION_TRANSPORT_LIMITS");
    cache:=NewDictionary(1,true); order:=[]; retained:=0; computed:=[];
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
            if value>=LB then Error("primary comparison: more vertex labels than the key base"); fi;
            AddDictionary(labels,Immutable(v),value);
        fi;
        return value;
    end;
    # The cache key of (i, anchored simplex) as one integer: the number of
    # vertices, then i, then the vertex labels, each below its base, so keys
    # of different dimensions or indices never coincide.
    keyOf:=function(i,anchor)
        local h,v;
        h:=Length(anchor)*4+i;
        for v in anchor do h:=h*LB+label(v); od;
        return h;
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
        local out,a,b,x,y,pl,pr;
        out:=[]; pl:=Length(left)-1; pr:=Length(right)-1;
        x:=Filtered(tr.f(left),a->IsOddInt(a[3]));
        y:=Filtered(tr.f(right),b->IsOddInt(b[3]));
        for a in x do
            for b in y do
                Add(out,tensor.pack(pl,a[1],elementIndex(a[2]),pr,b[1],elementIndex(b[2])));
            od;
        od;
        return tensor.reduce(out);
    end;
    comparison:=function(i,sigma)
        local anchor,key,value,old;
        if i<0 or ForAny([2..Length(sigma)],j->sigma[j]=sigma[j-1]) then return []; fi;
        anchor:=tr.normalizeSimplex(sigma);
        key:=keyOf(i,anchor);
        value:=LookupDictionary(cache,key);
        if value=fail then
            value:=compute(i,anchor);
            if failure then return []; fi;
            MakeImmutable(value);
            if Length(value)<=limits.comparisonCache[3] then
                while not IsEmpty(order) and (Length(order)>=limits.comparisonCache[1]
                      or retained+Length(value)>limits.comparisonCache[2]) do
                    old:=Remove(order,1); retained:=retained-Length(LookupDictionary(cache,old));
                    RemoveDictionary(cache,old);
                od;
                AddDictionary(cache,key,value); Add(order,key); retained:=retained+Length(value);
            fi;
        fi;
        return tensor.act(value,firstIndex(sigma));
    end;
    compute:=function(i,sigma)
        local n,rhs,t,faces,lower,j,value,checked;
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
            extend(rhs,List(lower,tensor.swap));
        fi;
        if n>0 then
            for j in [1..n+1] do
                extend(rhs,comparison(i,sigma{Filtered([1..n+1],q->q<>j)}));
            od;
        fi;
        if failure then return []; fi;
        rhs:=tensor.reduce(rhs);
        if not IsBound(computed[i+1]) then computed[i+1]:=[]; fi;
        if not IsBound(computed[i+1][n+1]) then computed[i+1][n+1]:=0; fi;
        checked:=limits.primaryAuditEvery<=1 or computed[i+1][n+1] mod limits.primaryAuditEvery=0;
        computed[i+1][n+1]:=computed[i+1][n+1]+1;
        if checked then
            value:=account(tensor.boundary(rhs));
            if failure then return []; fi;
            if value<>[] then Error("primary comparison source is not closed"); fi;
        fi;
        value:=account(tensor.contract(rhs,native.maxTotalDegree));
        if failure then return []; fi;
        if checked then
            lower:=account(tensor.boundary(value));
            if failure then return []; fi;
            if lower<>rhs then Error("primary comparison contraction failed"); fi;
        fi;
        return value;
    end;
    return function(i,sigma,patterns,split)
        local value,collapsed,t,u,out;
        cuts:=patterns; used:=0;
        value:=comparison(i,sigma);
        if failure then
            return rec(status:="unresolved",code:="primary-comparison-term-budget",
                reason:="native primary comparison exceeded its tensor term budget");
        fi;
        # Only here do the group actions disappear: both coefficients are
        # trivial over F_2. Equal native pairs cancel before serialization;
        # with a split [p,q] only that part of the tensor is returned.
        collapsed:=[];
        for t in value do
            u:=tensor.unpack(t);
            if split=fail or (u[1]=split[1] and u[4]=split[2]) then
                # [p,j,q,k] in the digits of the engine (no group elements)
                Add(collapsed,((u[1]*digits.J+u[2])*digits.Q+u[4])*digits.J+u[5]);
            fi;
        od;
        out:=[];
        for t in tensor.reduce(collapsed) do
            u:=[QuoInt(t,digits.J*digits.Q*digits.J),QuoInt(t,digits.Q*digits.J) mod digits.J,
                QuoInt(t,digits.J) mod digits.Q,t mod digits.J];
            Add(out,[1,u[1],u[2]-1,u[3],u[4]-1]);
        od;
        return rec(status:="computed",terms:=out);
    end;
end);
