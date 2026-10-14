# Free Z[W x C2]-resolution R (x) S of a crystallographic group W times a
# cyclic group C2 = <t>, with diagonal action and an integral contracting
# homotopy. HAP's ResolutionDirectProduct assumes one degree-zero generator,
# which fails for SpaceGroupCohomology resolutions of most wallpaper groups.
#
# S is the periodic resolution: S_q = Z[C2] f_q, d f_q = (t-(-1)^q) f_{q-1},
# with contraction h(f_q)=0, h(t f_q)=f_{q+1}. On R (x) S,
#   d(r(x)s) = dr(x)s + (-1)^|r| r(x)ds,
#   h(r(x)s) = h_R(r)(x)s + [|r|=0] e_1(x)h_S(s),
# so dh+hd=1 in positive degrees and dh = 1-[1,identity] in degree zero,
# the HAP convention used by the page engine's bar transport.
# Group elements are the 4x4 matrices blockdiag(A,eps): A in R!.group
# (an affine matrix group acting on the left) and eps=-1 for t.

BindGlobal("KOAHSSBatchProductWithC2",function(R,length)
    local W,n0,embed,T,gens,E,elts,eltPos,synced,sync,index,indexR,dimR,dim,offset,
        decode,boundary,homotopy,rlen,dimCache,boundaryCache,result;
    W:=R!.group; n0:=Length(One(W));
    rlen:=EvaluateProperty(R,"length");
    if rlen<length then Error("product resolution: factor resolution is too short"); fi;
    embed:=function(A,eps)
        local M,i,j;
        M:=NullMat(n0+1,n0+1);
        for i in [1..n0] do for j in [1..n0] do M[i][j]:=A[i][j]; od; od;
        M[n0+1][n0+1]:=eps;
        return M;
    end;
    T:=embed(One(W),-1);
    gens:=Concatenation(List(GeneratorsOfGroup(W),g->embed(g,1)),[T]);
    E:=Group(gens,embed(One(W),1));
    elts:=[Immutable(embed(One(W),1))]; eltPos:=NewDictionary(elts[1],true);
    AddDictionary(eltPos,elts[1],1); synced:=1;
    # HAP helpers such as ResolutionBoundaryOfWord append to elts directly;
    # index the first occurrence of every element, as Position would.
    sync:=function()
        while synced<Length(elts) do
            synced:=synced+1;
            if LookupDictionary(eltPos,elts[synced])=fail then
                AddDictionary(eltPos,Immutable(elts[synced]),synced);
            fi;
        od;
    end;
    index:=function(x)
        local i;
        sync(); i:=LookupDictionary(eltPos,x);
        if i=fail then
            Add(elts,Immutable(x)); i:=Length(elts); synced:=i;
            AddDictionary(eltPos,elts[i],i);
        fi;
        return i;
    end;
    indexR:=function(A)
        local i;
        i:=Position(R!.elts,A);
        if i=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(A); else Add(R!.elts,A); fi;
            i:=Position(R!.elts,A);
        fi;
        return i;
    end;
    dimR:=function(p)
        if p<0 or p>rlen then return 0; fi;
        return R!.dimension(p);
    end;
    # Degree-n generators are blocks q=0..n, block q holding e_i (x) f_q
    # with |e_i|=n-q. Degree 0 generator 1 is e_1 (x) f_0.
    offset:=function(n,q) return Sum([0..q-1],j->dimR(n-j)); end;
    dimCache:=[];
    dim:=function(n)
        if n<0 then return 0; fi;
        if n>length then Error("product resolution: degree beyond its length"); fi;
        if not IsBound(dimCache[n+1]) then dimCache[n+1]:=offset(n,n+1); fi;
        return dimCache[n+1];
    end;
    decode:=function(n,a)
        local q,i;
        i:=a;
        for q in [0..n] do
            if i<=dimR(n-q) then return [n-q,q,i]; fi;
            i:=i-dimR(n-q);
        od;
        Error("product resolution: generator index out of range");
    end;
    boundaryCache:=[];
    boundary:=function(n,a)
        local sgn,v,p,q,i,word,t,id;
        sgn:=SignInt(a); a:=AbsInt(a);
        if n<=0 then return []; fi;
        if not IsBound(boundaryCache[n]) then boundaryCache[n]:=[]; fi;
        if not IsBound(boundaryCache[n][a]) then
            v:=decode(n,a); p:=v[1]; q:=v[2]; i:=v[3]; word:=[];
            if p>=1 then
                for t in R!.boundary(p,i) do
                    Add(word,[SignInt(t[1])*(offset(n-1,q)+AbsInt(t[1])),
                        index(embed(R!.elts[t[2]],1))]);
                od;
            fi;
            if q>=1 then
                id:=offset(n-1,q-1)+i;
                Add(word,[(-1)^p*id,index(T)]);
                if IsOddInt(q) then Add(word,[-(-1)^p*id,1]);
                else Add(word,[(-1)^p*id,1]); fi;
            fi;
            boundaryCache[n][a]:=Immutable(word);
        fi;
        if sgn>0 then return StructuralCopy(boundaryCache[n][a]); fi;
        return List(boundaryCache[n][a],t->[-t[1],t[2]]);
    end;
    homotopy:=function(n,term)
        local a,sgn,x,A,eps,v,p,q,i,word,t;
        sgn:=SignInt(term[1]); a:=AbsInt(term[1]);
        x:=elts[term[2]]; eps:=x[n0+1][n0+1];
        A:=List(x{[1..n0]},row->row{[1..n0]});
        v:=decode(n,a); p:=v[1]; q:=v[2]; i:=v[3]; word:=[];
        if p<=rlen-2 then
            for t in R!.homotopy(p,[i,indexR(A)]) do
                Add(word,[sgn*SignInt(t[1])*(offset(n+1,q)+AbsInt(t[1])),
                    index(embed(R!.elts[t[2]],eps))]);
            od;
        fi;
        if p=0 and eps=-1 then
            Add(word,[sgn*(offset(n+1,q+1)+1),1]);
        fi;
        return word;
    end;
    result:=Objectify(HapResolution,rec(dimension:=dim,boundary:=boundary,
        homotopy:=homotopy,elts:=elts,group:=E,appendToElts:=function(x) index(x); end,
        properties:=[["length",length],["reduced",true],["type","resolution"],
            ["characteristic",0]]));
    return rec(resolution:=result,embed:=embed,t:=T,factorDegree:=n0);
end);

# Exact identity checks: d d = 0 on generators and the contraction identity
# on the supplied letters (generator, element index). Chains are collected
# directly as (generator, element index) -> coefficient.
BindGlobal("KOAHSSBatchCheckResolution",function(R,length,elementIndices)
    local index,collect,boundaryOf,homotopyOf,n,a,g,chain,count;
    index:=function(x)
        local i;
        i:=Position(R!.elts,x);
        if i=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(x); else Add(R!.elts,x); fi;
            i:=Position(R!.elts,x);
        fi;
        return i;
    end;
    collect:=function(terms)
        local keys,values,t,i,answer;
        keys:=[]; values:=[];
        for t in terms do
            i:=Position(keys,t{[1,2]});
            if i=fail then Add(keys,t{[1,2]}); Add(values,t[3]);
            else values[i]:=values[i]+t[3]; fi;
        od;
        answer:=[];
        for i in [1..Length(keys)] do
            if values[i]<>0 then Add(answer,[keys[i][1],keys[i][2],values[i]]); fi;
        od;
        return answer;
    end;
    # A chain is a list [generator, element index, coefficient].
    boundaryOf:=function(n,chain)
        local answer,t,u;
        answer:=[];
        if n<=0 then return []; fi;
        for t in chain do
            for u in R!.boundary(n,t[1]) do
                Add(answer,[AbsInt(u[1]),index(R!.elts[t[2]]*R!.elts[u[2]]),t[3]*SignInt(u[1])]);
            od;
        od;
        return collect(answer);
    end;
    homotopyOf:=function(n,chain)
        local answer,t,u;
        answer:=[];
        for t in chain do
            for u in R!.homotopy(n,[t[1],t[2]]) do
                Add(answer,[AbsInt(u[1]),u[2],t[3]*SignInt(u[1])]);
            od;
        od;
        return collect(answer);
    end;
    count:=0;
    for n in [2..length] do
        for a in [1..R!.dimension(n)] do
            if boundaryOf(n-1,boundaryOf(n,[[a,1,1]]))<>[] then
                Error("resolution check: d d is not zero in degree ",n);
            fi;
            count:=count+1;
        od;
    od;
    for n in [0..length-1] do
        for a in [1..R!.dimension(n)] do
            for g in elementIndices do
                chain:=[[a,g,1]];
                if n=0 then
                    chain:=collect(Concatenation(boundaryOf(1,homotopyOf(0,chain)),
                        [[a,g,-1],[1,index(One(R!.group)),1]]));
                else
                    chain:=collect(Concatenation(boundaryOf(n+1,homotopyOf(n,chain)),
                        homotopyOf(n-1,boundaryOf(n,chain)),[[a,g,-1]]));
                fi;
                if chain<>[] then
                    Error("resolution check: contraction identity fails in degree ",n);
                fi;
                count:=count+1;
            od;
        od;
    od;
    return count;
end);
