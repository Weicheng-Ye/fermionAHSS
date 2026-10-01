# Sparse comparison audit and native-coordinate extension transfer.
# The preflight alone does not certify a nonlinear differential or gauge action.
# Retraction on the integral group ring implies retraction for each coefficient
# character, whereas checking only the trivial/sign characters does not imply
# the integral group-ring identity.

# Comparison with the simplicial set EX on the free G-set X = G x J, J the
# degree-zero generators of R. Here g(e_j)=cone_b g(boundary e_j) and
# f(sigma)=K(f(boundary sigma)) with the modified contraction
#   K = h + boundary psi - psi boundary,  psi_n(x) = sum_j pi_j(x) h(e_j),
# where pi(x) are the coordinates of x against a right inverse over Z of
# the boundaries of the degree-n generators. Then K(boundary e_j)=e_j and,
# inductively, f g = 1 on every generator. A degree whose boundaries do not
# split off over Z instead gives each generator e_j a private cone vertex
# x_j with its own contraction k_j = h + boundary psi_j - psi_j boundary,
# psi_j(x) = pi_j(x) h(e_j), pi_j(boundary e_j) = 1. Vertices are lists
# [group element, type]; chains and conventions match KOAHSS_NaturalBarTransport.
BindGlobal("KOAHSS_ExtensionCellTransport",function(R,maxDegree)
    local identity,length,m0,typeData,privateOf,split,built,ensure,shiftOf,
        coordinates,liftOf,globalContraction,eltsIndex,synced,
        sync,index,reduceR,reduceBar,actR,actBar,degenerate,anchor,key,hChain,
        boundaryChain,functional,prepare,contraction,fCache,gCache,hCache,f,g,homotopy,
        character,vertexCharacter,characterCache,base,n,j,coefficients,result;
    identity:=One(R!.group);
    length:=ValueGlobal("EvaluateProperty")(R,"length");
    if not IsInt(length) or length<maxDegree then
        return rec(status:="unresolved",code:="resolution-length",
            reason:="the cell comparison needs the resolution through the requested degree");
    fi;
    m0:=R!.dimension(0);
    # First-occurrence index of R!.elts; R's own homotopy may append to it.
    eltsIndex:=NewDictionary(identity,true); synced:=0;
    sync:=function()
        while synced<Length(R!.elts) do
            synced:=synced+1;
            if LookupDictionary(eltsIndex,R!.elts[synced])=fail then
                AddDictionary(eltsIndex,R!.elts[synced],synced);
            fi;
        od;
    end;
    index:=function(element)
        local position;
        sync(); position:=LookupDictionary(eltsIndex,element);
        if position=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(element);
            else Add(R!.elts,element); fi;
            sync(); position:=LookupDictionary(eltsIndex,element);
            if position=fail then Error("cell transport could not index a group element"); fi;
        fi;
        return position;
    end;
    # R terms are [positive basis index, group element, integer coefficient].
    reduceR:=function(terms)
        local sorted,answer,term;
        sorted:=ShallowCopy(terms);
        Sort(sorted,function(a,b) return a{[1,2]}<b{[1,2]}; end);
        answer:=[];
        for term in sorted do
            if not IsEmpty(answer) and Last(answer){[1,2]}=term{[1,2]} then
                Last(answer)[3]:=Last(answer)[3]+term[3];
            else Add(answer,ShallowCopy(term)); fi;
        od;
        return Filtered(answer,term->term[3]<>0);
    end;
    degenerate:=function(simplex)
        local i;
        for i in [2..Length(simplex)] do
            if simplex[i]=simplex[i-1] then return true; fi;
        od;
        return false;
    end;
    # Bar terms are [integer coefficient, vertex list].
    reduceBar:=function(terms)
        local sorted,answer,term;
        sorted:=Filtered(terms,term->term[1]<>0 and not degenerate(term[2]));
        Sort(sorted,function(a,b) return a[2]<b[2]; end);
        answer:=[];
        for term in sorted do
            if not IsEmpty(answer) and Last(answer)[2]=term[2] then
                Last(answer)[1]:=Last(answer)[1]+term[1];
            else Add(answer,[term[1],ShallowCopy(term[2])]); fi;
        od;
        return Filtered(answer,term->term[1]<>0);
    end;
    actR:=function(terms,element,coefficient)
        return List(terms,t->[t[1],element*t[2],coefficient*t[3]]);
    end;
    actBar:=function(terms,element,coefficient)
        return List(terms,t->[coefficient*t[1],List(t[2],v->[element*v[1],v[2]])]);
    end;
    anchor:=function(simplex)
        local inverse;
        if simplex[1][1]=identity then return simplex; fi;
        inverse:=simplex[1][1]^-1;
        return List(simplex,v->[inverse*v[1],v[2]]);
    end;
    key:=simplex->JoinStringsWithSeparator(List(simplex,
        v->Concatenation(String(index(v[1])),":",String(v[2]))),"_");
    # h on an R chain of degree n, as an R chain of degree n+1.
    hChain:=function(n,chain)
        local answer,term,t;
        answer:=[];
        for term in chain do
            for t in R!.homotopy(n,[term[1],index(term[2])]) do
                Add(answer,[AbsInt(t[1]),R!.elts[t[2]],term[3]*SignInt(t[1])]);
            od;
        od;
        return reduceR(answer);
    end;
    boundaryChain:=function(n,chain)
        local answer,term,t;
        answer:=[];
        if n<=0 then return []; fi;
        for term in chain do
            for t in R!.boundary(n,term[1]) do
                Add(answer,[AbsInt(t[1]),term[2]*R!.elts[t[2]],term[3]*SignInt(t[1])]);
            od;
        od;
        return reduceR(answer);
    end;
    # typeData[t]=[degree, generator]: degree zero for the vertex orbits of
    # R0, positive for the private vertices of a non-split degree. A degree is
    # built when first used, after every lower degree, so the type indices,
    # chains and values are those of building all degrees at once.
    typeData:=List([1..m0],j->[0,j]); privateOf:=[]; split:=[]; built:=0;
    ensure:=function(degree)
        local n,chains,keys,M,snf,X,Nplus,i,j;
        if degree>maxDegree then Error("cell transport degree exceeds its capacity"); fi;
        while built<degree do
            n:=built+1;
            chains:=List([1..R!.dimension(n)],j->boundaryChain(n,[[j,identity,1]]));
            keys:=Set(Concatenation(List(chains,c->List(c,t->t{[1,2]}))));
            X:=fail;
            if IsEmpty(chains) then X:=[];
            elif not IsEmpty(keys) then
                M:=List(chains,c->List(keys,k0->Sum(Filtered(c,t->t{[1,2]}=k0),t->t[3])));
                snf:=SmithNormalFormIntegerMatTransforms(M);
                if snf.rank=Length(M) and ForAll([1..snf.rank],i->AbsInt(snf.normal[i][i])=1) then
                    Nplus:=NullMat(Length(keys),Length(M));
                    for i in [1..snf.rank] do Nplus[i][i]:=snf.normal[i][i]; od;
                    X:=snf.coltrans*Nplus*snf.rowtrans;
                    if M*X<>IdentityMat(Length(M)) then Error("cell transport right inverse failed"); fi;
                fi;
            fi;
            if X<>fail then
                split[n]:=rec(keys:=keys,inverse:=X,lift:=[],chains:=chains,shift:=[]);
            else
                privateOf[n]:=[];
                for j in [1..R!.dimension(n)] do
                    Add(typeData,[n,j]); privateOf[n][j]:=Length(typeData);
                od;
            fi;
            built:=n;
        od;
    end;
    # c_j = e_j - h(boundary e_j) = boundary h(e_j).
    shiftOf:=function(n,j)
        if not IsBound(split[n].shift[j]) then
            split[n].shift[j]:=reduceR(Concatenation([[j,identity,1]],
                List(hChain(n-1,split[n].chains[j]),t->[t[1],t[2],-t[3]])));
        fi;
        return split[n].shift[j];
    end;
    coordinates:=function(data,chain)
        local vector,t,i;
        vector:=List(data.keys,k0->0);
        for t in chain do
            i:=PositionSet(data.keys,t{[1,2]});
            if i<>fail then vector[i]:=vector[i]+t[3]; fi;
        od;
        if IsEmpty(vector) then return []; fi;
        return vector*data.inverse;
    end;
    liftOf:=function(n,j)
        if not IsBound(split[n].lift[j]) then
            split[n].lift[j]:=hChain(n,[[j,identity,1]]);
        fi;
        return split[n].lift[j];
    end;
    # K on an R chain of degree n-1, as an R chain of degree n.
    globalContraction:=function(n,chain)
        local answer,v,j;
        ensure(n);
        answer:=hChain(n-1,chain);
        if IsBound(split[n]) then
            v:=coordinates(split[n],chain);
            for j in [1..Length(v)] do
                if v[j]<>0 then
                    Append(answer,List(shiftOf(n,j),t->[t[1],t[2],v[j]*t[3]]));
                fi;
            od;
        fi;
        if n>=2 and IsBound(split[n-1]) then
            v:=coordinates(split[n-1],boundaryChain(n-1,chain));
            for j in [1..Length(v)] do
                if v[j]<>0 then
                    Append(answer,List(liftOf(n-1,j),t->[t[1],t[2],-v[j]*t[3]]));
                fi;
            od;
        fi;
        return reduceR(answer);
    end;
    # pi_j and the chains used by k_j, computed once per private type.
    functional:=[];
    prepare:=function(type)
        local data,m,j,eJ,boundaryJ,coefficients;
        if not IsBound(functional[type]) then
            data:=typeData[type]; m:=data[1]; j:=data[2];
            eJ:=[[j,identity,1]];
            boundaryJ:=boundaryChain(m,eJ);
            coefficients:=List(boundaryJ,t->t[3]);
            if IsEmpty(coefficients) or Gcd(coefficients)<>1 then
                functional[type]:=rec(status:="unresolved",
                    reason:=Concatenation("the boundary of resolution generator ",String(j),
                        " in degree ",String(m)," is not primitive"));
            else
                functional[type]:=rec(status:="computed",
                    keys:=List(boundaryJ,t->t{[1,2]}),
                    weights:=GcdRepresentation(coefficients),
                    # c_j = e_j - h(boundary e_j) = boundary h(e_j).
                    shift:=reduceR(Concatenation(eJ,List(hChain(m-1,boundaryJ),
                        t->[t[1],t[2],-t[3]]))),
                    lift:=fail,eJ:=eJ);
            fi;
        fi;
        return functional[type];
    end;
    # k_type on an R chain of degree n-1, as an R chain of degree n.
    contraction:=function(type,n,chain)
        local answer,m,pi,value,x;
        m:=typeData[type][1];
        if m=0 then return globalContraction(n,chain); fi;
        answer:=hChain(n-1,chain);
        if not n in [m,m+1] then return answer; fi;
        pi:=prepare(type);
        if pi.status<>"computed" then Error(pi.reason); fi;
        value:=function(terms)
            local total,t,i;
            total:=0;
            for t in terms do
                i:=Position(pi.keys,t{[1,2]});
                if i<>fail then total:=total+pi.weights[i]*t[3]; fi;
            od;
            return total;
        end;
        if n=m then
            x:=value(chain);
            if x<>0 then answer:=reduceR(Concatenation(answer,List(pi.shift,t->[t[1],t[2],x*t[3]]))); fi;
        else
            x:=value(boundaryChain(m,chain));
            if x<>0 then
                if pi.lift=fail then pi.lift:=hChain(m,pi.eJ); fi;
                answer:=reduceR(Concatenation(answer,List(pi.lift,t->[t[1],t[2],-x*t[3]])));
            fi;
        fi;
        return answer;
    end;
    fCache:=NewDictionary("",true); gCache:=[]; hCache:=NewDictionary("",true);
    f:=function(simplex)
        local normalized,cacheKey,n,rhs,i,face,answer,type;
        n:=Length(simplex)-1;
        if n>maxDegree then Error("cell transport degree exceeds its capacity"); fi;
        if degenerate(simplex) then return []; fi;
        normalized:=anchor(simplex); cacheKey:=key(normalized);
        answer:=LookupDictionary(fCache,cacheKey);
        if answer=fail then
            type:=normalized[1][2];
            if n=0 then
                if typeData[type][1]=0 then answer:=[[typeData[type][2],identity,1]];
                else answer:=[[1,identity,1]]; fi;
            else
                rhs:=[];
                for i in [1..n+1] do
                    face:=normalized{Filtered([1..n+1],x->x<>i)};
                    Append(rhs,actR(f(face),identity,(-1)^(i-1)));
                od;
                answer:=contraction(type,n,reduceR(rhs));
            fi;
            MakeImmutable(answer); AddDictionary(fCache,cacheKey,answer);
        fi;
        if simplex[1][1]=identity then return answer; fi;
        return actR(answer,simplex[1][1],1);
    end;
    g:=function(n,j)
        local rhs,t,answer,cone;
        if not IsBound(gCache[n+1]) then gCache[n+1]:=[]; fi;
        if not IsBound(gCache[n+1][j]) then
            if n=0 then answer:=[[1,[[identity,j]]]];
            else
                ensure(n);
                rhs:=[];
                for t in R!.boundary(n,j) do
                    Append(rhs,actBar(g(n-1,AbsInt(t[1])),R!.elts[t[2]],SignInt(t[1])));
                od;
                if IsBound(privateOf[n]) then cone:=[identity,privateOf[n][j]];
                else cone:=[identity,1]; fi;
                answer:=reduceBar(List(reduceBar(rhs),t->[t[1],Concatenation([cone],t[2])]));
            fi;
            MakeImmutable(answer); gCache[n+1][j]:=answer;
        fi;
        return gCache[n+1][j];
    end;
    base:=[identity,1];
    # 1-gf = boundary h + h boundary, with the cone at base on anchored simplices.
    homotopy:=function(simplex)
        local normalized,cacheKey,n,rhs,term,i,face,answer;
        n:=Length(simplex)-1;
        if degenerate(simplex) then return []; fi;
        normalized:=anchor(simplex); cacheKey:=key(normalized);
        answer:=LookupDictionary(hCache,cacheKey);
        if answer=fail then
            rhs:=[[1,normalized]];
            for term in f(normalized) do
                Append(rhs,actBar(g(n,term[1]),term[2],-term[3]));
            od;
            if n>0 then
                for i in [1..n+1] do
                    face:=normalized{Filtered([1..n+1],x->x<>i)};
                    Append(rhs,actBar(homotopy(face),identity,(-1)^i));
                od;
            fi;
            answer:=reduceBar(List(reduceBar(rhs),t->[t[1],Concatenation([base],t[2])]));
            MakeImmutable(answer); AddDictionary(hCache,cacheKey,answer);
        fi;
        if simplex[1][1]=identity then return answer; fi;
        return actBar(answer,simplex[1][1],1);
    end;
    # Local sign of a group element and of a vertex: the twist summed along
    # the contraction path from e_1 to g e_1, respectively to f(vertex).
    characterCache:=NewDictionary("",true);
    character:=function(sign,element)
        return vertexCharacter(sign,[element,1]);
    end;
    vertexCharacter:=function(sign,vertex)
        local generator,cacheKey,path,answer;
        if IsInt(sign) and sign=0 then return 1; fi;
        if typeData[vertex[2]][1]=0 then generator:=typeData[vertex[2]][2];
        else generator:=1; fi;
        cacheKey:=Concatenation(JoinStringsWithSeparator(List(sign,String),""),"_",
            String(index(vertex[1])),"_",String(generator));
        answer:=LookupDictionary(characterCache,cacheKey);
        if answer=fail then
            path:=R!.homotopy(0,[generator,index(vertex[1])]);
            answer:=(-1)^(Sum(path,t->SignInt(t[1])*sign[AbsInt(t[1])]) mod 2);
            AddDictionary(characterCache,cacheKey,answer);
        fi;
        return answer;
    end;
    result:=rec(status:="computed",resolution:=R,f:=f,g:=g,homotopy:=homotopy,
        character:=character,vertexCharacter:=vertexCharacter,
        normalizeSimplex:=anchor,actVertex:=function(element,v) return [element*v[1],v[2]]; end,
        vertexGroup:=v->v[1],baseVertex:=base,
        convention:="normalized-homogeneous-bar-local-cochains",comparison:="cells");
    # Refuse up front if some generator needs a private contraction that
    # does not exist; later requests would otherwise fail during a product.
    # The contraction is missing exactly for a boundary that is zero or whose
    # coefficients have gcd other than 1, and such a boundary has no right
    # inverse, so its degree is private: no degree needs to be built here.
    for n in [1..maxDegree] do
        for j in [1..R!.dimension(n)] do
            coefficients:=List(boundaryChain(n,[[j,identity,1]]),t->t[3]);
            if IsEmpty(coefficients) or Gcd(coefficients)<>1 then
                return rec(status:="unresolved",code:="non-primitive-boundary",
                    reason:=Concatenation("the boundary of resolution generator ",String(j),
                        " in degree ",String(n)," is not primitive"));
            fi;
        od;
    od;
    return result;
end);

# The sparse comparison audited one resolution basis element at a time, for
# any comparison record (group bar or cells). verify(n,j) returns the chain
# g(e_j) of the basis element j of R_n (0<=n<=k+2) once its g and f terms fit
# the term budget, its simplices fit the support budget of degree n (counted
# over the union of the chains verified in that degree) and f g(e_j)=e_j.
# Otherwise it returns the refusal record, which is final: every later call
# returns it. The audit counts only the basis elements verified so far.
BindGlobal("KOAHSS_ExtensionTransportVerifier",function(transport,k,limits)
    local R,unit,audit,verified,supports,refusal,refuse,eltsIndex,synced,sync,position,verify;
    R:=transport.resolution; unit:=One(R!.group);
    audit:=rec(status:="unresolved",kind:="transfer-preflight",packageDegree:=k,
        checkedThrough:=-1,requestedThrough:=k+2,transferReady:=false,
        chainRetractionVerified:=false,sideConditionsVerified:=false,
        nonlinearIdentitiesVerified:=false,limits:=limits,
        degrees:=List([0..k+2],n->rec(degree:=n,dimension:=R!.dimension(n),uniqueSupport:=0,
            gTerms:=0,comparisonTerms:=0,basesChecked:=0)),
        gTerms:=0,comparisonTerms:=0,processedTerms:=0,comparison:="group-bar");
    if IsBound(transport.comparison) then audit.comparison:=transport.comparison; fi;
    verified:=List([0..k+2],n->[]); supports:=List([0..k+2],n->NewDictionary([],true));
    refusal:=fail;
    refuse:=function(code,reason)
        audit.status:="unresolved"; audit.code:=code; audit.reason:=reason;
        refusal:=rec(status:="unresolved",code:=code,reason:=reason,failure:=audit.failure);
        return refusal;
    end;
    # The first position in R!.elts, as Position returns it; the contraction
    # of R may append elements.
    eltsIndex:=NewDictionary(unit,true); synced:=0;
    sync:=function()
        while synced<Length(R!.elts) do
            synced:=synced+1;
            if LookupDictionary(eltsIndex,R!.elts[synced])=fail then
                AddDictionary(eltsIndex,R!.elts[synced],synced);
            fi;
        od;
    end;
    position:=function(element) sync(); return LookupDictionary(eltsIndex,element); end;
    verify:=function(n,j)
        local degreeAudit,chain,coefficients,keys,term,simplex,lifted,entry,p,key,image,expected;
        if refusal<>fail then return refusal; fi;
        if not IsInt(n) or n<0 or n>k+2 or not IsInt(j) or j<1 or j>R!.dimension(n) then
            Error("extension transfer verification outside the audited resolution degrees");
        fi;
        if IsBound(verified[n+1][j]) then return transport.g(n,j); fi;
        degreeAudit:=audit.degrees[n+1];
        # The comparison constructs a chain before its size is available.
        # These limits bound accepted/processed support, not construction
        # time or the peak memory of a single pre-existing transport call.
        chain:=transport.g(n,j);
        if audit.processedTerms+Length(chain)>limits.maxTerms then
            audit.failure:=rec(degree:=n,basis:=j,stage:="g");
            return refuse("term-budget","sparse comparison exceeds the term budget");
        fi;
        audit.processedTerms:=audit.processedTerms+Length(chain);
        audit.gTerms:=audit.gTerms+Length(chain);
        degreeAudit.gTerms:=degreeAudit.gTerms+Length(chain);
        coefficients:=NewDictionary([],true); keys:=[];
        for term in chain do
            simplex:=transport.normalizeSimplex(term[2]);
            if LookupDictionary(supports[n+1],simplex)=fail then
                if degreeAudit.uniqueSupport>=limits.maxSupport then
                    audit.failure:=rec(degree:=n,basis:=j,stage:="support");
                    return refuse("support-budget","sparse comparison exceeds the per-degree support budget");
                fi;
                AddDictionary(supports[n+1],Immutable(simplex),true);
                degreeAudit.uniqueSupport:=degreeAudit.uniqueSupport+1;
            fi;
            # f accepts the original homogeneous simplex and therefore
            # retains its group action. Normalization is only for counting.
            lifted:=transport.f(term[2]);
            if audit.processedTerms+Length(lifted)>limits.maxTerms then
                audit.failure:=rec(degree:=n,basis:=j,stage:="f");
                return refuse("term-budget","sparse comparison exceeds the term budget");
            fi;
            audit.processedTerms:=audit.processedTerms+Length(lifted);
            audit.comparisonTerms:=audit.comparisonTerms+Length(lifted);
            degreeAudit.comparisonTerms:=degreeAudit.comparisonTerms+Length(lifted);
            for entry in lifted do
                p:=position(entry[2]);
                if p=fail then
                    Error("bar comparison returned an unindexed resolution group element");
                fi;
                key:=[entry[1],p];
                image:=LookupDictionary(coefficients,key);
                if image=fail then image:=0; Add(keys,key); fi;
                AddDictionary(coefficients,key,image+term[1]*entry[3]);
            od;
        od;
        Sort(keys);
        image:=List(Filtered(keys,key->LookupDictionary(coefficients,key)<>0),
            key->[key[1],key[2],LookupDictionary(coefficients,key)]);
        expected:=[[j,position(unit),1]];
        if image<>expected then
            audit.failure:=rec(degree:=n,basis:=j,image:=image,expected:=expected,
                coordinateConvention:="[basis,R.elts-index,integer-coefficient]");
            return refuse("chain-retraction-failed",
                "f composed with g is not the identity on this integral group-ring basis");
        fi;
        verified[n+1][j]:=true;
        degreeAudit.basesChecked:=degreeAudit.basesChecked+1;
        return chain;
    end;
    return rec(verify:=verify,audit:=audit);
end);

# The complete audit of a comparison record through degree k+2, every basis
# element in turn. The extension model itself verifies a basis element only
# when it first uses its chain (KOAHSS_ExtensionTransportVerifier).
BindGlobal("KOAHSS_ExtensionTransportPreflight",function(arg)
    local k,limits,name,audit,refuse,transport,R,length,verifier,n,j;
    transport:=arg[1]; k:=arg[2];
    if not IsInt(k) or not k in [1..6] then
        Error("extension transfer preflight supports package degrees 1..6");
    fi;
    limits:=rec(maxSupport:=8192,maxTerms:=2000000);
    if Length(arg)=3 then
        if not IsRecord(arg[3]) or
           ForAny(RecNames(arg[3]),name->not name in ["maxSupport","maxTerms"]) then
            Error("extension transfer limits are maxSupport and maxTerms");
        fi;
        for name in RecNames(arg[3]) do limits.(name):=arg[3].(name); od;
    fi;
    if not ForAll(RecNames(limits),name->IsInt(limits.(name)) and limits.(name)>0) then
        Error("extension transfer limits must be positive integers");
    fi;
    audit:=rec(status:="unresolved",kind:="transfer-preflight",packageDegree:=k,
        checkedThrough:=-1,requestedThrough:=k+2,transferReady:=false,
        chainRetractionVerified:=false,sideConditionsVerified:=false,
        nonlinearIdentitiesVerified:=false,limits:=limits,degrees:=[],
        gTerms:=0,comparisonTerms:=0,processedTerms:=0);
    refuse:=function(code,reason)
        audit.code:=code; audit.reason:=reason; return audit;
    end;
    if not IsRecord(transport) or not IsBound(transport.resolution) or
       not IsBound(transport.convention) or
       transport.convention<>"normalized-homogeneous-bar-local-cochains" or
       not IsBound(transport.f) or not IsFunction(transport.f) or
       not IsBound(transport.g) or not IsFunction(transport.g) or
       not IsBound(transport.normalizeSimplex) or not IsFunction(transport.normalizeSimplex) then
        return refuse("missing-bar-transport",
            "transfer preflight requires the fixed normalized group-bar comparison");
    fi;
    R:=transport.resolution;
    audit.comparison:="group-bar";
    if IsBound(transport.comparison) then audit.comparison:=transport.comparison; fi;
    length:=ValueGlobal("EvaluateProperty")(R,"length");
    # HAP may leave the last contraction level incomplete even when the
    # corresponding boundary module is present. Keep one spare degree.
    if not IsInt(length) or length<k+3 then
        audit.availableThrough:=length;
        return refuse("resolution-length",
            "preflight needs degrees through k+2 plus one spare HAP contraction degree");
    fi;
    verifier:=KOAHSS_ExtensionTransportVerifier(transport,k,limits);
    audit:=verifier.audit;
    for n in [0..k+2] do
        for j in [1..R!.dimension(n)] do
            if IsRecord(verifier.verify(n,j)) then return audit; fi;
        od;
        audit.checkedThrough:=n;
    od;
    audit.status:="checked"; audit.chainRetractionVerified:=true;
    audit.reason:="the bounded group-ring retraction audit passed; nonlinear transfer needs separate identities and certificates";
    return audit;
end);

# Whether the fixed group-bar comparison is a retraction through degree
# maxDegree, i.e. f g(e_j)=e_j for every basis element, decided on R alone.
# g(e_j) is the cone at the identity over the cycle g(boundary e_j), and f of
# such a cone over a normalized cycle c is h(f c). Hence f g(e_j)=h(boundary
# e_j) once f g is the identity below degree n, and in degree 0 f g sends
# every generator to the first one. The first failure is the one that the
# audit of all chains (KOAHSS_ExtensionTransportPreflight) reports.
BindGlobal("KOAHSS_ExtensionGroupBarRetraction",function(R,maxDegree)
    local unit,n,j,image,t,c,sorted,answer,term;
    unit:=One(R!.group);
    if R!.dimension(0)<>1 then return rec(retracts:=false,degree:=0,basis:=2); fi;
    for n in [1..maxDegree] do
        for j in [1..R!.dimension(n)] do
            image:=[];
            for t in R!.boundary(n,j) do
                for c in R!.homotopy(n-1,[AbsInt(t[1]),t[2]]) do
                    Add(image,[AbsInt(c[1]),R!.elts[c[2]],SignInt(t[1])*SignInt(c[1])]);
                od;
            od;
            sorted:=ShallowCopy(image);
            Sort(sorted,function(a,b) return a{[1,2]}<b{[1,2]}; end);
            answer:=[];
            for term in sorted do
                if not IsEmpty(answer) and Last(answer){[1,2]}=term{[1,2]} then
                    Last(answer)[3]:=Last(answer)[3]+term[3];
                else Add(answer,ShallowCopy(term)); fi;
            od;
            if Filtered(answer,term->term[3]<>0)<>[[j,unit,1]] then
                return rec(retracts:=false,degree:=n,basis:=j);
            fi;
        od;
    od;
    return rec(retracts:=true,degree:=maxDegree);
end);

BindGlobal("KOAHSS_ExtensionTransferPreflight",function(arg)
    local backend;
    if not Length(arg) in [2,3] then
        Error("KOAHSS_ExtensionTransferPreflight(backend,k[,limits])");
    fi;
    backend:=arg[1];
    if not IsRecord(backend) or not IsBound(backend.naturalBar) or
       not IsFunction(backend.naturalBar) then
        return rec(status:="unresolved",kind:="transfer-preflight",code:="missing-bar-transport",
            reason:="transfer preflight requires the fixed normalized group-bar comparison");
    fi;
    return CallFuncList(KOAHSS_ExtensionTransportPreflight,
        Concatenation([backend.naturalBar()],arg{[2..Length(arg)]}));
end);

# Normalize the comparison homotopy on sparse integral group-bar chains.
# With q=1-gf and u=q h q, the operator u boundary u has all SDR side
# conditions. Intermediate chains retain their first vertices and group
# actions; local coefficient frames are applied only at the RPC boundary.
# Internal transport options. cells:=true forces the cell comparison and
# labels:=true vertex labels without a multiplication table (testing only).
# maxSupport and maxTerms are the resource bounds of the sparse comparison
# and of the normalized homotopy expansion.
BindGlobal("KOAHSS_EXTENSION_TRANSPORT_OVERRIDE",rec(cells:=false,labels:=false,
    maxSupport:=8192,maxTerms:=2000000));

BindGlobal("KOAHSS_ExtensionNormalizedTransport",function(backend,k)
    local bar,override,limits,length,selection,cells,verifier,audit,tr,R,group,elements,unit,
        engine,maximum,used,failure,reduce,append,act,projectComplement,rawHomotopy,boundary,u,
        normalizedH,cache,order,retained,encode,decode,terms,actVertex,
        vertexCharacter,vertexGroup,vertexList,labels,label,tableMode,i;
    if not IsInt(k) or not k in [1..6] then
        Error("extension transfer preflight supports package degrees 1..6");
    fi;
    override:=KOAHSS_EXTENSION_TRANSPORT_OVERRIDE;
    limits:=rec(maxSupport:=override.maxSupport,maxTerms:=override.maxTerms);
    if not IsRecord(backend) or not IsBound(backend.naturalBar) or
       not IsFunction(backend.naturalBar) then
        return rec(status:="unresolved",kind:="transfer-preflight",code:="missing-bar-transport",
            reason:="transfer preflight requires the fixed normalized group-bar comparison");
    fi;
    bar:=backend.naturalBar();
    if not IsRecord(bar) or not IsBound(bar.resolution) or
       not IsBound(bar.convention) or
       bar.convention<>"normalized-homogeneous-bar-local-cochains" or
       not IsBound(bar.f) or not IsFunction(bar.f) or
       not IsBound(bar.g) or not IsFunction(bar.g) or
       not IsBound(bar.normalizeSimplex) or not IsFunction(bar.normalizeSimplex) then
        return rec(status:="unresolved",kind:="transfer-preflight",code:="missing-bar-transport",
            reason:="transfer preflight requires the fixed normalized group-bar comparison");
    fi;
    R:=bar.resolution;
    length:=ValueGlobal("EvaluateProperty")(R,"length");
    # HAP may leave the last contraction level incomplete even when the
    # corresponding boundary module is present. Keep one spare degree.
    if not IsInt(length) or length<k+3 then
        return rec(status:="unresolved",kind:="transfer-preflight",code:="resolution-length",
            reason:="preflight needs degrees through k+2 plus one spare HAP contraction degree",
            availableThrough:=length);
    fi;
    # The group bar retracts onto R only for one degree-zero generator and a
    # normalized contraction, which is decided on R without building any
    # chain. Otherwise compare with the cell complex, on which f g = 1 holds
    # by construction.
    if override.cells then selection:=rec(retracts:=false,method:="override");
    else
        selection:=KOAHSS_ExtensionGroupBarRetraction(R,k+2);
        selection.method:="retraction test on R";
    fi;
    if selection.retracts then tr:=bar;
    else
        cells:=KOAHSS_ExtensionCellTransport(R,k+2);
        if cells.status<>"computed" then return cells; fi;
        tr:=cells;
    fi;
    # A chain g(e_j) is built, checked against the budgets and checked for
    # f g(e_j)=e_j when it is first used, and kept.
    verifier:=KOAHSS_ExtensionTransportVerifier(tr,k,limits);
    audit:=verifier.audit; audit.status:="checking"; audit.mode:="on first use";
    audit.selection:=selection;
    group:=R!.group; unit:=One(group);
    actVertex:=function(element,v) return element*v; end;
    vertexCharacter:=tr.character; vertexGroup:=v->v;
    if IsBound(tr.actVertex) then actVertex:=tr.actVertex; fi;
    if IsBound(tr.vertexCharacter) then vertexCharacter:=tr.vertexCharacter; fi;
    if IsBound(tr.vertexGroup) then vertexGroup:=tr.vertexGroup; fi;
    # Vertex labels. A finite group on the group bar keeps the multiplication
    # table, so labels are group elements in a fixed order; otherwise labels
    # are assigned on demand and GAP normalizes every simplex it receives.
    tableMode:=not IsBound(tr.comparison) and not override.labels and
        (IsPcGroup(group) or IsPermGroup(group) or (HasIsFinite(group) and IsFinite(group)));
    if tableMode then
        elements:=AsList(group);
        vertexList:=Concatenation([unit],Filtered(elements,g->g<>unit));
    else
        elements:=fail;
        if IsBound(tr.baseVertex) then vertexList:=[tr.baseVertex]; else vertexList:=[unit]; fi;
    fi;
    labels:=NewDictionary(vertexList[1],true);
    for i in [1..Length(vertexList)] do AddDictionary(labels,vertexList[i],i-1); od;
    label:=function(v)
        local answer;
        answer:=LookupDictionary(labels,v);
        if answer=fail then
            if tableMode then Error("transport vertex outside the finite group"); fi;
            Add(vertexList,Immutable(v)); answer:=Length(vertexList)-1;
            AddDictionary(labels,vertexList[answer+1],answer);
        fi;
        return answer;
    end;
    engine:=rec(status:="computed",audit:=audit,elements:=vertexList,tableMode:=tableMode,
        comparison:=audit.comparison,vertexLabel:=label,
        maxDegree:=k+2,stats:=rec(hRequests:=0,hCacheHits:=0,maxExpansionTerms:=0,gRequests:=0));
    maximum:=override.maxTerms; cache:=NewDictionary("",true); order:=[]; retained:=0;
    used:=0; failure:=fail;
    append:=function(target,source,coefficient)
        local t;
        if source=fail then return false; fi;
        if used+Length(source)>maximum then
            failure:=rec(status:="unresolved",code:="homotopy-term-budget",
                reason:="normalized sparse homotopy exceeded its expansion term budget");
            return false;
        fi;
        used:=used+Length(source);
        for t in source do Add(target,[coefficient*t[1],t[2]]); od;
        return true;
    end;
    reduce:=function(chain)
        local sorted,answer,t;
        if chain=fail then return fail; fi;
        sorted:=Filtered(chain,t->t[1]<>0 and
            not ForAny([2..Length(t[2])],i->t[2][i]=t[2][i-1]));
        Sort(sorted,function(a,b) return a[2]<b[2]; end); answer:=[];
        for t in sorted do
            if not IsEmpty(answer) and Last(answer)[2]=t[2] then
                Last(answer)[1]:=Last(answer)[1]+t[1];
            else Add(answer,[t[1],t[2]]); fi;
        od;
        return Filtered(answer,t->t[1]<>0);
    end;
    act:=function(chain,g)
        return List(chain,t->[t[1],List(t[2],v->actVertex(g,v))]);
    end;
    projectComplement:=function(chain)
        local answer,t,v,image;
        if chain=fail then return fail; fi;
        answer:=[];
        if not append(answer,chain,1) then return fail; fi;
        for t in chain do
            for v in tr.f(t[2]) do
                image:=verifier.verify(Length(t[2])-1,v[1]);
                if IsRecord(image) then
                    failure:=rec(status:="unresolved",code:=image.code,reason:=image.reason);
                    return fail;
                fi;
                if not append(answer,act(image,v[2]),-t[1]*v[3]) then return fail; fi;
            od;
        od;
        return reduce(answer);
    end;
    rawHomotopy:=function(chain)
        local answer,t;
        if chain=fail then return fail; fi;
        answer:=[];
        for t in chain do
            if not append(answer,tr.homotopy(t[2]),t[1]) then return fail; fi;
        od;
        return reduce(answer);
    end;
    boundary:=function(chain)
        local answer,t,i,vertices;
        if chain=fail then return fail; fi;
        answer:=[];
        for t in chain do
            vertices:=t[2];
            if Length(vertices)>1 then
                for i in [1..Length(vertices)] do
                    if not append(answer,[[t[1]*(-1)^(i-1),
                        vertices{Filtered([1..Length(vertices)],j->j<>i)}]],1) then
                        return fail;
                    fi;
                od;
            fi;
        od;
        return reduce(answer);
    end;
    u:=chain->projectComplement(rawHomotopy(projectComplement(chain)));
    normalizedH:=function(vertices)
        local key,answer,old;
        key:=JoinStringsWithSeparator(List(vertices,v->String(label(v))),"_");
        engine.stats.hRequests:=engine.stats.hRequests+1;
        answer:=LookupDictionary(cache,key);
        if answer<>fail then
            engine.stats.hCacheHits:=engine.stats.hCacheHits+1; return answer;
        fi;
        used:=0; failure:=fail; answer:=u(boundary(u([[1,vertices]])));
        engine.stats.maxExpansionTerms:=Maximum(engine.stats.maxExpansionTerms,used);
        if answer=fail then return fail; fi;
        MakeImmutable(answer);
        if Length(answer)<=200000 then
            while not IsEmpty(order) and (Length(order)>=128 or retained+Length(answer)>200000) do
                old:=Remove(order,1); retained:=retained-Length(LookupDictionary(cache,old));
                RemoveDictionary(cache,old);
            od;
            AddDictionary(cache,key,answer); Add(order,key); retained:=retained+Length(answer);
        fi;
        return answer;
    end;
    encode:=function(chain)
        return List(chain,t->[t[1],t[1]*vertexCharacter(backend.twists.s,t[2][1]),
            List(tr.normalizeSimplex(t[2]),label)]);
    end;
    decode:=function(indices)
        if not IsList(indices) or IsEmpty(indices) or
           not ForAll(indices,x->IsInt(x) and x>=0 and x<Length(vertexList)) then
            Error("invalid normalized transport vertex indices");
        fi;
        return tr.normalizeSimplex(List(indices,x->vertexList[x+1]));
    end;
    # Local values live in the fiber at the first vertex of the anchored
    # simplex; its sign is 1 on the group bar and may differ on cells.
    terms:=function(kind,indices)
        local simplex,n,answer,first;
        simplex:=decode(indices); n:=Length(simplex)-1;
        if not kind in ["f","h"] or n>engine.maxDegree or
           (kind="h" and n>=engine.maxDegree) then
            Error("normalized transport request exceeds the available degree");
        fi;
        first:=vertexCharacter(backend.twists.s,simplex[1]);
        if kind="f" then
            answer:=List(tr.f(simplex),t->[t[1]-1,t[3],
                first*t[3]*tr.character(backend.twists.s,t[2])]);
            if Length(answer)>maximum then
                return rec(status:="unresolved",code:="lift-term-budget",
                    reason:="sparse cochain lift exceeded its term budget");
            fi;
        else
            answer:=normalizedH(simplex);
            if answer=fail then return failure; fi;
            answer:=List(encode(answer),t->[t[1],first*t[2],t[3]]);
        fi;
        return rec(status:="computed",terms:=answer);
    end;
    engine.terms:=terms;
    # The comparison chain of one basis element (basis is 0-based), verified
    # and encoded when the worker first needs it.
    engine.chain:=function(n,basis)
        local chain;
        if not IsInt(n) or n<0 or n>engine.maxDegree or not IsInt(basis)
           or basis<0 or basis>=R!.dimension(n) then
            Error("comparison chain request outside the resolution basis");
        fi;
        engine.stats.gRequests:=engine.stats.gRequests+1;
        chain:=verifier.verify(n,basis+1);
        if IsRecord(chain) then
            return rec(status:="unresolved",code:=chain.code,reason:=chain.reason);
        fi;
        return rec(status:="computed",terms:=encode(chain));
    end;
    engine.homotopy:=function(simplex)
        local answer;
        answer:=normalizedH(tr.normalizeSimplex(simplex));
        if answer=fail then return fail; fi;
        return act(answer,vertexGroup(simplex[1]));
    end;
    engine.normalization:="q h q boundary q h q, q=1-gf";
    return engine;
end);

CallFuncList(function()
    local path,slash;
    path:=INPUT_FILENAME();
    if path[1]<>'/' then path:=Filename(DirectoryCurrent(),path); fi;
    slash:=Last(Positions(path,'/'));
    BindGlobal("KOAHSS_EXTENSION_TRANSFER_WORKER",Concatenation(path{[1..slash]},"../python/extension_transfer.py"));
end,[]);

# Native coordinates and exact native linear algebra; nonlinear values are
# evaluated by the unchanged calibrated formulas through sparse bar callbacks.
BindGlobal("KOAHSS_ExtensionTransferredModel",function(backend,k)
    local transport,elements,model,setup,matrices,executable,stream,cache,order,
        readAnswer,request,n,sign,j,vector,fields,checkState,limited;
    transport:=KOAHSS_ExtensionNormalizedTransport(backend,k);
    if transport.status<>"computed" then return transport; fi;
    if LoadPackage("json")=fail then Error("koFull: transferred stacking requires GAP JSON"); fi;
    elements:=transport.elements;
    model:=rec(status:="computed",modelId:="transferred-normalized-bar",maxDegree:=k,
        comparison:=transport.comparison,
        certificateLevel:="transfer-R",transportAudit:=transport.audit,
        transportNormalization:=transport.normalization,transportStats:=transport.stats,
        gaugeCompletenessAssumed:=true,sideConditionsByConstruction:=true,
        s:=ShallowCopy(backend.twists.s),omega:=ShallowCopy(backend.twists.omega));
    model.supports:=degree->degree=k and degree in [1..6];
    # The stacking operations accept a last layer index (see below).
    model.layerLimited:=true;
    model.dimension:=function(n)
        if n<0 then return 0; fi;
        if n>k+2 then Error("transferred cochain degree exceeds model capacity"); fi;
        return backend.dimension(n);
    end;
    model.coboundary:=function(n,vector,signed)
        KOAHSS_CC_CheckVector(vector,model.dimension(n),"native extension cochain");
        if n<0 then return List([1..model.dimension(n+1)],i->0); fi;
        return backend.coboundary(n,vector,signed);
    end;
    matrices:=rec();
    model.matrix:=function(n,signed)
        local key,matrix,j,vector;
        key:=Concatenation(String(n),"_",String(signed));
        if not IsBound(matrices.(key)) then
            matrix:=[];
            for j in [1..model.dimension(n)] do
                vector:=List([1..model.dimension(n)],i->0); vector[j]:=1;
                Add(matrix,model.coboundary(n,vector,signed));
            od;
            MakeImmutable(matrix); matrices.(key):=matrix;
        fi;
        return matrices.(key);
    end;
    model.lift:=function(n,vector,signed)
        KOAHSS_CC_CheckVector(vector,model.dimension(n),"native extension lift");
        if signed then return ShallowCopy(vector); fi;
        return List(vector,x->x mod 2);
    end;
    model.project:=function(n,vector,signed)
        KOAHSS_CC_CheckVector(vector,model.dimension(n),"native extension projection");
        return ShallowCopy(vector);
    end;
    model.zero:=degree->rec(A:=List([1..model.dimension(degree-3)],i->0),
        B:=List([1..model.dimension(degree-2)],i->0),C:=List([1..model.dimension(degree-1)],i->0),
        D:=List([1..model.dimension(degree+1)],i->0));
    setup:=rec(schema:=1,operation:="setup",model:="transferred-normalized-bar",k:=k,
        maxDegree:=k+2,ranks:=List([0..k+2],model.dimension),
        ordinary:=List([0..k+1],n->model.matrix(n,false)),
        signed:=List([0..k+1],n->model.matrix(n,true)),
        s:=model.s,omega:=model.omega,gMode:="lazy");
    if transport.tableMode then
        setup.vertexMode:="table"; setup.elements:=Length(elements);
        setup.multiplication:=List(elements,g->List(elements,h->transport.vertexLabel(g*h)));
    else
        setup.vertexMode:="labels";
    fi;
    executable:=Filename(DirectoriesSystemPrograms(),"python3");
    if executable=fail then Error("koFull: transferred stacking requires Python 3"); fi;
    stream:=fail; cache:=NewDictionary("",true); order:=[];
    model.close:=function()
        if stream<>fail then CloseStream(stream); stream:=fail; fi;
    end;
    readAnswer:=function()
        local line,answer,response;
        while true do
            line:=KOAHSS_ReadWorkerLine(stream);
            if line=fail then
                model.lastFailure:=rec(status:="unresolved",reason:="transferred stacking worker stopped before returning a result");
                model.close(); Error("koFull: transferred stacking worker stopped");
            fi;
            answer:=JsonStringToGap(line);
            if IsBound(answer.operation) and answer.operation="transport" then
                if not IsBound(answer.kind) then
                    model.close(); Error("malformed sparse transport callback");
                fi;
                if answer.kind="g" then
                    if not IsBound(answer.degree) or not IsBound(answer.basis) then
                        model.close(); Error("malformed sparse transport callback");
                    fi;
                    response:=transport.chain(answer.degree,answer.basis);
                else
                    if not IsBound(answer.vertices) then
                        model.close(); Error("malformed sparse transport callback");
                    fi;
                    response:=transport.terms(answer.kind,answer.vertices);
                fi;
                WriteLine(stream,GapToJsonString(response));
            else return answer; fi;
        od;
    end;
    request:=function(input)
        local encoded,found,answer,old;
        encoded:=GapToJsonString(input); found:=LookupDictionary(cache,encoded);
        if found<>fail then return found; fi;
        if stream=fail then
            stream:=InputOutputLocalProcess(DirectoryCurrent(),executable,["-u",KOAHSS_EXTENSION_TRANSFER_WORKER]);
            WriteLine(stream,GapToJsonString(setup)); answer:=readAnswer();
            if not IsBound(answer.status) or answer.status<>"computed" then
                model.lastFailure:=answer; model.close(); Error("koFull: transferred worker setup failed: ",answer);
            fi;
        fi;
        WriteLine(stream,encoded); answer:=readAnswer();
        if not IsBound(answer.status) or answer.status<>"computed" then
            model.lastFailure:=answer; model.close(); Error("koFull: transferred stacking operation failed: ",answer);
        fi;
        MakeImmutable(answer);
        if Length(order)>=128 then old:=Remove(order,1); RemoveDictionary(cache,old); fi;
        Add(order,encoded); AddDictionary(cache,encoded,answer);
        return answer;
    end;
    fields:=["A","B","C","D"];
    checkState:=function(degree,state)
        local ns,j;
        ns:=[degree-3,degree-2,degree-1,degree+1];
        if not IsRecord(state) then Error("transferred result must be a state"); fi;
        for j in [1..4] do
            if not IsBound(state.(fields[j])) then Error("transferred result is missing a state layer"); fi;
            KOAHSS_CC_CheckVector(state.(fields[j]),model.dimension(ns[j]),"transferred result");
            if j in [2,3] and not ForAll(state.(fields[j]),x->x in [0,1]) then
                Error("transferred binary result is not reduced");
            fi;
        od;
        return state;
    end;
    # Every operation takes an optional last layer index ``upto`` (A=0 to
    # D=3, default 3): the worker computes the layers through it and returns
    # zeros above it, so a relation measured only down to its target layer
    # never evaluates the D-layer stacking correction.
    limited:=function(arg,position,requestRecord)
        local upto;
        upto:=3;
        if Length(arg)>=position then
            upto:=arg[position];
            if not IsInt(upto) or upto<0 or upto>3 then Error("koFull: upto must be a layer index from zero through three"); fi;
        fi;
        if upto<3 then requestRecord.upto:=upto; fi;
        return requestRecord;
    end;
    model.d:=function(arg)
        local degree,state;
        degree:=arg[1]; state:=arg[2];
        return checkState(degree+1,request(limited(arg,3,rec(operation:="d",degree:=degree,state:=state))).state);
    end;
    model.xtimes:=function(arg)
        local degree,x,y;
        degree:=arg[1]; x:=arg[2]; y:=arg[3];
        return checkState(degree,request(limited(arg,4,rec(operation:="xtimes",degree:=degree,state:=x,other:=y))).state);
    end;
    model.act:=function(arg)
        local degree,gauge,canonical;
        degree:=arg[1]; gauge:=arg[2]; canonical:=arg[3];
        return checkState(degree,request(limited(arg,4,rec(operation:="act",degree:=degree,state:=canonical,gauge:=gauge))).state);
    end;
    model.divideLeft:=function(arg)
        local degree,left,total;
        degree:=arg[1]; left:=arg[2]; total:=arg[3];
        return checkState(degree,request(limited(arg,4,rec(operation:="divide_left",degree:=degree,state:=left,other:=total))).state);
    end;
    model.transferCertificate:=function(arg)
        local degree,gauge,canonical,target,upto,certificate;
        degree:=arg[1]; gauge:=arg[2]; canonical:=arg[3]; target:=arg[4];
        upto:=3; if Length(arg)>=5 then upto:=arg[5]; fi;
        if model.act(degree,gauge,canonical,upto)<>target then Error("transferred gauge equality failed"); fi;
        certificate:=rec(certificateLevel:="transfer-R",nativeEqualityVerified:=true,
            gaugeCompletenessAssumed:=true,
            gauge:=StructuralCopy(gauge),canonical:=StructuralCopy(canonical),target:=StructuralCopy(target),
            comparisonAudit:=StructuralCopy(transport.audit),homotopyNormalization:=transport.normalization,
            searchComplete:=false);
        if upto<3 then certificate.comparedLayers:=fields{[1..upto+1]}; fi;
        return certificate;
    end;
    # The two-layer three-local model on the same worker process: its
    # requests carry the prime, the worker answers them with
    # extension_three_local, and its states are the four-layer records with
    # B=C=0 (doc/extensions.md, "Localization at the primes").
    model.primeLocal:=function(prime)
        local local3,withPrime;
        if prime<>3 then Error("koFull: the transferred model localizes at the prime three only"); fi;
        if IsBound(model.localModels) then return model.localModels.three; fi;
        withPrime:=function(requestRecord) requestRecord.prime:=prime; return requestRecord; end;
        local3:=rec(status:="computed",modelId:=Concatenation(model.modelId,"/three-local"),
            maxDegree:=k,localPrime:=prime,certificateLevel:="transfer-R",
            gaugeCompletenessAssumed:=true,layerLimited:=true,
            s:=model.s,omega:=model.omega,supports:=degree->degree=k and k in [5,6],
            dimension:=model.dimension,coboundary:=model.coboundary,matrix:=model.matrix,
            lift:=model.lift,project:=model.project,zero:=model.zero);
        local3.d:=function(arg)
            local degree,state;
            degree:=arg[1]; state:=arg[2];
            return checkState(degree+1,request(withPrime(limited(arg,3,
                rec(operation:="d",degree:=degree,state:=state)))).state);
        end;
        local3.xtimes:=function(arg)
            local degree,x,y;
            degree:=arg[1]; x:=arg[2]; y:=arg[3];
            return checkState(degree,request(withPrime(limited(arg,4,
                rec(operation:="xtimes",degree:=degree,state:=x,other:=y)))).state);
        end;
        local3.act:=function(arg)
            local degree,gauge,canonical;
            degree:=arg[1]; gauge:=arg[2]; canonical:=arg[3];
            return checkState(degree,request(withPrime(limited(arg,4,
                rec(operation:="act",degree:=degree,state:=canonical,gauge:=gauge)))).state);
        end;
        local3.divideLeft:=function(arg)
            local degree,left,total;
            degree:=arg[1]; left:=arg[2]; total:=arg[3];
            return checkState(degree,request(withPrime(limited(arg,4,
                rec(operation:="divide_left",degree:=degree,state:=left,other:=total)))).state);
        end;
        local3.transferCertificate:=function(arg)
            local degree,gauge,canonical,target,upto,certificate;
            degree:=arg[1]; gauge:=arg[2]; canonical:=arg[3]; target:=arg[4];
            upto:=3; if Length(arg)>=5 then upto:=arg[5]; fi;
            if local3.act(degree,gauge,canonical,upto)<>target then Error("three-local gauge equality failed"); fi;
            certificate:=rec(certificateLevel:="transfer-R",nativeEqualityVerified:=true,
                gaugeCompletenessAssumed:=true,localPrime:=prime,
                gauge:=StructuralCopy(gauge),canonical:=StructuralCopy(canonical),target:=StructuralCopy(target),
                comparisonAudit:=StructuralCopy(transport.audit),homotopyNormalization:=transport.normalization,
                searchComplete:=false);
            if upto<3 then certificate.comparedLayers:=fields{[1..upto+1]}; fi;
            return certificate;
        end;
        model.localModels:=rec(three:=local3);
        return local3;
    end;
    # A light task (python/extension_light.py, doc/extensions.md "Light
    # rows"). The worker answers in band: a resource limit is an unresolved
    # record, any other failure an error; the stream stays open.
    model.light:=function(task,data)
        local answer;
        answer:=request(rec(operation:="light",degree:=k,task:=task,data:=data));
        if IsBound(answer.refused) then
            return rec(status:="unresolved",code:="resource-limit",reason:=answer.refused.reason);
        fi;
        if IsBound(answer.error) then
            return rec(status:="error",reason:=Concatenation(answer.error.exception,": ",
                answer.error.reason),traceback:=answer.error.traceback);
        fi;
        return rec(status:="computed",result:=answer.result);
    end;
    model.debugRequest:=request;
    model.transportData:=transport.terms;
    # Explicit export for developer reference comparisons only. koFull never
    # calls this helper or constructs a complete-bar model for certification.
    model.barState:=function(degree,state,barModel)
        local samples,ns,j,result;
        samples:=rec(); ns:=[degree-3,degree-2,degree-1,degree+1];
        for j in [1..4] do
            samples.(fields[j]):=List(barModel.simplices(ns[j]),sigma->
                List(sigma,transport.vertexLabel));
        od;
        result:=request(rec(operation:="phi_values",degree:=degree,state:=state,simplices:=samples));
        return result.state;
    end;
    return model;
end);
