# Sparse comparison audit and native-coordinate extension transfer.
# The preflight alone does not certify a nonlinear differential or gauge action.
# Retraction on the integral group ring implies retraction for each coefficient
# character, whereas checking only the trivial/sign characters does not imply
# the integral group-ring identity.
BindGlobal("KOAHSS_ExtensionTransferPreflight",function(arg)
    local backend,k,limits,name,audit,refuse,transport,R,length,unit,n,j,
        degreeAudit,support,chain,term,simplex,image,lifted,entry,key,position,
        coefficients,keys,expected;
    if not Length(arg) in [2,3] then
        Error("KOAHSS_ExtensionTransferPreflight(backend,k[,limits])");
    fi;
    backend:=arg[1]; k:=arg[2];
    if not IsInt(k) or not k in [3..5] then
        Error("extension transfer preflight supports package degrees 3..5");
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
    if not IsRecord(backend) or not IsBound(backend.naturalBar) or
       not IsFunction(backend.naturalBar) then
        return refuse("missing-bar-transport",
            "transfer preflight requires the fixed normalized group-bar comparison");
    fi;
    transport:=backend.naturalBar();
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
    if R!.dimension(0)<>1 then
        audit.degreeZeroRank:=R!.dimension(0);
        return refuse("degree-zero-rank",
            "this transfer preflight requires one degree-zero resolution generator");
    fi;
    length:=ValueGlobal("EvaluateProperty")(R,"length");
    # HAP may leave the last contraction level incomplete even when the
    # corresponding boundary module is present. Keep one spare degree.
    if not IsInt(length) or length<k+3 then
        audit.availableThrough:=length;
        return refuse("resolution-length",
            "preflight needs degrees through k+2 plus one spare HAP contraction degree");
    fi;
    unit:=One(R!.group);
    for n in [0..k+2] do
        degreeAudit:=rec(degree:=n,dimension:=R!.dimension(n),uniqueSupport:=0,
            gTerms:=0,comparisonTerms:=0,basesChecked:=0);
        Add(audit.degrees,degreeAudit); support:=NewDictionary([],true);
        for j in [1..degreeAudit.dimension] do
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
                if LookupDictionary(support,simplex)=fail then
                    if degreeAudit.uniqueSupport>=limits.maxSupport then
                        audit.failure:=rec(degree:=n,basis:=j,stage:="support");
                        return refuse("support-budget","sparse comparison exceeds the per-degree support budget");
                    fi;
                    AddDictionary(support,Immutable(simplex),true);
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
                    position:=Position(R!.elts,entry[2]);
                    if position=fail then
                        Error("bar comparison returned an unindexed resolution group element");
                    fi;
                    key:=[entry[1],position];
                    image:=LookupDictionary(coefficients,key);
                    if image=fail then image:=0; Add(keys,key); fi;
                    AddDictionary(coefficients,key,image+term[1]*entry[3]);
                od;
            od;
            Sort(keys);
            image:=List(Filtered(keys,key->LookupDictionary(coefficients,key)<>0),
                key->[key[1],key[2],LookupDictionary(coefficients,key)]);
            expected:=[[j,Position(R!.elts,unit),1]];
            if image<>expected then
                audit.failure:=rec(degree:=n,basis:=j,image:=image,expected:=expected,
                    coordinateConvention:="[basis,R.elts-index,integer-coefficient]");
                return refuse("chain-retraction-failed",
                    "f composed with g is not the identity on this integral group-ring basis");
            fi;
            degreeAudit.basesChecked:=degreeAudit.basesChecked+1;
        od;
        audit.checkedThrough:=n;
    od;
    audit.status:="checked"; audit.chainRetractionVerified:=true;
    audit.reason:="the bounded group-ring retraction audit passed; nonlinear transfer needs separate identities and certificates";
    return audit;
end);

# Normalize the comparison homotopy on sparse integral group-bar chains.
# With q=1-gf and u=q h q, the operator u boundary u has all SDR side
# conditions. Intermediate chains retain their first vertices and group
# actions; local coefficient frames are applied only at the RPC boundary.
BindGlobal("KOAHSS_ExtensionNormalizedTransport",function(backend,k)
    local audit,tr,R,group,elements,unit,engine,maximum,used,failure,
        reduce,append,act,projectComplement,rawHomotopy,boundary,u,normalizedH,
        cache,order,retained,encode,decode,terms;
    audit:=KOAHSS_ExtensionTransferPreflight(backend,k);
    if audit.status<>"checked" then return audit; fi;
    tr:=backend.naturalBar(); R:=tr.resolution; group:=R!.group;
    if not IsFinite(group) then
        return rec(status:="unresolved",reason:="transferred stacking currently needs a finite group for exact global flags");
    fi;
    unit:=One(group); elements:=AsList(group);
    elements:=Concatenation([unit],Filtered(elements,g->g<>unit));
    engine:=rec(status:="computed",audit:=audit,elements:=elements,
        maxDegree:=k+2,stats:=rec(hRequests:=0,hCacheHits:=0,maxExpansionTerms:=0));
    maximum:=2000000; cache:=NewDictionary("",true); order:=[]; retained:=0;
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
        return List(chain,t->[t[1],List(t[2],v->g*v)]);
    end;
    projectComplement:=function(chain)
        local answer,t,v;
        if chain=fail then return fail; fi;
        answer:=[];
        if not append(answer,chain,1) then return fail; fi;
        for t in chain do
            for v in tr.f(t[2]) do
                if not append(answer,act(tr.g(Length(t[2])-1,v[1]),v[2]),-t[1]*v[3]) then
                    return fail;
                fi;
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
        key:=JoinStringsWithSeparator(List(vertices,v->String(Position(elements,v)-1)),"_");
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
        return List(chain,t->[t[1],t[1]*tr.character(backend.twists.s,t[2][1]),
            List(tr.normalizeSimplex(t[2]),v->Position(elements,v)-1)]);
    end;
    decode:=function(vertices)
        if not IsList(vertices) or IsEmpty(vertices) or
           not ForAll(vertices,x->IsInt(x) and x>=0 and x<Length(elements)) then
            Error("invalid normalized transport vertex indices");
        fi;
        return tr.normalizeSimplex(List(vertices,x->elements[x+1]));
    end;
    terms:=function(kind,vertices)
        local simplex,n,answer;
        simplex:=decode(vertices); n:=Length(simplex)-1;
        if not kind in ["f","h"] or n>engine.maxDegree or
           (kind="h" and n>=engine.maxDegree) then
            Error("normalized transport request exceeds the available degree");
        fi;
        if kind="f" then
            answer:=List(tr.f(simplex),t->[t[1]-1,t[3],t[3]*tr.character(backend.twists.s,t[2])]);
            if Length(answer)>maximum then
                return rec(status:="unresolved",code:="lift-term-budget",
                    reason:="sparse cochain lift exceeded its term budget");
            fi;
        else
            answer:=normalizedH(simplex);
            if answer=fail then return failure; fi;
            answer:=encode(answer);
        fi;
        return rec(status:="computed",terms:=answer);
    end;
    engine.terms:=terms;
    engine.homotopy:=function(simplex)
        local answer;
        answer:=normalizedH(tr.normalizeSimplex(simplex));
        if answer=fail then return fail; fi;
        return act(answer,simplex[1]);
    end;
    engine.g:=List([0..k+2],n->List([1..R!.dimension(n)],j->encode(tr.g(n,j))));
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
        readAnswer,request,n,sign,j,vector,fields,checkState;
    transport:=KOAHSS_ExtensionNormalizedTransport(backend,k);
    if transport.status<>"computed" then return transport; fi;
    if LoadPackage("json")=fail then Error("koFull: transferred stacking requires GAP JSON"); fi;
    elements:=transport.elements;
    model:=rec(status:="computed",modelId:="transferred-normalized-bar",maxDegree:=k,
        certificateLevel:="transfer-R",transportAudit:=transport.audit,
        transportNormalization:=transport.normalization,transportStats:=transport.stats,
        classEquivalenceVerified:=false,sideConditionsByConstruction:=true,
        s:=ShallowCopy(backend.twists.s),omega:=ShallowCopy(backend.twists.omega));
    # BC2 has exactly one normalized simplex orbit in every degree. If g
    # sends each sole native generator to that literal bar generator, the
    # checked inverse f makes this an actual isomorphism, not a proper
    # retract. Its complement and hence the normalized homotopy vanish.
    if Length(elements)=2 and ForAll([0..k+2],n->
       Length(transport.g[n+1])=1 and Length(transport.g[n+1][1])=1 and
       transport.g[n+1][1][1]=[1,1,List([0..n],i->i mod 2)]) then
        model.classEquivalenceVerified:=true;
        model.classEquivalenceReason:="literal normalized C2 bar basis isomorphism through the required degrees";
    fi;
    model.supports:=degree->degree=k and degree in [3..5];
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
        maxDegree:=k+2,elements:=Length(elements),
        multiplication:=List(elements,g->List(elements,h->Position(elements,g*h)-1)),
        ranks:=List([0..k+2],model.dimension),
        ordinary:=List([0..k+1],n->model.matrix(n,false)),
        signed:=List([0..k+1],n->model.matrix(n,true)),
        s:=model.s,omega:=model.omega,g:=transport.g,maxFlagSimplices:=8192);
    executable:=Filename(DirectoriesSystemPrograms(),"python3");
    if executable=fail then Error("koFull: transferred stacking requires Python 3"); fi;
    stream:=fail; cache:=NewDictionary("",true); order:=[];
    model.close:=function()
        if stream<>fail then CloseStream(stream); stream:=fail; fi;
    end;
    readAnswer:=function()
        local line,answer,response;
        while true do
            line:=ReadLine(stream);
            if line=fail then
                model.lastFailure:=rec(status:="unresolved",reason:="transferred stacking worker stopped before returning a result");
                model.close(); Error("koFull: transferred stacking worker stopped");
            fi;
            answer:=JsonStringToGap(line);
            if IsBound(answer.operation) and answer.operation="transport" then
                if not IsBound(answer.vertices) or not IsBound(answer.kind) then
                    model.close(); Error("malformed sparse transport callback");
                fi;
                response:=transport.terms(answer.kind,answer.vertices);
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
    model.d:=function(degree,state)
        return checkState(degree+1,request(rec(operation:="d",degree:=degree,state:=state)).state);
    end;
    model.xtimes:=function(degree,x,y)
        return checkState(degree,request(rec(operation:="xtimes",degree:=degree,state:=x,other:=y)).state);
    end;
    model.act:=function(degree,gauge,canonical)
        return checkState(degree,request(rec(operation:="act",degree:=degree,state:=canonical,gauge:=gauge)).state);
    end;
    model.divideLeft:=function(degree,left,total)
        return checkState(degree,request(rec(operation:="divide_left",degree:=degree,state:=left,other:=total)).state);
    end;
    model.transferCertificate:=function(degree,gauge,canonical,target)
        if model.act(degree,gauge,canonical)<>target then Error("transferred gauge equality failed"); fi;
        return rec(certificateLevel:="transfer-R",nativeEqualityVerified:=true,
            gauge:=StructuralCopy(gauge),canonical:=StructuralCopy(canonical),target:=StructuralCopy(target),
            comparisonAudit:=transport.audit,homotopyNormalization:=transport.normalization,
            searchComplete:=false);
    end;
    model.debugRequest:=request;
    model.transportData:=transport.terms;
    model.barState:=function(degree,state,barModel)
        local samples,ns,j,result;
        samples:=rec(); ns:=[degree-3,degree-2,degree-1,degree+1];
        for j in [1..4] do
            samples.(fields[j]):=List(barModel.simplices(ns[j]),sigma->
                List(sigma,v->Position(elements,v)-1));
        od;
        result:=request(rec(operation:="phi_values",degree:=degree,state:=state,simplices:=samples));
        return result.state;
    end;
    model.certifyPresentation:=function(degree,layers,result)
        if not IsBoundGlobal("KOAHSS_ExtensionTransferCertify") then
            return rec(status:="unresolved",reason:="complete-bar transfer certificate is unavailable");
        fi;
        return CallFuncList(ValueGlobal("KOAHSS_ExtensionTransferCertify"),[backend,model,degree,layers,result]);
    end;
    return model;
end);
