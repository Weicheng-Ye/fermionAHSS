# Sparse comparison audit for a future transferred extension model.
# This does not implement the nonlinear differential, product, or gauge action.
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
    audit.reason:="the bounded group-ring retraction audit passed; nonlinear transfer is not implemented";
    return audit;
end);
