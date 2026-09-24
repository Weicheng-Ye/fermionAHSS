# The final Danus d5, for input degrees zero through three. All local linear
# solves below construct the allowed defining cochains b,c. The universal R
# phase is evaluated by the bundled exact chain kernel, never solved on X.
CallFuncList(function()
    local name, slashes, directory;
    name:=INPUT_FILENAME();
    if not IsEmpty(name) and name[1]<>'/' then name:=Filename(DirectoryCurrent(),name); fi;
    slashes:=Positions(name,'/');
    directory:=Directory(name{[1..Last(slashes)]});
    BindGlobal("KOAHSS_DANUS_KERNEL",Filename(directory,"../python/worker.py"));
end,[]);

BindGlobal("KOAHSS_DanusBatch",function(n,samples)
    local executable, input, output, text, status, answer, serialized;
    if LoadPackage("json")=fail then Error("koAHSS: natural d5 requires the GAP json package"); fi;
    executable:=Filename(DirectoriesSystemPrograms(),"python3");
    if executable=fail then Error("koAHSS: natural d5 requires Python 3.10 or newer"); fi;
    serialized:=CallFuncList(ValueGlobal("GapToJsonString"),[rec(schema:=1,degree:=n,samples:=samples)]);
    input:=InputTextString(serialized); text:=""; output:=OutputTextString(text,true);
    SetPrintFormattingStatus(output,false);
    status:=Process(DirectoryCurrent(),executable,input,output,[KOAHSS_DANUS_KERNEL]);
    CloseStream(input); CloseStream(output);
    if status<>0 then Error("koAHSS: exact Danus d5 kernel failed: ",text); fi;
    answer:=CallFuncList(ValueGlobal("JsonStringToGap"),[text]);
    if not IsRecord(answer) or not IsBound(answer.status) or answer.status<>"computed"
       or not IsBound(answer.audit) or not IsRecord(answer.audit)
       or not IsBound(answer.phases) or Length(answer.phases)<>Length(samples) then
        Error("koAHSS: malformed or incomplete Danus d5 kernel result");
    fi;
    if not ForAll(answer.phases,p->IsList(p) and Length(p)=2 and
        ForAll(p,IsInt) and p[2]>0) then Error("koAHSS: invalid exact Danus phases"); fi;
    return answer;
end);

InstallGlobalFunction(koAHSSNaturalTertiary,function(arg)
    local backend,n,A,options,definingOptions,node,data,ctx,c,transport,degree,
          samples,sampleIndex,chains,j,term,sigma,key,position,fields,name,cochain,
          faces,sample,answer,phase,denominator,numerator,value,result,auditNode,
          canonicalSamples,canonicalPhase,remainder,projectionCarry,sourceCarry,
          closedHomotopyValue;
    if not Length(arg) in [3,4] then Error("koAHSSNaturalTertiary(backend,n,A[,options])"); fi;
    backend:=arg[1]; n:=arg[2]; A:=arg[3]; options:=rec();
    if Length(arg)=4 then options:=arg[4]; fi;
    if not IsInt(n) or n<0 or n>3 then
        Error("koAHSS: natural d5 is implemented only for input degrees 0,1,2,3");
    fi;
    if not IsRecord(options) then Error("koAHSS: natural d5 options must be a record"); fi;
    if ForAny(RecNames(options),name->not name in ["b","c"]) then
        Error("koAHSS: natural d5 accepts only the supplied defining cochains b and c");
    fi;
    KOAHSS_EquationBackend(backend,n);
    KOAHSS_CC_CheckVector(A,backend.dimension(n),"natural d5 input A");
    if ForAny(backend.coboundary(n,A,true),x->x<>0) then
        Error("koAHSS: natural d5 input must be a sign-integral cocycle");
    fi;
    if not IsBound(backend.naturalTransport) and not IsBound(backend.naturalBar) then return rec(status:="unavailable",
        reason:="the fixed natural d5 requires a simplicial transport"); fi;
    definingOptions:=ShallowCopy(options);
    definingOptions.secondaryEvaluator:=koAHSSNaturalSecondary;
    node:=koAHSSDefiningSystem(backend,n,A,definingOptions);
    if node.status="unsatisfied-identity" then
        Error("koAHSS: natural d5 defining-system identity failed: ",node.reason);
    fi;
    if node.status<>"candidate" then return node; fi;
    if not IsBound(node.secondary.modelCochains) then
        Error("koAHSS: natural d5 requires the actual calibrated simplicial defining system");
    fi;
    data:=node.secondary.modelCochains; ctx:=data.context; transport:=ctx.transport;
    closedHomotopyValue:=transport.homotopyValue;
    if IsBound(transport.closedMod2HomotopyValue) then
        closedHomotopyValue:=transport.closedMod2HomotopyValue;
    fi;
    # psi is closed on the envelope modulo two. K and u agree after
    # pairing with it modulo two; the following reduce keeps its literal lift.
    if ctx.isZero(data.psi) then
        c:=ctx.reduce(ctx.lift(n+2,node.c,0));
    else
        c:=ctx.reduce(ctx.add([ctx.lift(n+2,node.c,0),ctx.make(n+2,
            sigma->closedHomotopyValue(n+3,data.psi.value,sigma))]));
    fi;
    degree:=n+4; samples:=[]; sampleIndex:=NewDictionary([],true); chains:=[];
    fields:=rec(s:=data.s,A:=data.A,omega:=data.omega,b:=data.b,c:=c);
    for j in [1..backend.dimension(degree)] do
        chains[j]:=[];
        for term in transport.g(degree,j) do
            sigma:=term[2];
            key:=transport.normalizeSimplex(sigma);
            position:=LookupDictionary(sampleIndex,key);
            if position=fail then
                sample:=rec();
                for name in ["s","A","omega","b","c"] do
                    cochain:=fields.(name);
                    faces:=Combinations([1..degree+1],cochain.degree+1);
                    sample.(name):=List(faces,face->cochain.value(sigma{face}));
                od;
                Add(samples,sample); position:=Length(samples);
                AddDictionary(sampleIndex,ShallowCopy(key),position);
            fi;
            if IsBound(transport.frameSign) then
                value:=transport.frameSign(backend.twists.s,sigma);
            else value:=transport.character(backend.twists.s,sigma[1]); fi;
            Add(chains[j],[position,term[1]*value]);
        od;
    od;
    answer:=KOAHSS_DanusBatch(n,samples);
    phase:=List(chains,chain->Sum(chain,term->term[2]*
        answer.phases[term[1]][1]/answer.phases[term[1]][2]));
    # Audit both carries without changing the literal rational phase.
    # Local canonical residues are framed only during linear projection.
    canonicalSamples:=List(answer.phases,p->(p[1] mod p[2])/p[2]);
    canonicalPhase:=List(chains,chain->Sum(chain,
        term->term[2]*canonicalSamples[term[1]]));
    remainder:=List(phase,p->(NumeratorRat(p) mod DenominatorRat(p))/DenominatorRat(p));
    projectionCarry:=canonicalPhase-remainder;
    sourceCarry:=phase-canonicalPhase;
    if not ForAll(Concatenation(projectionCarry,sourceCarry),IsInt) or
       phase<>remainder+projectionCarry+sourceCarry then
        Error("koAHSS: rational phase carry audit failed");
    fi;
    denominator:=Lcm(Concatenation([1],List(phase,DenominatorRat)));
    numerator:=List(phase,x->NumeratorRat(x*denominator));
    value:=backend.coboundary(degree,numerator,true);
    if ForAny(value,x->x mod denominator<>0) then
        Error("koAHSS: the transported natural d5 phase has a nonintegral boundary");
    fi;
    value:=List(value,x->QuoInt(x,denominator));
    if ForAny(backend.coboundary(n+5,value,true),x->x<>0) then
        Error("koAHSS: natural d5 output is not a sign-integral cocycle");
    fi;
    # Saved evaluations need the exact defining vectors and equation audits,
    # not the lazy simplicial expression graph.  Keeping that
    # graph here retains every per-cochain memo table for each source generator.
    # Copy only the records being pruned: the calculation's defining system
    # and its transport remain intact until this function returns.
    auditNode:=ShallowCopy(node);
    auditNode.secondary:=ShallowCopy(node.secondary);
    Unbind(auditNode.secondary.modelCochains);
    result:=rec(status:="computed",cochain:=value,inputDegree:=n,
        definingSystem:=auditNode,phase:=phase,phaseNumerator:=numerator,modulus:=denominator,
        referenceConvention:="Danus-degree0to3-chi7_tail-epsilon100-eta101-mu0",
        isLocalChoice:=false,isNormalizedOperation:=true,stableDegreeRange:=[0,3],
        oldCorrectionApplied:=false,muR:=0,oddPrimaryCoefficient:=2,
        modelSamples:=Length(samples),kernelAudit:=answer.audit,
        carryAudit:=rec(remainder:=remainder,projectionCarry:=projectionCarry,
            sourceLiftCarry:=sourceCarry,identityVerified:=true));
    result.transportModel:="custom";
    if IsBound(transport.convention) and
         transport.convention="normalized-homogeneous-bar-local-cochains" then
        result.transportModel:="bar";
    fi;
    if result.transportModel="bar" then
        result.inputGroupBarUsed:=true;
        result.barSamples:=Length(samples);
    fi;
    return result;
end);

InstallGlobalFunction(koAHSSNaturalTCallback,function(arg)
    local options;
    if Length(arg)>1 then Error("koAHSSNaturalTCallback([options])"); fi;
    options:=rec(); if Length(arg)=1 then options:=arg[1]; fi;
    if not IsRecord(options) or ForAny(RecNames(options),n->n<>"evaluations")
       or (IsBound(options.evaluations) and not IsList(options.evaluations)) then
        Error("koAHSS: natural T callback accepts an optional evaluations list");
    fi;
    return function(ctx)
        local result;
        result:=koAHSSNaturalTertiary(ctx.backend,ctx.degree,ctx.cochain);
        if not IsBound(ctx.backend.tertiaryEvaluations) then ctx.backend.tertiaryEvaluations:=[]; fi;
        Add(ctx.backend.tertiaryEvaluations,result);
        if IsBound(options.evaluations) then Add(options.evaluations,result); fi;
        if result.status<>"computed" then
            if not result.status in ["unavailable","obstructed"] then
                Error("koAHSS: invalid natural d5 result: ",result.status);
            fi;
            return rec(status:="unresolved",operations:=["T"],
                reasons:=[Concatenation("natural d5: ",result.reason)]);
        fi;
        return result.cochain;
    end;
end);
