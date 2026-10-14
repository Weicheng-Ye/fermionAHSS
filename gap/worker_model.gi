CallFuncList(function()
    local path,slash;
    path:=INPUT_FILENAME();
    if path[1]<>'/' then path:=Filename(DirectoryCurrent(),path); fi;
    slash:=Last(Positions(path,'/'));
    BindGlobal("KOAHSS_EXTENSION_TRANSFER_WORKER",Concatenation(path{[1..slash]},"../python/stacking_worker.py"));
end,[]);

# Native coordinates and exact native linear algebra; nonlinear values are
# evaluated by the unchanged calibrated formulas through sparse bar callbacks.
BindGlobal("KOAHSS_ExtensionTransferredModel",function(arg)
    local backend,k,mode,transport,elements,model,setup,matrices,executable,stream,cache,order,
        readAnswer,request,n,sign,j,vector,fields,checkState,limited;
    backend:=arg[1]; k:=arg[2]; mode:="extension";
    if Length(arg)>=3 then mode:=arg[3]; fi;
    transport:=KOAHSS_ExtensionNormalizedTransport(backend,k,mode);
    if transport.status<>"computed" then return transport; fi;
    if LoadPackage("json")=fail then Error("koFull: transferred stacking requires GAP JSON"); fi;
    elements:=transport.elements;
    model:=rec(status:="computed",modelId:="transferred-normalized-bar",maxDegree:=k,
        comparison:=transport.comparison,
        certificateLevel:="transfer-R",transportAudit:=transport.audit,
        transportNormalization:=transport.normalization,transportStats:=transport.stats,
        gaugeCompletenessAssumed:=true,sideConditionsByConstruction:=true,
        s:=ShallowCopy(backend.twists.s),omega:=ShallowCopy(backend.twists.omega));
    model.supports:=degree->degree in [1..6] and degree<=k;
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
        s:=model.s,omega:=model.omega,gMode:="lazy",mode:=mode,
        strictIdentities:=transport.strictIdentities,binaryHomotopy:=transport.binaryHomotopy);
    setup.primaryComparison:=transport.primaryComparison;
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
                elif answer.kind="cupR" then
                    response:=rec(status:="computed",terms:=backend.cupMod2(answer.degree,
                        answer.vertices.degrees[1],answer.vertices.vectors[1],
                        answer.vertices.degrees[2],answer.vertices.vectors[2]));
                elif StartsWith(answer.kind,"cupK") then
                    response:=transport.cupHomotopy(Int(answer.kind{[5..Length(answer.kind)]}),
                        answer.vertices,answer.cuts);
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
    # B=C=0 (docs/extension_assembly.md, "Localization at the primes").
    model.primeLocal:=function(prime)
        local local3,withPrime;
        if prime<>3 then Error("koFull: the transferred model localizes at the prime three only"); fi;
        if IsBound(model.localModels) then return model.localModels.three; fi;
        withPrime:=function(requestRecord) requestRecord.prime:=prime; return requestRecord; end;
        local3:=rec(status:="computed",modelId:=Concatenation(model.modelId,"/three-local"),
            maxDegree:=k,localPrime:=prime,certificateLevel:="transfer-R",
            gaugeCompletenessAssumed:=true,layerLimited:=true,
            s:=model.s,omega:=model.omega,supports:=degree->degree in [5,6] and degree<=k,
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
    # A light task (python/light_rows.py, docs/natural_evaluation.md "Light
    # rows"). The worker answers in band: a resource limit is an unresolved
    # record, any other failure an error; the stream stays open.
    model.light:=function(degree,task,data)
        local answer;
        if not IsInt(degree) or degree<1 or degree>k then
            Error("koFull: a light task needs a package degree within the worker capacity");
        fi;
        answer:=request(rec(operation:="light",degree:=degree,task:=task,data:=data));
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
