# Abelian extensions in a fixed, marked presentation.  Rows are relations;
# c -> c V and rows of V^-1 express the Smith basis in the named generators.
BindGlobal("KOAHSS_ExtensionPresentationSerial",rec(next:=0));
BindGlobal("KOAHSS_ExtensionSmith",function(matrix,ids)
    local n, m, snf, orders, active, inverse, group, basis;
    n:=Length(ids); m:=Length(matrix);
    if not ForAll(matrix,r->IsList(r) and Length(r)=n and ForAll(r,IsInt)) then
        Error("koAHSS: malformed extension relation matrix");
    fi;
    if n=0 or m=0 then
        snf:=rec(normal:=StructuralCopy(matrix),rowtrans:=IdentityMat(m),
            coltrans:=IdentityMat(n),rank:=0);
    else
        snf:=SmithNormalFormIntegerMatTransforms(matrix);
        if snf.rowtrans*matrix*snf.coltrans<>snf.normal then
            Error("koAHSS: extension Smith identity failed");
        fi;
    fi;
    orders:=List([1..n],i->0);
    for n in [1..snf.rank] do orders[n]:=AbsInt(snf.normal[n][n]); od;
    active:=Filtered([1..Length(ids)],i->orders[i]<>1);
    if IsEmpty(ids) then inverse:=[]; else inverse:=Inverse(snf.coltrans); fi;
    group:=AbelianGroup(IsPcpGroup,orders{active});
    basis:=rec(orders:=orders{active},expressions:=inverse{active},
        generatorIds:=ShallowCopy(ids));
    return rec(group:=group,invariants:=AbelianInvariants(group),basis:=basis,
        U:=snf.rowtrans,V:=snf.coltrans,inverseV:=inverse,S:=snf.normal,
        rank:=snf.rank,activeIndices:=active,diagonal:=orders);
end);

BindGlobal("KOAHSS_ExtensionGroupCoordinates",function(group,generators,x)
    local independent, coordinates, i, result;
    independent:=IndependentGeneratorsOfAbelianGroup(group);
    if generators<>independent then
        Error("koAHSS: extension coordinates require the recorded independent basis");
    fi;
    if not x in group then Error("koAHSS: extension element is outside its group"); fi;
    if IsEmpty(generators) then return []; fi;
    coordinates:=IndependentGeneratorExponents(group,x);
    result:=One(group);
    for i in [1..Length(generators)] do result:=result*generators[i]^coordinates[i]; od;
    if result<>x then Error("koAHSS: invalid extension group coordinates"); fi;
    return coordinates;
end);

InstallGlobalFunction(koAHSSExtensionFromLayers,function(arg)
    local input, oracle, layers, name, layer, ids, matrix, vectors, stages,
        model, smith, oldCount, count, i, j, order, response, row, newIds,
        previous, inclusion, projection, stage, options, reason, presentationPrefix;
    if not Length(arg) in [2,3] or not IsRecord(arg[1]) then
        Error("usage: koAHSSExtensionFromLayers(layers, oracle[, options])");
    fi;
    input:=arg[1]; oracle:=arg[2]; options:=rec();
    if Length(arg)=3 then options:=arg[3]; fi;
    if not IsRecord(options) or not IsEmpty(RecNames(options)) then
        Error("koAHSS: extension options must currently be an empty record");
    fi;
    if oracle<>fail and not IsFunction(oracle) then
        Error("koAHSS: extension oracle must be a function or fail");
    fi;
    if not ForAll(RecNames(input),n->n in ["A","B","C","D"]) then
        Error("koAHSS: extension layers must be named A, B, C, D");
    fi;
    layers:=[];
    for name in ["D","C","B","A"] do
        if IsBound(input.(name)) then
            if IsList(input.(name)) then layer:=rec(orders:=ShallowCopy(input.(name)));
            elif IsRecord(input.(name)) then layer:=ShallowCopy(input.(name));
            else Error("koAHSS: a layer must contain independent generator orders"); fi;
        else layer:=rec(orders:=[]); fi;
        layer.name:=name; Add(layers,layer);
    od;
    KOAHSS_ExtensionPresentationSerial.next:=KOAHSS_ExtensionPresentationSerial.next+1;
    presentationPrefix:=Concatenation("extension:",String(KOAHSS_ExtensionPresentationSerial.next),":");
    ids:=[]; matrix:=[]; vectors:=[]; stages:=[];
    smith:=KOAHSS_ExtensionSmith(matrix,ids);
    model:=rec(generatorCount:=0,generatorIds:=ids,relationMatrix:=matrix,
        layers:=[],presentationId:=Concatenation(presentationPrefix,"0"),smith:=smith);
    for layer in layers do
        reason:=fail;
        if IsBound(layer.status) and layer.status<>"computed" then
            reason:="associated-graded layer unresolved";
            if IsBound(layer.reason) then reason:=layer.reason; fi;
        elif not IsBound(layer.orders) or not IsList(layer.orders)
            or not ForAll(layer.orders,o->IsInt(o) and (o=0 or o>1)) then
            Error("koAHSS: layer orders must be zero (free) or integers greater than one");
        fi;
        if reason<>fail then
            return rec(status:="unresolved",reason:=reason,pendingLayer:=layer.name,
                completedStages:=stages,extensionVectors:=vectors,lowerModel:=model);
        fi;
        oldCount:=Length(ids); count:=Length(layer.orders);
        layer.startColumn:=oldCount+1; layer.endColumn:=oldCount+count;
        newIds:=List([1..count],i->Concatenation(layer.name,":",String(i)));
        previous:=smith;
        # Query in the immutable common lower presentation, before adding any
        # generators from this layer; correlations between all rows survive.
        for i in [1..count] do
            order:=layer.orders[i];
            if order<>0 then
                if oldCount=0 then
                    response:=rec(status:="computed",lowerPresentationId:=model.presentationId,
                        lowerCoordinates:=[],witness:=rec(kind:="zero-lower-group"));
                elif oracle=fail then
                    response:=rec(status:="unresolved",reason:="no relation-vector oracle for this layer");
                else response:=oracle(layer,i,order,model); fi;
                if not IsRecord(response) or not IsBound(response.status) then
                    Error("koAHSS: malformed extension oracle response");
                fi;
                if response.status<>"computed" then
                    if response.status<>"unresolved" or not IsBound(response.reason) then
                        Error("koAHSS: invalid extension oracle status");
                    fi;
                    return rec(status:="unresolved",reason:=response.reason,
                        pendingLayer:=layer.name,pendingGenerator:=i,completedStages:=stages,
                        extensionVectors:=vectors,lowerModel:=model,pendingRelation:=response);
                fi;
                if not IsBound(response.lowerPresentationId)
                    or response.lowerPresentationId<>model.presentationId then
                    Error("koAHSS: extension oracle returned a different lower basis");
                fi;
                if not IsBound(response.lowerCoordinates)
                    or not IsList(response.lowerCoordinates)
                    or Length(response.lowerCoordinates)<>oldCount
                    or not ForAll(response.lowerCoordinates,IsInt)
                    or not IsBound(response.witness) then
                    Error("koAHSS: relation vector needs integer lower coordinates and a witness");
                fi;
                Add(vectors,rec(layer:=layer.name,generatorId:=newIds[i],order:=order,
                    lowerGeneratorIds:=ShallowCopy(ids),result:=response));
            fi;
        od;
        matrix:=List(matrix,r->Concatenation(r,List([1..count],j->0)));
        for response in Filtered(vectors,v->v.layer=layer.name) do
            row:=Concatenation(-response.result.lowerCoordinates,List([1..count],j->0));
            i:=Position(newIds,response.generatorId); row[oldCount+i]:=response.order;
            Add(matrix,row);
        od;
        Append(ids,newIds);
        smith:=KOAHSS_ExtensionSmith(matrix,ids);
        inclusion:=List(previous.basis.expressions,r->Concatenation(r,List([1..count],j->0)));
        if not IsEmpty(inclusion) then inclusion:=List(inclusion,r->(r*smith.V){smith.activeIndices}); fi;
        projection:=List(smith.basis.expressions,r->r{[oldCount+1..oldCount+count]});
        stage:=rec(layer:=layer.name,group:=smith.group,invariants:=smith.invariants,
            relationMatrix:=StructuralCopy(matrix),smith:=smith,
            inclusionMatrix:=inclusion,quotientMatrix:=projection,quotientOrders:=layer.orders);
        Add(stages,stage); Add(model.layers,layer);
        model:=rec(generatorCount:=Length(ids),generatorIds:=ShallowCopy(ids),
            relationMatrix:=StructuralCopy(matrix),layers:=ShallowCopy(model.layers),
            presentationId:=Concatenation(presentationPrefix,String(Length(stages))),smith:=smith);
    od;
    return rec(status:="computed",group:=smith.group,invariants:=smith.invariants,
        basis:=smith.basis,relationMatrix:=matrix,extensionVectors:=vectors,
        smith:=smith,filtration:=stages,lowerModel:=model);
end);

BindGlobal("KOAHSS_ExtensionLayers",function(context,degree)
    local result, name, q, p, cell, layer, generators, data;
    result:=rec();
    for name in ["D","C","B","A"] do
        q:=rec(D:=-4,C:=-2,B:=-1,A:=0).(name); p:=degree-q-3;
        if p<0 then
            layer:=rec(name:=name,p:=p,q:=q,orders:=[],generators:=[],cochains:=[]);
        else
            cell:=context.getCell(6,p,q);
            if KOAHSS_IsUnresolved(cell) then
                layer:=rec(name:=name,p:=p,q:=q,status:="unresolved",cell:=cell);
            else
                generators:=IndependentGeneratorsOfAbelianGroup(cell.group);
                data:=context.backend.cohomologyData(p,q);
                layer:=rec(name:=name,p:=p,q:=q,cell:=cell,group:=cell.group,
                    generators:=generators,orders:=List(generators,function(g)
                        if Order(g)=infinity then return 0; else return Order(g); fi;
                    end),cochains:=List(generators,g->data.represent(cell.lift(g))));
            fi;
        fi;
        result.(name):=layer;
    od;
    return result;
end);

InstallGlobalFunction(koFull,function(arg)
    local ahss, context, degrees, results, invariants, degree, layers, oracle,
        result, status, model, computeDegree, runDegree, options, requested,
        models,closeModels,getModel,selected,transferAttempt,selection;
    options:=rec();
    if Length(arg) in [4,5] then
        if Length(arg)=5 then options:=arg[5]; fi;
    elif Length(arg) in [1,2] and IsRecord(arg[1]) then
        if Length(arg)=2 then options:=arg[2]; fi;
    else
        Error("usage: koFull(group or HAP resolution,s,omega,k[,options]) or koFull(detailedE6Result[,options])");
    fi;
    if not IsRecord(options) or ForAny(RecNames(options),name->name<>"extensionModel") then
        Error("koFull: options must be a record containing only extensionModel");
    fi;
    requested:="bar";
    if IsBound(options.extensionModel) then requested:=options.extensionModel; fi;
    if not IsString(requested) or not requested in ["bar","transfer"] then
        Error("koFull: extensionModel must be bar or transfer");
    fi;
    if Length(arg) in [4,5] then
        ahss:=koAHSS(arg[1],arg[2],arg[3],arg[4],rec(details:=true));
    else ahss:=arg[1]; fi;
    if not IsBound(ahss.kind) or ahss.kind<>"koAHSSResult"
        or not IsBound(ahss.computedThrough) or ahss.computedThrough<>6
        or not IsBound(ahss._context) or not IsBound(ahss._context.getCell)
        or not IsBound(ahss._context.backend.cohomologyData)
        or not IsBound(ahss.pages) or not 6 in ahss.pages.pageNumbers then
        Error("koFull: a detailed E6 result with retained cochain context is required");
    fi;
    context:=ahss._context; degrees:=[-1..ahss.maxDegree]; results:=[]; invariants:=[];
    model:=fail; models:=rec(bar:=rec(),transfer:=rec());
    closeModels:=function()
        local kind,key,stored;
        for kind in ["bar","transfer"] do
            for key in RecNames(models.(kind)) do
                stored:=models.(kind).(key);
                if IsBound(stored.close) then stored.close(); fi;
            od;
        od;
    end;
    getModel:=function(kind)
        local key,factory;
        key:=String(degree);
        # The complete bar worker already supports all later degrees and
        # retains expensive universal-formula caches across those degrees.
        if kind="bar" then key:="shared"; fi;
        if not IsBound(models.(kind).(key)) then
            factory:="KOAHSS_ExtensionBarModel";
            if kind="transfer" then factory:="KOAHSS_ExtensionTransferredModel"; fi;
            if not IsBoundGlobal(factory) then
                models.(kind).(key):=rec(status:="unresolved",
                    reason:="the requested extension model is unavailable");
            else
                models.(kind).(key):=CallFuncList(ValueGlobal(factory),[context.backend,degree]);
            fi;
        fi;
        return models.(kind).(key);
    end;
    computeDegree:=function(kind)
        local candidate,certification;
        layers:=KOAHSS_ExtensionLayers(context,degree);
        if degree in [3..6] then
            model:=getModel(kind);
            if IsBound(model.lastFailure) then Unbind(model.lastFailure); fi;
            if model.status="computed" and model.supports(degree) then
                oracle:=CallFuncList(ValueGlobal("KOAHSS_ExtensionHigherOracle"),
                    [context.backend,degree,layers,model]);
            else
                if kind="transfer" then
                    if IsBound(model.reason) then
                        return rec(status:="unresolved",reason:=model.reason,pendingLayer:="model-setup");
                    fi;
                    return rec(status:="unresolved",reason:=
                        "transfer extension resource budget exceeded in this degree",pendingLayer:="model-setup");
                fi;
                # Preserve the reference path: zero layers and free quotient
                # layers need no torsion query, so setup refusal alone cannot
                # turn a previously computed zero/split group into unknown.
                oracle:=function(l,i,m,h)
                    if IsBound(model.reason) then return rec(status:="unresolved",reason:=model.reason); fi;
                    return rec(status:="unresolved",reason:="complete bar extension resource budget exceeded in this degree");
                end;
            fi;
        else
            oracle:=CallFuncList(ValueGlobal("KOAHSS_StackingExtensionOracle"),
                [context.backend,degree,layers]);
        fi;
        candidate:=koAHSSExtensionFromLayers(layers,oracle);
        if degree in [3..6] and candidate.status="computed"
            and model<>fail and model.status="computed" and model.supports(degree) then
            candidate.algebraAudit:=CallFuncList(ValueGlobal("KOAHSS_ExtensionFiniteAudit"),
                [model,degree,layers,candidate]);
            if candidate.algebraAudit.status<>"computed" then
                return rec(status:="unresolved",reason:=candidate.algebraAudit.reason,
                    pendingLayer:="quotient-audit",candidatePresentation:=candidate);
            fi;
        fi;
        if kind="transfer" and candidate.status="computed" then
            candidate.certificateLevel:="transfer-R";
            if not IsBound(model.classEquivalenceVerified) or model.classEquivalenceVerified<>true then
                certification:=rec(status:="unresolved",
                    reason:="the transferred presentation has no verified bar-equivalence certificate");
                if IsBound(model.certifyPresentation) then
                    certification:=model.certifyPresentation(degree,layers,candidate);
                fi;
                candidate.barCertification:=certification;
                if certification.status<>"computed" then
                    return rec(status:="unresolved",reason:=certification.reason,
                        pendingLayer:="transfer-certification",candidatePresentation:=candidate);
                fi;
                candidate.certificateLevel:="complete-bar-certified-transfer";
            fi;
        fi;
        return candidate;
    end;
    runDegree:=function(kind)
        local attempt,oldBreak,answer;
        model:=fail;
        oldBreak:=BreakOnError; BreakOnError:=false;
        attempt:=CALL_WITH_CATCH(computeDegree,[kind]);
        BreakOnError:=oldBreak;
        if attempt[1] then answer:=attempt[2];
        else
            if model<>fail and IsBound(model.lastFailure) and model.lastFailure.status="unresolved" then
                if IsBound(model.close) then model.close(); fi;
                answer:=rec(status:="unresolved",reason:=model.lastFailure.reason,pendingLayer:="resource-limit");
            else
                closeModels();
                Error("koFull: exact extension calculation failed; the original error is reported above");
            fi;
        fi;
        answer.degree:=degree; answer.layers:=layers;
        if model<>fail and IsBound(model.modelId) then answer.modelId:=model.modelId; fi;
        if kind="transfer" and model<>fail then
            if IsBound(model.close) then model.close(); fi;
            # Transfer models have degree-specific supports and are never
            # reused in a later degree. Release their worker and RPC caches.
            Unbind(models.transfer.(String(degree)));
        fi;
        return answer;
    end;
    for degree in degrees do
        selected:="low";
        if degree in [3..6] then selected:="bar"; fi;
        selection:=rec(requested:=requested,selected:=selected,fallback:=false);
        if requested="transfer" and degree in [3..5] then
            transferAttempt:=runDegree("transfer");
            if transferAttempt.status="computed" then
                result:=transferAttempt; selection.selected:="transfer";
            else
                result:=runDegree("bar"); result.transferAttempt:=transferAttempt;
                selection.fallback:=true; selection.reason:=transferAttempt.reason;
            fi;
        else
            result:=runDegree(selected);
            if requested="transfer" and degree=6 then
                selection.fallback:=true;
                selection.reason:="degree six uses the complete bar section";
            fi;
        fi;
        result.modelSelection:=selection;
        Add(results,result);
        if result.status="computed" then Add(invariants,result.invariants);
        else Add(invariants,rec(status:=result.status,reason:=result.reason)); fi;
    od;
    closeModels();
    status:="unresolved";
    if ForAll(results,r->r.status="computed") then status:="computed";
    elif ForAny(results,r->r.status="computed") then status:="partial"; fi;
    return rec(kind:="koFullResult",maxDegree:=ahss.maxDegree,degrees:=degrees,
        invariants:=invariants,degreeResults:=results,ahss:=ahss,
        pages:=ahss.pages,extensionModel:=requested,
        status:=status,scope:="five-row-stacking-model",certified_ko:=false);
end);
