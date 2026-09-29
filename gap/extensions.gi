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

# Localization at the primes (doc/extensions.md, "Localization at the primes").
# The prime localization of the relation measurement is on unless
# FERMIONAHSS_PRIME_LOCAL=0.
BindGlobal("KOAHSS_PrimeLocalizationEnabled",function()
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_PRIME_LOCAL) then
        return GAPInfo.SystemEnvironment.FERMIONAHSS_PRIME_LOCAL<>"0";
    fi;
    return true;
end);

# The prime of a relation of prime-power order, or fail.
BindGlobal("KOAHSS_RelationPrime",function(order)
    if not IsInt(order) or order<2 or not IsPrimePowerInt(order) then return fail; fi;
    return SmallestRootInt(order);
end);

# The stacking model measuring the relation of a generator of the given order
# in the layer name of package degree k: "complete" (the transferred model),
# "three-local" (the two-layer model of the prime three in degrees five and
# six) or "split" (no k-invariant of ko links the rows q=0 and q=-4 at the
# prime: the primes at least five in every degree, and the prime three below
# degree five, where P^1 vanishes). Only A generators have odd-primary relations.
BindGlobal("KOAHSS_ExtensionRelationModel",function(k,name,order)
    local prime;
    prime:=KOAHSS_RelationPrime(order);
    if order=0 or prime=fail or prime=2 or name<>"A" or not KOAHSS_PrimeLocalizationEnabled() then
        return rec(model:="complete",prime:=prime);
    fi;
    if prime>=5 or k<=4 then return rec(model:="split",prime:=prime); fi;
    return rec(model:="three-local",prime:=prime);
end);

BindGlobal("KOAHSS_ExtensionSplitResponse",function(order,prime,lower)
    local reason;
    if prime=3 then reason:="P^1 vanishes on classes of degree below two";
    else reason:=Concatenation("the rows q=0 and q=-4 lie in different Adams summands of ko localized at ",String(prime)); fi;
    return rec(status:="computed",lowerPresentationId:=lower.presentationId,
        lowerCoordinates:=List([1..lower.generatorCount],j->0),
        witness:=rec(operation:="split",power:=order,prime:=prime,model:="split",
            measuredLayers:=[],truncatedBelow:=fail,modelId:="prime-split",
            certificate:=Concatenation("no stacking correction at the prime ",String(prime),": ",reason,
                "; the relation is ",String(order),"*g=0")));
end);

# Per prime: the number of measured relations and the models that measured them.
BindGlobal("KOAHSS_ExtensionPrimeSummary",function(result)
    local summary,v,w,prime,model,entry;
    summary:=[];
    if not IsBound(result.extensionVectors) then return summary; fi;
    for v in result.extensionVectors do
        w:=v.result.witness;
        # A generator over the zero lower group measures nothing.
        if IsBound(w.kind) and w.kind="zero-lower-group" then continue; fi;
        if IsBound(w.prime) then prime:=w.prime; else prime:=KOAHSS_RelationPrime(v.order); fi;
        if IsBound(w.model) then model:=w.model; else model:="complete"; fi;
        entry:=First(summary,e->e.prime=prime);
        if entry=fail then entry:=rec(prime:=prime,relations:=0,models:=[]); Add(summary,entry); fi;
        entry.relations:=entry.relations+1; AddSet(entry.models,model);
    od;
    return summary;
end);

# The generators of the resolved layers (those below the first unresolved
# layer), grouped by the prime of their order, right after the layers are
# read. Layer generators have prime-power or infinite order. A free generator
# (order 0) belongs to no prime: it is a lower generator of the relations of
# every prime, and carries no relation itself.
BindGlobal("KOAHSS_ExtensionPrimeParts",function(layers)
    local names,name,orders,free,primes,i,prime,parts,part;
    names:=[];
    for name in ["D","C","B","A"] do
        if IsBound(layers.(name).status) then break; fi;
        Add(names,name);
    od;
    free:=rec(); primes:=[];
    for name in names do
        orders:=layers.(name).orders;
        free.(name):=Filtered([1..Length(orders)],i->orders[i]=0);
        for i in [1..Length(orders)] do
            if orders[i]<>0 then
                prime:=KOAHSS_RelationPrime(orders[i]);
                if prime=fail then
                    Error("koFull: an E6 generator order is not a prime power: ",orders[i]);
                fi;
                AddSet(primes,prime);
            fi;
        od;
    od;
    parts:=[];
    for prime in primes do
        part:=rec(prime:=prime,generators:=rec());
        for name in names do
            orders:=layers.(name).orders;
            part.generators.(name):=Filtered([1..Length(orders)],
                i->orders[i]<>0 and KOAHSS_RelationPrime(orders[i])=prime);
        od;
        Add(parts,part);
    od;
    return rec(names:=names,free:=free,parts:=parts);
end);

# The relation oracle that returns recorded rows, in the presentation of the
# lower layers that koAHSSExtensionFromLayers passes to it.
BindGlobal("KOAHSS_ExtensionReplayOracle",function(rows)
    return function(layer,index,order,lower)
        local recorded,response;
        if not IsBound(rows.(layer.name)[index]) then
            return rec(status:="unresolved",
                reason:="the relation was not measured: a relation below it is unresolved");
        fi;
        recorded:=rows.(layer.name)[index];
        if recorded.status<>"computed" then return recorded; fi;
        if Length(recorded.lowerCoordinates)<>lower.generatorCount then
            Error("koFull: a recorded relation does not match its lower presentation");
        fi;
        response:=ShallowCopy(recorded);
        response.lowerPresentationId:=lower.presentationId;
        return response;
    end;
end);

# The lower presentation of the layers names (from D upward) with their
# recorded relation rows.
BindGlobal("KOAHSS_ExtensionLowerPresentation",function(layers,names,oracle)
    local input,name,result;
    input:=rec();
    for name in names do input.(name):=layers.(name); od;
    result:=koAHSSExtensionFromLayers(input,oracle);
    if result.status<>"computed" then
        Error("koFull: a lower presentation needs every recorded lower relation");
    fi;
    return result.lowerModel;
end);

# The relation rows of a degree, prime by prime (two first) and, for each
# prime, one generator at a time from D upward. A relation of the prime p is
# measured only if it does not split and some lower generator has p-power or
# infinite order: otherwise its row is zero, because the lower group has no
# p-primary part and no free part, so it lies in m*H for the order m.
# With prime localization on, a row keeps only the coordinates on
# generators of its prime and free generators. The other coordinates are on
# torsion of order prime to m that no kept row refers to (the odd torsion of
# the rows B, C and A lies in D, below every two-primary relation, and an odd
# relation belongs to A and has no B or C coordinate), so they vanish after
# localization at p and the rows of all primes together present the group.
BindGlobal("KOAHSS_ExtensionPrimeRows",function(engine,layers,parts,degree)
    local rows,restrict,failure,replay,part,prime,position,name,below,fullLower,
        localLower,lower,i,order,kind,response,failed,keep,j,n,o,coordinates;
    rows:=rec(D:=[],C:=[],B:=[],A:=[]);
    restrict:=KOAHSS_PrimeLocalizationEnabled();
    failure:=Length(parts.names)+1;
    replay:=KOAHSS_ExtensionReplayOracle(rows);
    for part in parts.parts do
        prime:=part.prime; failed:=false;
        for position in [1..Length(parts.names)] do
            if failed or position>failure then break; fi;
            name:=parts.names[position]; below:=parts.names{[1..position-1]};
            fullLower:=Sum(below,n->Length(layers.(n).orders));
            # Over the zero lower group there is no relation to record.
            if fullLower=0 or IsEmpty(part.generators.(name)) then continue; fi;
            keep:=Concatenation(List(below,n->List(layers.(n).orders,
                o->o=0 or KOAHSS_RelationPrime(o)=prime)));
            localLower:=Number(keep,x->x);
            lower:=fail;
            for i in part.generators.(name) do
                order:=layers.(name).orders[i];
                kind:=KOAHSS_ExtensionRelationModel(degree,name,order);
                if restrict and kind.model<>"split" and localLower=0 then
                    response:=rec(status:="computed",lowerCoordinates:=List([1..fullLower],j->0),
                        witness:=rec(kind:="zero-local-lower-group",operation:="xtimes",
                            power:=order,prime:=prime,model:=kind.model,measuredLayers:=[],
                            certificate:="the lower group has no generator of this prime and no free generator"));
                else
                    if lower=fail then lower:=KOAHSS_ExtensionLowerPresentation(layers,below,replay); fi;
                    response:=engine.answer(layers.(name),i,order,lower);
                    if restrict and response.status="computed" and not ForAll(keep,x->x) then
                        coordinates:=List([1..fullLower],j->0);
                        for j in [1..fullLower] do
                            if keep[j] then coordinates[j]:=response.lowerCoordinates[j]; fi;
                        od;
                        if coordinates<>response.lowerCoordinates then
                            response:=ShallowCopy(response);
                            response.witness:=ShallowCopy(response.witness);
                            response.witness.fullLowerCoordinates:=response.lowerCoordinates;
                            response.witness.localizedAt:=prime;
                            response.lowerCoordinates:=coordinates;
                        fi;
                    fi;
                fi;
                rows.(name)[i]:=response;
                if response.status<>"computed" then
                    failed:=true; failure:=Minimum(failure,position); break;
                fi;
            od;
        od;
    od;
    return rows;
end);

# A detailed koAHSS_batch result computed through E6, with its cochain context.
BindGlobal("KOAHSS_DetailedE6Context",function(ahss)
    if not IsBound(ahss.kind) or ahss.kind<>"koAHSSResult"
        or not IsBound(ahss.computedThrough) or ahss.computedThrough<>6
        or not IsBound(ahss._context) or not IsBound(ahss._context.getCell)
        or not IsBound(ahss._context.backend.cohomologyData)
        or not IsBound(ahss.pages) or not 6 in ahss.pages.pageNumbers then
        Error("koFull: a detailed E6 result with retained cochain context is required");
    fi;
    return ahss._context;
end);

# The extension problem of one degree: its E6 layers A (degree-3,0),
# B (degree-2,-1), C (degree-1,-2) and D (degree+1,-4). Right after the layers
# are read, their generators are grouped by prime (KOAHSS_ExtensionPrimeParts);
# the relations of each prime are recorded from D upward, one generator at a
# time (KOAHSS_ExtensionPrimeRows), and the rows of all primes are assembled
# into one presentation. The stacking model, its worker, the flat lifts and
# the bar transport are built only when a relation is measured.
BindGlobal("KOAHSS_FullDegree",function(context,degree)
    local layers,engine,model,computeDegree,attempt,oldBreak,answer;
    layers:=fail; engine:=fail;
    computeDegree:=function()
        local parts,rows,candidate;
        layers:=KOAHSS_ExtensionLayers(context,degree);
        if degree<1 then
            # Degrees -1 and 0 have the single layer D: no relation is measured.
            return koAHSSExtensionFromLayers(layers,fail);
        fi;
        parts:=KOAHSS_ExtensionPrimeParts(layers);
        # Degrees 1-6: the native transferred model on R, when needed.
        engine:=CallFuncList(ValueGlobal("KOAHSS_ExtensionRelationEngine"),
            [context.backend,degree,layers,function()
                if not IsBoundGlobal("KOAHSS_ExtensionTransferredModel") then
                    return rec(status:="unresolved",reason:="the native extension model is unavailable");
                fi;
                return CallFuncList(ValueGlobal("KOAHSS_ExtensionTransferredModel"),
                    [context.backend,degree]);
            end,rec(layerLimited:=true)]);
        rows:=KOAHSS_ExtensionPrimeRows(engine,layers,parts,degree);
        candidate:=koAHSSExtensionFromLayers(layers,KOAHSS_ExtensionReplayOracle(rows));
        if candidate.status<>"computed" then
            if IsBound(candidate.pendingRelation) and IsBound(candidate.pendingRelation.code)
               and candidate.pendingRelation.code="model-setup" then
                return rec(status:="unresolved",reason:=candidate.reason,pendingLayer:="model-setup");
            fi;
            return candidate;
        fi;
        if engine.model()<>fail then
            # Gauge completeness is an assumption of this transferred model,
            # as are commutativity and associativity of stacking on gauge
            # classes. No finite multiplication table is audited; the
            # accepted relation equations are checked on R.
            candidate.certificateLevel:="transfer-R";
            candidate.gaugeCompletenessAssumed:=true;
            candidate.abelianQuotientAssumed:=true;
        elif ForAny(candidate.extensionVectors,v->IsBound(v.result.witness.model)
                and v.result.witness.model="split") then
            candidate.certificateLevel:="prime-split";
        else
            # No relation needed a measurement: every row is zero by the
            # orders of the lower layers or by a divisibility certificate,
            # and the group is the direct sum of its layers.
            candidate.certificateLevel:="direct-sum";
        fi;
        candidate.primes:=KOAHSS_ExtensionPrimeSummary(candidate);
        candidate.primeParts:=List(parts.parts,part->rec(prime:=part.prime,
            generators:=part.generators,invariants:=Filtered(candidate.invariants,
                x->x<>0 and KOAHSS_RelationPrime(x)=part.prime)));
        if Number(["D","C","B","A"],name->not IsEmpty(layers.(name).orders))<=1 then
            candidate.singleLayer:=true;
        fi;
        return candidate;
    end;
    oldBreak:=BreakOnError; BreakOnError:=false;
    attempt:=CALL_WITH_CATCH(computeDegree,[]);
    BreakOnError:=oldBreak;
    model:=fail;
    if engine<>fail then model:=engine.builtModel(); engine.close(); fi;
    if attempt[1] then answer:=attempt[2];
    elif model<>fail and IsRecord(model) and IsBound(model.lastFailure)
         and model.lastFailure.status="unresolved" then
        answer:=rec(status:="unresolved",reason:=model.lastFailure.reason,pendingLayer:="resource-limit");
    else
        Error("koFull: exact extension calculation failed; the original error is reported above");
    fi;
    answer.degree:=degree; answer.layers:=layers;
    if engine<>fail and engine.model()<>fail and IsBound(engine.model().modelId) then
        answer.modelId:=engine.model().modelId;
    fi;
    return answer;
end);

InstallGlobalFunction(koFull_batch,function(arg)
    local ahss, context, degrees, results, invariants, degree, result, status;
    if Length(arg)=4 then
        ahss:=koAHSS_batch(arg[1],arg[2],arg[3],arg[4],rec(details:=true));
    elif Length(arg)=1 and IsRecord(arg[1]) then
        ahss:=arg[1];
    else
        Error("usage: koFull_batch(group or HAP resolution,s,omega,k) or koFull_batch(detailedE6Result)");
    fi;
    context:=KOAHSS_DetailedE6Context(ahss);
    degrees:=[-1..ahss.maxDegree]; results:=[]; invariants:=[];
    for degree in degrees do
        # Each native worker supports one degree and is released afterwards.
        result:=KOAHSS_FullDegree(context,degree);
        Add(results,result);
        if result.status="computed" then Add(invariants,result.invariants);
        else Add(invariants,rec(status:=result.status,reason:=result.reason)); fi;
    od;
    status:="unresolved";
    if ForAll(results,r->r.status="computed") then status:="computed";
    elif ForAny(results,r->r.status="computed") then status:="partial"; fi;
    return rec(kind:="koFullResult",maxDegree:=ahss.maxDegree,degrees:=degrees,
        invariants:=invariants,degreeResults:=results,ahss:=ahss,
        pages:=ahss.pages,gaugeCompletenessAssumed:=true,
        status:=status,scope:="five-row-stacking-model",certified_ko:=false);
end);

# Degree k only: the pages through E6, then the extension problem of degree k.
InstallGlobalFunction(koFull,function(arg)
    local ahss, context, k, answer, result;
    if Length(arg)=4 then
        ahss:=koAHSS_batch(arg[1],arg[2],arg[3],arg[4],rec(details:=true));
    elif Length(arg)=1 and IsRecord(arg[1]) then
        ahss:=arg[1];
    else
        Error("usage: koFull(group or HAP resolution,s,omega,k) or koFull(detailedE6Result)");
    fi;
    context:=KOAHSS_DetailedE6Context(ahss);
    k:=ahss.maxDegree;
    answer:=KOAHSS_FullDegree(context,k);
    result:=rec(kind:="koFullDegreeResult",k:=k,status:=answer.status,
        line:=rec(kind:="koAHSSLine",k:=k,pageNumbers:=[6],lines:=[KOAHSS_PageLine(
            ahss.pages.tables[Position(ahss.pages.pageNumbers,6)],k)]),
        degreeResult:=answer,scope:="five-row-stacking-model",certified_ko:=false);
    if answer.status="computed" then result.invariants:=answer.invariants;
    else result.reason:=answer.reason; fi;
    return result;
end);
