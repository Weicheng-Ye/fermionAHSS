# Measure powers of fixed flat representatives, retaining every boundary carry.
# Every division below verifies an actual ordered cochain equation. No lower
# coordinate is identified merely from the order of a stacked element.
BindGlobal("KOAHSS_ExtensionStateIsZero",state->ForAll(["A","B","C","D"],
    name->ForAll(state.(name),x->x=0)));

# A state measured only through the layer index upto (A=0 to D=3) has zeros
# above it; only its measured layers are compared.
BindGlobal("KOAHSS_ExtensionStateIsZeroThrough",function(state,upto)
    return ForAll(["A","B","C","D"]{[1..upto+1]},name->ForAll(state.(name),x->x=0));
end);

BindGlobal("KOAHSS_ExtensionStatesAgreeThrough",function(x,y,upto)
    return ForAll(["A","B","C","D"]{[1..upto+1]},name->x.(name)=y.(name));
end);

# The target-layer shortcut of the relation measurement is on unless
# FERMIONAHSS_LAYERED_RELATIONS=0.
BindGlobal("KOAHSS_LayeredRelationsEnabled",function()
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_LAYERED_RELATIONS) then
        return GAPInfo.SystemEnvironment.FERMIONAHSS_LAYERED_RELATIONS<>"0";
    fi;
    return true;
end);

# Every lower generator e of the recorded presentation with e not in n*G_low
# is free for a relation of order n: its component of the relation changes
# the class in Ext(Z/n,G_low)=G_low/nG_low. The target layer is the lowest
# layer with a free generator; the components of the relation in the layers
# below it cannot change the group, and the certificate lists, for every
# generator below the target layer, the integer combination showing that it
# lies in n*G_low modulo the relations.
BindGlobal("KOAHSS_ExtensionTargetLayer",function(lower,n)
    local fields,m,rows,j,unit,solution,free,certificate,layerOf,position,entry,check,combination;
    fields:=["A","B","C","D"];
    m:=lower.generatorCount;
    if m=0 then return rec(layer:=fail,index:=fail,freeGenerators:=[],certificate:=[]); fi;
    if not IsInt(n) or n<2 then Error("koFull: a relation order must be an integer greater than one"); fi;
    rows:=Concatenation(n*IdentityMat(m),lower.relationMatrix);
    layerOf:=function(column)
        local stored;
        for stored in lower.layers do
            if column>=stored.startColumn and column<stored.startColumn+Length(stored.orders) then
                return stored.name;
            fi;
        od;
        Error("koFull: a lower generator lies in no recorded layer");
    end;
    free:=[]; certificate:=[];
    for j in [1..m] do
        unit:=List([1..m],i->0); unit[j]:=1;
        solution:=koAHSSSolveIntegerSystem(rows,unit);
        if solution=fail then Add(free,j);
        else
            combination:=solution.particular{[1..m]};
            check:=n*combination;
            if Length(rows)>m then check:=check+solution.particular{[m+1..Length(rows)]}*lower.relationMatrix; fi;
            if check<>unit then Error("koFull: the target-layer certificate failed its exact equation"); fi;
            Add(certificate,rec(generatorId:=lower.generatorIds[j],column:=j,layer:=layerOf(j),
                multiple:=n,combination:=combination,
                relationCoefficients:=solution.particular{[m+1..Length(rows)]}));
        fi;
    od;
    if IsEmpty(free) then
        return rec(layer:=fail,index:=fail,freeGenerators:=[],certificate:=certificate);
    fi;
    position:=Maximum(List(free,j->Position(fields,layerOf(j))));
    return rec(layer:=fields[position],index:=position-1,
        freeGenerators:=List(free,j->lower.generatorIds[j]),
        certificate:=Filtered(certificate,entry->Position(fields,entry.layer)>position));
end);

BindGlobal("KOAHSS_ExtensionDivideLeft",function(arg)
    local model,k,left,total,upto,right,product,name;
    model:=arg[1]; k:=arg[2]; left:=arg[3]; total:=arg[4];
    upto:=3; if Length(arg)>=5 then upto:=arg[5]; fi;
    if IsBound(model.divideLeft) then
        if not IsFunction(model.divideLeft) then
            Error("koFull: divideLeft must be a function");
        fi;
        if upto<3 then
            right:=model.divideLeft(k,left,total,upto);
            if not KOAHSS_ExtensionStatesAgreeThrough(model.xtimes(k,left,right,upto),total,upto) then
                Error("koFull: transferred stacking division failed its exact equation");
            fi;
            if not KOAHSS_ExtensionStateIsZeroThrough(model.d(k,right,upto),upto) then
                Error("koFull: transferred division produced a nonflat state");
            fi;
            MakeImmutable(right); return right;
        fi;
        right:=model.divideLeft(k,left,total);
        if model.xtimes(k,left,right)<>total then
            Error("koFull: transferred stacking division failed its exact equation");
        fi;
        if not KOAHSS_ExtensionStateIsZero(model.d(k,right)) then
            Error("koFull: transferred division produced a nonflat state");
        fi;
        MakeImmutable(right); return right;
    fi;
    if upto<3 then Error("koFull: a layer-limited division needs the transferred model"); fi;
    right:=StructuralCopy(model.zero(k)); right.A:=total.A-left.A;
    product:=model.xtimes(k,left,right); right.B:=List(total.B-product.B,x->x mod 2);
    product:=model.xtimes(k,left,right); right.C:=List(total.C-product.C,x->x mod 2);
    product:=model.xtimes(k,left,right); right.D:=total.D-product.D;
    product:=model.xtimes(k,left,right);
    if product<>total then Error("koFull: triangular stacking division failed its exact equation"); fi;
    if not KOAHSS_ExtensionStateIsZero(model.d(k,right)) then
        Error("koFull: division of flat states produced a nonflat state");
    fi;
    MakeImmutable(right); return right;
end);

BindGlobal("KOAHSS_ExtensionPower",function(arg)
    local model,k,state,m,upto,answer,base;
    model:=arg[1]; k:=arg[2]; state:=arg[3]; m:=arg[4];
    upto:=3; if Length(arg)>=5 then upto:=arg[5]; fi;
    if not IsInt(m) then Error("koFull: a stacking power must be integral"); fi;
    if m<0 then
        state:=KOAHSS_ExtensionDivideLeft(model,k,state,model.zero(k),upto); m:=-m;
    fi;
    answer:=model.zero(k); base:=state;
    while m>0 do
        if m mod 2=1 then
            if KOAHSS_ExtensionStateIsZero(answer) then answer:=base;
            elif upto<3 then answer:=model.xtimes(k,answer,base,upto);
            else answer:=model.xtimes(k,answer,base); fi;
        fi;
        m:=QuoInt(m,2);
        if m>0 then
            if upto<3 then base:=model.xtimes(k,base,base,upto);
            else base:=model.xtimes(k,base,base); fi;
        fi;
    od;
    if upto<3 then
        if not KOAHSS_ExtensionStateIsZeroThrough(model.d(k,answer,upto),upto) then
            Error("koFull: the power of a flat representative is not flat");
        fi;
    elif not KOAHSS_ExtensionStateIsZero(model.d(k,answer)) then
        Error("koFull: the power of a flat representative is not flat");
    fi;
    return answer;
end);

# Compare with marked lower normal forms when ordinary coboundaries do not
# generate the E6 equivalence relation. Every candidate uses the original
# complete lower lifts, and every accepted vector has one full gauge witness.
BindGlobal("KOAHSS_ExtensionGaugeReduce",function(arg)
    local model,k,state,layers,lower,preferred,ranges,stored,i,count,choices,
        canonicals,coordinates,canonical,power,j,stage,stageIndex,stages,
        proof,attempts,used,budget,allowance,remaining,coordinateMethod;
    if not Length(arg) in [5,6] then Error("koFull: invalid gauge reduction arguments"); fi;
    model:=arg[1]; k:=arg[2]; state:=arg[3]; layers:=arg[4]; lower:=arg[5];
    preferred:=fail; if Length(arg)=6 then preferred:=arg[6]; fi;
    if preferred<>fail then
        if not IsList(preferred) or Length(preferred)<>lower.generatorCount
            or not ForAll(preferred,IsInt) then
            Error("koFull: projected lower coordinates have the wrong marked basis");
        fi;
        choices:=[ShallowCopy(preferred)]; coordinateMethod:="E6 projection of the final D residual";
    else
        ranges:=[];
        for stored in lower.layers do
            for i in [1..Length(stored.orders)] do
                if stored.orders[i]=0 then
                    return rec(status:="unresolved",reason:="gauge reduction needs projected coordinates for a free lower generator");
                fi;
                Add(ranges,[0..stored.orders[i]-1]);
            od;
        od;
        count:=Product(List(ranges,Length));
        if count>32 then return rec(status:="unresolved",
            reason:="gauge reduction exceeds 32 marked lower normal forms"); fi;
        choices:=Cartesian(ranges); coordinateMethod:="marked finite lower normal forms";
    fi;
    canonicals:=[];
    for coordinates in choices do
        canonical:=model.zero(k);
        for stored in lower.layers do
            for i in [1..Length(stored.orders)] do
                j:=stored.startColumn+i-1;
                power:=KOAHSS_ExtensionPower(model,k,layers.(stored.name).fullLifts[i].state,coordinates[j]);
                if not KOAHSS_ExtensionStateIsZero(power) then
                    if KOAHSS_ExtensionStateIsZero(canonical) then canonical:=power;
                    else canonical:=model.xtimes(k,canonical,power); fi;
                fi;
            od;
        od;
        Add(canonicals,canonical);
    od;
    stages:=["D","CD","BCD","ABCD"]; attempts:=[]; used:=0; budget:=4096;
    for stageIndex in [1..Length(stages)] do
        stage:=stages[stageIndex];
        for j in [1..Length(choices)] do
            if used>=budget then break; fi;
            remaining:=(Length(stages)-stageIndex+1)*Length(choices)-j+1;
            allowance:=Maximum(1,QuoInt(budget-used,remaining));
            proof:=KOAHSS_ExtensionGaugeCompare(model,k,state,canonicals[j],
                rec(stages:=[stage],maxChoices:=allowance,maxLeadingChoices:=64));
            used:=used+proof.search.differentialEvaluations;
            Add(attempts,rec(stage:=stage,lowerCoordinates:=choices[j],
                status:=proof.status,search:=proof.search,attemptedStages:=proof.attemptedStages));
            if proof.status="computed" then
                return rec(status:="computed",lowerCoordinates:=ShallowCopy(choices[j]),
                    lowerPresentationId:=lower.presentationId,canonicalLowerProduct:=canonicals[j],
                    canonicalComparison:=proof,orderedReductionVerified:=true,
                    coordinateMethod:=coordinateMethod,gaugeSearchAttempts:=attempts,
                    differentialEvaluations:=used,residualState:=model.zero(k));
            fi;
        od;
    od;
    return rec(status:="unresolved",
        reason:="no exact lower-coordinate comparison was found with D, C/D, B/C/D or A/B/C/D gauges within the search bounds",
        coordinateMethod:=coordinateMethod,gaugeSearchAttempts:=attempts,
        differentialEvaluations:=used,residualState:=state);
end);

BindGlobal("KOAHSS_ExtensionHigherOracle",function(backend,k,layers,model)
    local lift,layer,name,i,value,flat,times,flatProduct,boundary,reduce,answer,kernelCache;
    kernelCache:=rec();
    model.gaugeKernelRepresentatives:=function(n,signed,family)
        local key,q,data,generators;
        key:=Concatenation(String(n),"_",String(signed));
        if not IsBound(kernelCache.(key)) then
            q:=-1; if signed then q:=0; fi;
            data:=backend.cohomologyData(n,q);
            generators:=IndependentGeneratorsOfAbelianGroup(data.group);
            kernelCache.(key):=List(generators,g->model.lift(n,data.represent(g),signed));
        fi;
        return kernelCache.(key);
    end;
    lift:=KOAHSS_ExtensionLiftSolver(backend,k,layers,model);
    # Prepare even free generators and generators over a zero lower group.
    # They can be needed in later comparisons, despite requiring no power row.
    for name in ["D","C","B","A"] do
        layer:=layers.(name);
        if not IsBound(layer.fullLifts) then layer.fullLifts:=[]; fi;
        if not IsBound(layer.status) then
            for i in [1..Length(layer.generators)] do
                value:=lift(layer,i);
                if value.status="obstructed" and not IsBound(model.act) then
                    Error("koFull: an E6 survivor has no full flat lift: ",name," ",i," ",value);
                elif value.status<>"computed" then
                    layer.status:="unresolved"; layer.reason:=value.reason; break;
                fi;
            od;
        fi;
    od;
    flat:=function(state,upto)
        if upto<3 then
            if not KOAHSS_ExtensionStateIsZeroThrough(model.d(k,state,upto),upto) then
                Error("koFull: an extension reduction used a nonflat state");
            fi;
        elif not KOAHSS_ExtensionStateIsZero(model.d(k,state)) then
            Error("koFull: an extension reduction used a nonflat state");
        fi;
    end;
    times:=function(x,y,upto)
        if upto<3 then return model.xtimes(k,x,y,upto); fi;
        return model.xtimes(k,x,y);
    end;
    flatProduct:=function(name,coefficients,upto)
        local product,j,power;
        product:=model.zero(k);
        for j in [1..Length(coefficients)] do
            power:=KOAHSS_ExtensionPower(model,k,layers.(name).fullLifts[j].state,coefficients[j],upto);
            if not KOAHSS_ExtensionStateIsZeroThrough(power,upto) then
                if KOAHSS_ExtensionStateIsZeroThrough(product,upto) then product:=power;
                else product:=times(product,power,upto); fi;
            fi;
        od;
        return product;
    end;
    boundary:=function(name,primitive,canonical,upto)
        local gauge,output;
        gauge:=StructuralCopy(model.zero(k-1)); gauge.(name):=primitive;
        if IsBound(model.act) then
            if upto<3 then output:=model.act(k,gauge,canonical,upto);
            else output:=model.act(k,gauge,canonical); fi;
            flat(output,upto);
            return rec(gauge:=gauge,state:=output,canonical:=canonical,
                certificateLevel:="transfer-R",equation:="state = act(gauge, canonical)");
        fi;
        if upto<3 then Error("koFull: a layer-limited reduction needs the transferred model"); fi;
        output:=model.d(k-1,gauge); flat(output,3);
        return rec(gauge:=gauge,state:=output);
    end;
    # Reduce the layers A.. of a stacked power through the layer index upto
    # (D=3 is the complete reduction). Below the target layer of a relation
    # the state is never read, and every operation stays layer-limited.
    reduce:=function(state,lower,upto)
        local current,coefficients,steps,name,n,signed,gens,rows,solution,
            count,values,primitive,gauge,chosen,left,next,stored,j,coordinates,
            canonical,comparison,comparisonSteps,preferred,data,class,projected,fallback,fields;
        fields:=["A","B","C","D"];
        current:=state; coefficients:=rec(A:=[],B:=[],C:=[],D:=[]); steps:=[];
        for name in fields{[1..upto+1]} do
            n:=k+rec(A:=-3,B:=-2,C:=-1,D:=1).(name); signed:=name in ["A","D"];
            stored:=First(lower.layers,l->l.name=name); gens:=[];
            if stored<>fail then gens:=List(layers.(name).fullLifts,x->x.state.(name)); fi;
            count:=Length(gens); rows:=Concatenation(gens,model.matrix(n-1,signed));
            if signed then solution:=koAHSSSolveIntegerSystem(rows,current.(name));
            else solution:=koAHSSSolveMod2System(rows,current.(name)); fi;
            if solution=fail then
                if upto<3 then
                    # The complete measurement handles E6 quotients that
                    # ordinary coboundaries do not generate.
                    return rec(status:="fallback",layer:=name,residualState:=current,reductionSteps:=steps);
                fi;
                preferred:=fail;
                # Once only D remains, its E6 projection determines the
                # marked coordinates without enumerating possible groups.
                # It does not establish the relation: the full comparison
                # below still has to solve and verify the higher gauge.
                if name="D" and IsBound(layers.D.cell) then
                    data:=backend.cohomologyData(n,-4);
                    class:=layers.D.cell.project(data.class(model.project(n,current.D,true)));
                    projected:=KOAHSS_ExtensionGroupCoordinates(layers.D.group,layers.D.generators,class);
                    coefficients.D:=projected;
                    preferred:=List([1..lower.generatorCount],j->0);
                    for stored in lower.layers do
                        for j in [1..Length(stored.orders)] do
                            preferred[stored.startColumn+j-1]:=coefficients.(stored.name)[j];
                        od;
                    od;
                fi;
                fallback:=KOAHSS_ExtensionGaugeReduce(model,k,state,layers,lower,preferred);
                fallback.reductionSteps:=steps;
                fallback.ordinaryReductionFailure:=rec(layer:=name,residualState:=current);
                return fallback;
            fi;
            values:=solution.particular{[1..count]}; coefficients.(name):=values;
            primitive:=solution.particular{[count+1..Length(solution.particular)]};
            chosen:=flatProduct(name,values,upto); gauge:=boundary(name,primitive,chosen,upto);
            if IsBound(model.act) then left:=gauge.state;
            elif KOAHSS_ExtensionStateIsZero(gauge.state) then left:=chosen;
            elif KOAHSS_ExtensionStateIsZero(chosen) then left:=gauge.state;
            else left:=model.xtimes(k,gauge.state,chosen); fi;
            if KOAHSS_ExtensionStateIsZeroThrough(left,upto) then next:=current;
            elif upto<3 then next:=KOAHSS_ExtensionDivideLeft(model,k,left,current,upto);
            else next:=KOAHSS_ExtensionDivideLeft(model,k,left,current); fi;
            if ForAny(next.(name),x->x<>0) then
                Error("koFull: leading component survived its measured reduction");
            fi;
            Add(steps,rec(layer:=name,coordinates:=values,boundary:=gauge,
                chosenProduct:=chosen,before:=current,after:=next,equationVerified:=true));
            current:=next;
        od;
        if not KOAHSS_ExtensionStateIsZeroThrough(current,upto) then
            Error("koFull: nonzero final extension residual");
        fi;
        coordinates:=List([1..lower.generatorCount],j->0);
        canonical:=model.zero(k);
        for stored in lower.layers do
            values:=coefficients.(stored.name);
            for j in [1..Length(values)] do coordinates[stored.startColumn+j-1]:=values[j]; od;
            chosen:=flatProduct(stored.name,values,upto);
            if not KOAHSS_ExtensionStateIsZeroThrough(chosen,upto) then
                if KOAHSS_ExtensionStateIsZeroThrough(canonical,upto) then canonical:=chosen;
                else canonical:=times(canonical,chosen,upto); fi;
            fi;
        od;
        # The preceding equations use an explicit ordered word. Preserve that
        # word and its boundary witnesses; do not assume strict associativity
        # or silently replace it by componentwise addition of the chosen lifts.
        if upto<3 then comparison:=KOAHSS_ExtensionGaugeCompare(model,k,state,canonical,rec(upto:=upto));
        else comparison:=KOAHSS_ExtensionGaugeCompare(model,k,state,canonical); fi;
        if comparison.status<>"computed" then
            return rec(status:="unresolved",reason:="the measured power lacks an exact gauge comparison with the fixed lower product",
                comparison:=comparison,reductionSteps:=steps,residualState:=current);
        fi;
        return rec(status:="computed",lowerCoordinates:=coordinates,
            reductionSteps:=steps,canonicalLowerProduct:=canonical,
            canonicalComparison:=comparison,measuredThrough:=fields[upto+1],
            orderedReductionVerified:=true,residualState:=current);
    end;
    answer:=function(layer,index,order,lower)
        local marked,power,result,fields,target,upto,partial,position,below,lowerNames;
        fields:=["A","B","C","D"];
        marked:=layers.(layer.name);
        if not IsBound(marked.fullLifts) or not IsBound(marked.fullLifts[index])
            or marked.fullLifts[index].status<>"computed" then
            return rec(status:="unresolved",reason:="the fixed full flat representative is unavailable");
        fi;
        position:=Position(fields,layer.name);
        lowerNames:=fields{[position+1..4]};
        # The target-layer shortcut: only the layers down to the lowest free
        # lower generator decide the group (KOAHSS_ExtensionTargetLayer).
        if KOAHSS_LayeredRelationsEnabled() and IsBound(model.layerLimited) and model.layerLimited=true
            and lower.generatorCount>0 then
            target:=KOAHSS_ExtensionTargetLayer(lower,order);
            if target.layer=fail then
                return rec(status:="computed",lowerPresentationId:=lower.presentationId,
                    lowerCoordinates:=List([1..lower.generatorCount],j->0),
                    witness:=rec(operation:="xtimes",power:=order,flatLift:=marked.fullLifts[index],
                        measuredLayers:=[],truncatedBelow:=lowerNames[1],
                        sufficiency:=target,modelId:=model.modelId));
            fi;
            upto:=target.index;
            if upto<3 then
                power:=KOAHSS_ExtensionPower(model,k,marked.fullLifts[index].state,order,upto);
                partial:=reduce(power,lower,upto);
                if partial.status="computed" then
                    below:=fail; if upto+2<=4 then below:=fields[upto+2]; fi;
                    return rec(status:="computed",lowerPresentationId:=lower.presentationId,
                        lowerCoordinates:=partial.lowerCoordinates,witness:=rec(operation:="xtimes",
                            power:=order,flatLift:=marked.fullLifts[index],stackedState:=power,
                            reduction:=partial,measuredLayers:=fields{[position+1..upto+1]},
                            truncatedBelow:=below,sufficiency:=target,modelId:=model.modelId));
                fi;
                # Otherwise the complete measurement below takes over.
            fi;
        fi;
        power:=KOAHSS_ExtensionPower(model,k,marked.fullLifts[index].state,order);
        result:=reduce(power,lower,3);
        if result.status<>"computed" then return result; fi;
        return rec(status:="computed",lowerPresentationId:=lower.presentationId,
            lowerCoordinates:=result.lowerCoordinates,witness:=rec(operation:="xtimes",
                power:=order,flatLift:=marked.fullLifts[index],stackedState:=power,
                reduction:=result,measuredLayers:=lowerNames,truncatedBelow:=fail,
                modelId:=model.modelId));
    end;
    return answer;
end);
