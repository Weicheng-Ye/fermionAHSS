# Measure powers of fixed flat representatives, retaining every boundary carry.
# Every division below verifies an actual ordered cochain equation. No lower
# coordinate is identified merely from the order of a stacked element.
BindGlobal("KOAHSS_ExtensionStateIsZero",state->ForAll(["A","B","C","D"],
    name->ForAll(state.(name),x->x=0)));

BindGlobal("KOAHSS_ExtensionDivideLeft",function(model,k,left,total)
    local right,product,name;
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

BindGlobal("KOAHSS_ExtensionPower",function(model,k,state,m)
    local answer,base;
    if not IsInt(m) then Error("koFull: a stacking power must be integral"); fi;
    if m<0 then
        state:=KOAHSS_ExtensionDivideLeft(model,k,state,model.zero(k)); m:=-m;
    fi;
    answer:=model.zero(k); base:=state;
    while m>0 do
        if m mod 2=1 then
            if KOAHSS_ExtensionStateIsZero(answer) then answer:=base;
            else answer:=model.xtimes(k,answer,base); fi;
        fi;
        m:=QuoInt(m,2);
        if m>0 then base:=model.xtimes(k,base,base); fi;
    od;
    if not KOAHSS_ExtensionStateIsZero(model.d(k,answer)) then
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
    local lift,layer,name,i,value,flat,flatProduct,boundary,reduce,answer,kernelCache;
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
                if value.status="obstructed" then
                    Error("koFull: an E6 survivor has no full flat lift: ",name," ",i," ",value);
                elif value.status<>"computed" then
                    layer.status:="unresolved"; layer.reason:=value.reason; break;
                fi;
            od;
        fi;
    od;
    flat:=function(state)
        if not KOAHSS_ExtensionStateIsZero(model.d(k,state)) then
            Error("koFull: an extension reduction used a nonflat state");
        fi;
    end;
    flatProduct:=function(name,coefficients)
        local product,j,power;
        product:=model.zero(k);
        for j in [1..Length(coefficients)] do
            power:=KOAHSS_ExtensionPower(model,k,layers.(name).fullLifts[j].state,coefficients[j]);
            if not KOAHSS_ExtensionStateIsZero(power) then
                if KOAHSS_ExtensionStateIsZero(product) then product:=power;
                else product:=model.xtimes(k,product,power); fi;
            fi;
        od;
        return product;
    end;
    boundary:=function(name,primitive)
        local gauge,output;
        gauge:=StructuralCopy(model.zero(k-1)); gauge.(name):=primitive;
        output:=model.d(k-1,gauge); flat(output);
        return rec(gauge:=gauge,state:=output);
    end;
    reduce:=function(state,lower)
        local current,coefficients,steps,name,n,signed,gens,rows,solution,
            count,values,primitive,gauge,chosen,left,next,stored,j,coordinates,
            canonical,comparison,comparisonSteps,preferred,data,class,projected,fallback;
        current:=state; coefficients:=rec(A:=[],B:=[],C:=[],D:=[]); steps:=[];
        for name in ["A","B","C","D"] do
            n:=k+rec(A:=-3,B:=-2,C:=-1,D:=1).(name); signed:=name in ["A","D"];
            stored:=First(lower.layers,l->l.name=name); gens:=[];
            if stored<>fail then gens:=List(layers.(name).fullLifts,x->x.state.(name)); fi;
            count:=Length(gens); rows:=Concatenation(gens,model.matrix(n-1,signed));
            if signed then solution:=koAHSSSolveIntegerSystem(rows,current.(name));
            else solution:=koAHSSSolveMod2System(rows,current.(name)); fi;
            if solution=fail then
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
            gauge:=boundary(name,primitive); chosen:=flatProduct(name,values);
            if KOAHSS_ExtensionStateIsZero(gauge.state) then left:=chosen;
            elif KOAHSS_ExtensionStateIsZero(chosen) then left:=gauge.state;
            else left:=model.xtimes(k,gauge.state,chosen); fi;
            if KOAHSS_ExtensionStateIsZero(left) then next:=current;
            else next:=KOAHSS_ExtensionDivideLeft(model,k,left,current); fi;
            if ForAny(next.(name),x->x<>0) then
                Error("koFull: leading component survived its measured reduction");
            fi;
            Add(steps,rec(layer:=name,coordinates:=values,boundary:=gauge,
                chosenProduct:=chosen,before:=current,after:=next,equationVerified:=true));
            current:=next;
        od;
        if not KOAHSS_ExtensionStateIsZero(current) then Error("koFull: nonzero final extension residual"); fi;
        coordinates:=List([1..lower.generatorCount],j->0);
        canonical:=model.zero(k);
        for stored in lower.layers do
            values:=coefficients.(stored.name);
            for j in [1..Length(values)] do coordinates[stored.startColumn+j-1]:=values[j]; od;
            chosen:=flatProduct(stored.name,values);
            if not KOAHSS_ExtensionStateIsZero(chosen) then
                if KOAHSS_ExtensionStateIsZero(canonical) then canonical:=chosen;
                else canonical:=model.xtimes(k,canonical,chosen); fi;
            fi;
        od;
        # The preceding equations use an explicit ordered word. Preserve that
        # word and its boundary witnesses; do not assume strict associativity
        # or silently replace it by componentwise addition of the chosen lifts.
        comparison:=KOAHSS_ExtensionGaugeCompare(model,k,state,canonical);
        if comparison.status<>"computed" then
            return rec(status:="unresolved",reason:="the measured power lacks an exact gauge comparison with the fixed lower product",
                comparison:=comparison,reductionSteps:=steps,residualState:=current);
        fi;
        return rec(status:="computed",lowerCoordinates:=coordinates,
            reductionSteps:=steps,canonicalLowerProduct:=canonical,
            canonicalComparison:=comparison,
            orderedReductionVerified:=true,residualState:=current);
    end;
    answer:=function(layer,index,order,lower)
        local marked,power,result;
        marked:=layers.(layer.name);
        if not IsBound(marked.fullLifts) or not IsBound(marked.fullLifts[index])
            or marked.fullLifts[index].status<>"computed" then
            return rec(status:="unresolved",reason:="the fixed full flat representative is unavailable");
        fi;
        power:=KOAHSS_ExtensionPower(model,k,marked.fullLifts[index].state,order);
        result:=reduce(power,lower);
        if result.status<>"computed" then return result; fi;
        return rec(status:="computed",lowerPresentationId:=lower.presentationId,
            lowerCoordinates:=result.lowerCoordinates,witness:=rec(operation:="xtimes",
                power:=order,flatLift:=marked.fullLifts[index],stackedState:=power,
                reduction:=result,modelId:=model.modelId));
    end;
    return answer;
end);
