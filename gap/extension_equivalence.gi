# Exact comparison with a marked, ordered product of flat representatives.
# The verified equation is target = d(gauge) xtimes canonical. No reordering
# or associativity of the cochain product is assumed by this certificate.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.

BindGlobal("KOAHSS_ExtensionGaugeCompare",function(arg)
    local model,k,target,canonical,options,maxChoices,maxLeadingChoices,
        radius,fields,degrees,index,name,n,signed,matrix,rhs,family,synthetic,
        zeroState,zeroBoundary,checkState,isZero,initial,choices,attempts,
        exhausted,result,lastFailure,chosen,coefficients,values,walk,try,
        verify,leadingEquation,search,stages,stageNumber,stage,diagnostics,
        stageChoices,stageAttempts,stageMaxChoices,stageMaxLeading,remaining,
        anyLimited,allExact,stageName,j,stageNames,requestedStages,preferred,
        activeGenerators,kernelSource,seen,independent;
    if not Length(arg) in [4,5] then
        Error("koFull: gauge comparison needs model, degree, target, canonical[, options]");
    fi;
    model:=arg[1]; k:=arg[2]; target:=arg[3]; canonical:=arg[4];
    options:=rec(); if Length(arg)=5 then options:=arg[5]; fi;
    if not IsRecord(model) or not IsBound(model.dimension) or not IsFunction(model.dimension)
        or not IsBound(model.coboundary) or not IsFunction(model.coboundary)
        or not IsBound(model.zero) or not IsFunction(model.zero)
        or not IsBound(model.d) or not IsFunction(model.d)
        or not IsBound(model.xtimes) or not IsFunction(model.xtimes) then
        Error("koFull: gauge comparison needs a complete cochain model with differential and product");
    fi;
    if not IsInt(k) or not k in [0..6] then
        Error("koFull: gauge comparison degree must be in zero through six");
    fi;
    if not IsRecord(options) or ForAny(RecNames(options),x->
        not x in ["maxChoices","maxLeadingChoices","integerRadius","stages"]) then
        Error("koFull: gauge comparison options are maxChoices, maxLeadingChoices, integerRadius, stages");
    fi;
    stageNames:=["D","CD","BCD","ABCD"];
    requestedStages:=stageNames;
    if IsBound(options.stages) then
        if not IsList(options.stages) or Length(options.stages)=0
            or not ForAll(options.stages,x->IsString(x) and x in stageNames) then
            Error("koFull: gauge stages must be a nonempty ordered sublist of D, CD, BCD, ABCD");
        fi;
        requestedStages:=List(options.stages,x->Position(stageNames,x));
        if requestedStages<>Set(requestedStages) then
            Error("koFull: gauge stages must be distinct and ordered D, CD, BCD, ABCD");
        fi;
        requestedStages:=options.stages;
    fi;
    maxChoices:=4096; maxLeadingChoices:=64; radius:=1;
    if IsBound(options.maxChoices) then maxChoices:=options.maxChoices; fi;
    if IsBound(options.maxLeadingChoices) then maxLeadingChoices:=options.maxLeadingChoices; fi;
    if IsBound(options.integerRadius) then radius:=options.integerRadius; fi;
    if not IsInt(maxChoices) or maxChoices<1 or not IsInt(maxLeadingChoices)
        or maxLeadingChoices<1 or not IsInt(radius) or radius<0 then
        Error("koFull: gauge comparison budgets must be positive and integerRadius nonnegative");
    fi;
    fields:=["A","B","C","D"];
    degrees:=[k-4,k-3,k-2,k];
    checkState:=function(degree,state)
        local ns,j;
        if not IsRecord(state) then Error("koFull: gauge comparison needs cochain state records"); fi;
        ns:=[degree-3,degree-2,degree-1,degree+1];
        for j in [1..4] do
            if not IsBound(state.(fields[j])) then
                Error("koFull: gauge comparison state is missing a layer");
            fi;
            KOAHSS_CC_CheckVector(state.(fields[j]),model.dimension(ns[j]),
                "gauge comparison cochain");
            if j in [2,3] and not ForAll(state.(fields[j]),x->x in [0,1]) then
                Error("koFull: gauge comparison B and C cochains must be binary");
            fi;
        od;
    end;
    isZero:=state->ForAll(fields,f->ForAll(state.(f),x->x=0));
    checkState(k,target); checkState(k,canonical);
    initial:=model.d(k,target); checkState(k+1,initial);
    if not isZero(initial) then Error("koFull: gauge comparison target is not flat"); fi;
    initial:=model.d(k,canonical); checkState(k+1,initial);
    if not isZero(initial) then Error("koFull: gauge comparison canonical product is not flat"); fi;
    zeroState:=model.zero(k-1); checkState(k-1,zeroState);
    if not isZero(zeroState) then Error("koFull: gauge model zero is not zero"); fi;
    zeroBoundary:=model.d(k-1,zeroState); checkState(k,zeroBoundary);
    if not isZero(zeroBoundary) then Error("koFull: gauge differential is not pointed at zero"); fi;
    if model.xtimes(k,zeroBoundary,canonical)<>canonical then
        Error("koFull: the zero boundary is not a left unit for the canonical product");
    fi;
    synthetic:=ShallowCopy(model);
    synthetic.d:=function(degree,gauge)
        local boundary,product,difference,j;
        if degree<>k-1 then Error("koFull: gauge equation evaluated in the wrong degree"); fi;
        boundary:=model.d(k-1,gauge); checkState(k,boundary);
        product:=model.xtimes(k,boundary,canonical); checkState(k,product);
        difference:=rec();
        for j in [1..4] do
            difference.(fields[j]):=product.(fields[j])-target.(fields[j]);
            if j in [2,3] then
                difference.(fields[j]):=List(difference.(fields[j]),x->x mod 2);
            fi;
        od;
        return difference;
    end;
    choices:=0; attempts:=0; result:=fail; lastFailure:=fail;
    diagnostics:=[]; anyLimited:=false;
    verify:=function(lift)
        local boundary,product,curvature,answer;
        boundary:=model.d(k-1,lift.state); checkState(k,boundary);
        curvature:=model.d(k,boundary); checkState(k+1,curvature);
        if not isZero(curvature) then
            Error("koFull: gauge boundary failed the exact square-zero identity");
        fi;
        product:=model.xtimes(k,boundary,canonical); checkState(k,product);
        if product<>target then Error("koFull: gauge witness failed its ordered product equality"); fi;
        answer:=rec(status:="computed",gauge:=StructuralCopy(lift.state),
            boundary:=StructuralCopy(boundary),product:=StructuralCopy(product),
            target:=StructuralCopy(target),canonical:=StructuralCopy(canonical),
            equalityVerified:=true,boundaryFlatnessVerified:=true,
            targetFlatnessVerified:=true,canonicalFlatnessVerified:=true,
            equation:="target = d(gauge) xtimes canonical",
            definingSystemWitness:=ShallowCopy(lift.witness));
        # This is a zero of the comparison equation, not a claim that the
        # gauge itself is flat under the original nonlinear differential.
        answer.definingSystemWitness.comparisonEquationVerified:=true;
        Unbind(answer.definingSystemWitness.flatnessVerified);
        if IsBound(model.modelId) then answer.modelId:=model.modelId; fi;
        return answer;
    end;
    # Try the smallest gauge support first. Earlier components are literal
    # zero throughout each stage, so a D-only solution cannot accidentally
    # rely on a B or C representative selected by an earlier calculation.
    stages:=Filtered(List(requestedStages,x->5-Position(stageNames,x)),j->degrees[j]>=0);
    for stageNumber in [1..Length(stages)] do
        index:=stages[stageNumber]; name:=fields[index]; n:=degrees[index];
        signed:=index in [1,4]; stageName:=Concatenation(fields{[index..4]});
        stage:=rec(stage:=stageName,leadingLayer:=name,
            fixedZeroLayers:=fields{[1..index-1]},status:="unresolved",
            leadingChoices:=0,differentialEvaluations:=0);
        Add(diagnostics,stage);
        # Triangularity forces all earlier entries to agree before any
        # cochain in this stage can help. The final certificate uses the
        # complete nonlinear product, rather than this necessary test alone.
        if ForAny([1..index-1],j->target.(fields[j])<>canonical.(fields[j])) then
            stage.reason:="an earlier layer differs while its gauge component is fixed to zero";
            stage.searchComplete:=true;
            continue;
        fi;
        if IsBound(model.matrix) then matrix:=model.matrix(n,signed);
        else
            matrix:=List([1..model.dimension(n)],function(j)
                local e;
                e:=List([1..model.dimension(n)],i->0); e[j]:=1;
                return model.coboundary(n,e,signed);
            end);
        fi;
        rhs:=target.(name)-canonical.(name);
        if signed then family:=koAHSSSolveIntegerSystem(matrix,rhs);
        else
            rhs:=List(rhs,x->x mod 2);
            family:=koAHSSSolveMod2System(matrix,rhs);
        fi;
        if family=fail then
            stage.reason:="the leading difference has no primitive in this gauge stage";
            stage.leadingDifference:=ShallowCopy(rhs); stage.searchComplete:=true;
            continue;
        fi;
        if choices>=maxChoices or attempts>=maxLeadingChoices then
            stage.reason:="the shared gauge comparison budget is exhausted";
            stage.searchComplete:=false; stage.resourceLimit:=true;
            anyLimited:=true; continue;
        fi;
        # Reserve a share for higher support: a large C cocycle kernel must
        # not prevent the B/C/D or A/B/C/D stages from being attempted.
        remaining:=Length(stages)-stageNumber+1;
        stageMaxChoices:=Maximum(1,QuoInt(maxChoices-choices,remaining));
        stageMaxLeading:=Maximum(1,QuoInt(maxLeadingChoices-attempts,remaining));
        if index=4 then stageMaxLeading:=1; fi;
        stageChoices:=0; stageAttempts:=0; exhausted:=false;
        chosen:=ShallowCopy(family.particular);
        coefficients:=List(family.homogeneousGenerators,x->0);
        seen:=[];
        try:=function()
            local lift,used;
            if result<>fail or exhausted then return; fi;
            if chosen in seen then return; fi;
            if stageAttempts>=stageMaxLeading or stageChoices>=stageMaxChoices then
                exhausted:=true; return;
            fi;
            Add(seen,ShallowCopy(chosen));
            attempts:=attempts+1; stageAttempts:=stageAttempts+1;
            lift:=KOAHSS_ExtensionFlatLift(synthetic,k-1,name,chosen,
                rec(maxChoices:=stageMaxChoices-stageChoices));
            if lift.status="computed" then
                used:=lift.witness.differentialEvaluations;
                result:=verify(lift);
                leadingEquation:=rec(layer:=name,cochainDegree:=n,rhs:=ShallowCopy(rhs),
                    particular:=ShallowCopy(family.particular),chosenPrimitive:=ShallowCopy(chosen),
                    kernelCoefficients:=ShallowCopy(coefficients),
                    kernelCoefficientBasis:=StructuralCopy(activeGenerators),
                    kernelSearchSource:=kernelSource,
                    affineKernelRank:=Length(family.homogeneousGenerators));
                if signed then leadingEquation.modulus:=0; else leadingEquation.modulus:=2; fi;
                result.leadingEquation:=leadingEquation;
                result.winningStage:=stageName;
            else
                lastFailure:=lift;
                if IsBound(lift.differentialEvaluations) then used:=lift.differentialEvaluations;
                else used:=1; fi;
                if lift.status="unresolved" then exhausted:=true; fi;
            fi;
            choices:=choices+used; stageChoices:=stageChoices+used;
        end;
        values:=[0,1];
        if signed then
            values:=[0];
            for j in [1..radius] do Add(values,j); Add(values,-j); od;
        fi;
        # Stream affine candidates. There is no later equation after D, so
        # its homogeneous primitives cannot improve the solved comparison.
        walk:=function(position)
            local coefficient,old;
            if result<>fail or exhausted then return; fi;
            if index=4 or position>Length(activeGenerators) then try(); return; fi;
            old:=chosen;
            for coefficient in values do
                coefficients[position]:=coefficient;
                chosen:=old+coefficient*activeGenerators[position];
                if not signed then chosen:=List(chosen,x->x mod 2); fi;
                walk(position+1);
                if result<>fail or exhausted then break; fi;
            od;
            chosen:=old;
        end;
        # Native cohomology lifts may provide useful first candidates, but
        # they do not replace the complete cocycle kernel: coboundary gauge
        # directions can carry nontrivial corrections in lower layers.
        preferred:=[];
        if index<>4 and IsBound(model.gaugeKernelRepresentatives) then
            if not IsFunction(model.gaugeKernelRepresentatives) then
                Error("koFull: gaugeKernelRepresentatives must be a function");
            fi;
            preferred:=model.gaugeKernelRepresentatives(n,signed,family);
            if not IsList(preferred) then Error("koFull: preferred gauge kernel directions must be a list"); fi;
            for chosen in preferred do
                KOAHSS_CC_CheckVector(chosen,model.dimension(n),"preferred gauge kernel direction");
                if not signed and not ForAll(chosen,x->x in [0,1]) then
                    Error("koFull: a preferred binary gauge direction must be binary");
                fi;
                rhs:=model.coboundary(n,chosen,signed);
                if not signed then rhs:=List(rhs,x->x mod 2); fi;
                if ForAny(rhs,x->x<>0) then
                    Error("koFull: a preferred gauge direction is not a cocycle");
                fi;
            od;
            # Remove redundant directions before recursive enumeration;
            # otherwise many zero/dependent hints could consume unbounded
            # traversal time without spending a differential evaluation.
            independent:=[];
            for chosen in preferred do
                if ForAny(chosen,x->x<>0) and not chosen in independent then
                    if signed or koAHSSSolveMod2System(independent,chosen)=fail then
                        Add(independent,ShallowCopy(chosen));
                    fi;
                fi;
            od;
            preferred:=independent;
            rhs:=target.(name)-canonical.(name);
            if not signed then rhs:=List(rhs,x->x mod 2); fi;
        fi;
        if Length(preferred)>0 then
            activeGenerators:=preferred; kernelSource:="preferred";
            chosen:=ShallowCopy(family.particular);
            coefficients:=List(activeGenerators,x->0); walk(1);
        fi;
        activeGenerators:=family.homogeneousGenerators; kernelSource:="complete-kernel";
        chosen:=ShallowCopy(family.particular);
        coefficients:=List(activeGenerators,x->0); walk(1);
        stage.preferredKernelDirections:=Length(preferred);
        stage.leadingChoices:=stageAttempts; stage.differentialEvaluations:=stageChoices;
        stage.maxLeadingChoices:=stageMaxLeading; stage.maxChoices:=stageMaxChoices;
        stage.leadingKernelRank:=Length(family.homogeneousGenerators);
        if result<>fail then
            stage.status:="computed"; stage.searchComplete:=true;
            break;
        fi;
        stage.resourceLimit:=exhausted; stage.searchComplete:=not exhausted;
        # An integral kernel is deliberately bounded; even a completed
        # radius search does not exhaust the integral affine family.
        if signed and index<>4 and Length(family.homogeneousGenerators)>0 then
            stage.searchComplete:=false; stage.integerSearchBounded:=true;
        fi;
        if exhausted then
            stage.reason:="the gauge stage reached its share of the search budget";
            anyLimited:=true;
        else stage.reason:="no exact comparison exists among the searched affine primitives"; fi;
        stage.lastFailure:=lastFailure;
    od;
    search:=rec(leadingChoices:=attempts,differentialEvaluations:=choices,
        maxLeadingChoices:=maxLeadingChoices,maxChoices:=maxChoices,integerRadius:=radius);
    if result<>fail then
        result.search:=search; result.attemptedStages:=diagnostics;
        MakeImmutable(result); return result;
    fi;
    allExact:=ForAll(diagnostics,stage->stage.searchComplete);
    return rec(status:="unresolved",
        reason:="no exact boundary witness was found in the staged affine search",
        equalityVerified:=false,resourceLimit:=anyLimited,search:=search,
        attemptedStages:=diagnostics,searchComplete:=allExact,lastFailure:=lastFailure);
end);
