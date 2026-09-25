# Flat representatives for the marked generators used in extension assembly.
# Every vector in a state is expressed in one complete finite cochain model.
# In particular, a chosen B or C is never lifted again from its cohomology
# class after the lower defining equations have been solved.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.

BindGlobal("KOAHSS_ExtensionFlatLift",function(arg)
    local model,degree,leadingName,leading,options,maxChoices,fields,degrees,
        leadingIndex,state,matrices,matrix,checkState,evaluate,zeroCurvature,
        choices,exhausted,solution,lastFailure,search,initial;
    if not Length(arg) in [4,5] then
        Error("koFull: flat lift needs model, degree, layer, cochain[, options]");
    fi;
    model:=arg[1]; degree:=arg[2]; leadingName:=arg[3]; leading:=arg[4];
    options:=rec(); if Length(arg)=5 then options:=arg[5]; fi;
    if not IsRecord(model) or not IsBound(model.dimension)
        or not IsFunction(model.dimension) or not IsBound(model.coboundary)
        or not IsFunction(model.coboundary) or not IsBound(model.d)
        or not IsFunction(model.d) or not IsBound(model.zero)
        or not IsFunction(model.zero) then
        Error("koFull: flat lifts require a complete cochain model and differential");
    fi;
    if not IsInt(degree) or not degree in [-1..6]
        or not IsString(leadingName) or not leadingName in ["A","B","C","D"] then
        Error("koFull: invalid degree or leading layer for a flat lift");
    fi;
    if not IsRecord(options) or ForAny(RecNames(options),n->n<>"maxChoices") then
        Error("koFull: flat lift options accept only maxChoices");
    fi;
    maxChoices:=4096;
    if IsBound(options.maxChoices) then maxChoices:=options.maxChoices; fi;
    if not IsInt(maxChoices) or maxChoices<1 then
        Error("koFull: maxChoices must be a positive integer");
    fi;
    fields:=["A","B","C","D"];
    degrees:=[degree-3,degree-2,degree-1,degree+1];
    leadingIndex:=Position(fields,leadingName);
    if degrees[leadingIndex]<0 then
        Error("koFull: a generator cannot occupy an absent negative-degree layer");
    fi;
    checkState:=function(k,value)
        local j,ns;
        if not IsRecord(value) then Error("koFull: cochain differential must return a state"); fi;
        ns:=[k-3,k-2,k-1,k+1];
        for j in [1..4] do
            if not IsBound(value.(fields[j])) then
                Error("koFull: cochain state is missing a layer");
            fi;
            KOAHSS_CC_CheckVector(value.(fields[j]),model.dimension(ns[j]),
                "complete-model cochain state");
            if j in [2,3] and not ForAll(value.(fields[j]),x->x in [0,1]) then
                Error("koFull: B and C state coordinates must be binary");
            fi;
        od;
    end;
    state:=StructuralCopy(model.zero(degree)); checkState(degree,state);
    if ForAny(fields,f->ForAny(state.(f),x->x<>0)) then
        Error("koFull: the model zero state is not zero");
    fi;
    KOAHSS_CC_CheckVector(leading,model.dimension(degrees[leadingIndex]),
        "flat lift leading cochain");
    if leadingIndex in [2,3] and not ForAll(leading,x->x in [0,1]) then
        Error("koFull: leading B and C cochains must be binary");
    fi;
    state.(leadingName):=ShallowCopy(leading);
    matrices:=rec();
    matrix:=function(n,signed)
        local key,value,j,e;
        key:=Concatenation(String(n),"_",String(signed));
        if IsBound(matrices.(key)) then return matrices.(key); fi;
        if IsBound(model.matrix) then value:=model.matrix(n,signed);
        else
            value:=[];
            for j in [1..model.dimension(n)] do
                e:=List([1..model.dimension(n)],i->0); e[j]:=1;
                Add(value,model.coboundary(n,e,signed));
            od;
        fi;
        if not IsList(value) or Length(value)<>model.dimension(n)
            or not ForAll(value,r->IsList(r) and Length(r)=model.dimension(n+1)
                and ForAll(r,IsInt)) then
            Error("koFull: invalid complete-model coboundary matrix");
        fi;
        matrices.(key):=value; return value;
    end;
    choices:=0; exhausted:=false; solution:=fail; lastFailure:=fail;
    evaluate:=function(value)
        local answer;
        if choices>=maxChoices then exhausted:=true; return fail; fi;
        choices:=choices+1;
        answer:=model.d(degree,value); checkState(degree+1,answer);
        return answer;
    end;
    zeroCurvature:=value->ForAll(fields,f->ForAll(value.(f),x->x=0));
    # Search affine primitive families lazily. Different B choices may be
    # needed before C is soluble, and different C choices may be needed
    # before the integral equation for D is soluble. Keeping the whole
    # affine kernel is exact; a resource limit is never treated as vanishing.
    search:=function(index,value,curvature,witnesses)
        local name,n,rhs,family,j,chosen,child,childCurvature,entry,walk;
        if solution<>fail or exhausted then return; fi;
        if index=5 then
            if not zeroCurvature(curvature) then
                Error("koFull: flat-lift solve left nonzero curvature");
            fi;
            solution:=rec(status:="computed",state:=StructuralCopy(value),
                witness:=rec(leadingLayer:=leadingName,
                    leadingCochain:=ShallowCopy(leading),
                    definingEquations:=StructuralCopy(witnesses),
                    curvature:=StructuralCopy(curvature),flatnessVerified:=true));
            return;
        fi;
        name:=fields[index]; n:=degrees[index];
        rhs:=-curvature.(name);
        if index in [2,3] then
            rhs:=List(rhs,x->x mod 2);
            family:=koAHSSSolveMod2System(matrix(n,false),rhs);
        else family:=koAHSSSolveIntegerSystem(matrix(n,true),rhs); fi;
        if family=fail then
            lastFailure:=rec(layer:=name,rhs:=ShallowCopy(rhs),
                reason:="the current defining system has no primitive");
            return;
        fi;
        chosen:=ShallowCopy(family.particular);
        walk:=function(position)
            local old;
            if solution<>fail or exhausted then return; fi;
            if position<=Length(family.homogeneousGenerators) and index in [2,3] then
                walk(position+1);
                if solution<>fail or exhausted then return; fi;
                old:=chosen;
                chosen:=List(chosen+family.homogeneousGenerators[position],x->x mod 2);
                walk(position+1); chosen:=old;
                return;
            fi;
            child:=StructuralCopy(value); child.(name):=ShallowCopy(chosen);
            # Choosing the already present zero primitive does not change
            # any equation. Reuse its exact audit, rather than spending
            # another nonlinear evaluation (or its resource allowance).
            if child=value then childCurvature:=curvature;
            else childCurvature:=evaluate(child); fi;
            if childCurvature=fail then return; fi;
            if ForAny([1..index],j->ForAny(childCurvature.(fields[j]),x->x<>0)) then
                Error("koFull: the full differential disagrees with its triangular defining equation");
            fi;
            entry:=rec(layer:=name,cochainDegree:=n,rhs:=ShallowCopy(rhs),
                particular:=ShallowCopy(family.particular),
                chosenPrimitive:=ShallowCopy(chosen),
                adjustment:=chosen-family.particular,
                affineKernelRank:=Length(family.homogeneousGenerators));
            if index in [2,3] then
                entry.modulus:=2;
                entry.adjustment:=List(entry.adjustment,x->x mod 2);
            else entry.modulus:=0; fi;
            search(index+1,child,childCurvature,Concatenation(witnesses,[entry]));
        end;
        walk(1);
    end;
    initial:=evaluate(state);
    if ForAny([1..leadingIndex],j->ForAny(initial.(fields[j]),x->x<>0)) then
        return rec(status:="obstructed",reason:="the leading cochain is not closed in its layer",
            failedCurvature:=initial,leadingLayer:=leadingName);
    fi;
    search(leadingIndex+1,state,initial,[]);
    if solution<>fail then
        solution.witness.differentialEvaluations:=choices;
        solution.witness.searchStoppedAtFlatRepresentative:=true;
        if IsBound(model.modelId) then solution.witness.modelId:=model.modelId; fi;
        MakeImmutable(solution);
        return solution;
    fi;
    if exhausted then
        return rec(status:="unresolved",reason:="the flat-lift affine search reached maxChoices",
            differentialEvaluations:=choices,maxChoices:=maxChoices,lastFailure:=lastFailure);
    fi;
    return rec(status:="obstructed",
        reason:="no consistent lower defining cochains satisfy the full differential",
        differentialEvaluations:=choices,lastFailure:=lastFailure);
end);

# The layer basis, complete-model basis, and chosen defining systems remain
# marked together. A later relation query receives these very same objects.
BindGlobal("KOAHSS_ExtensionLiftSolver",function(arg)
    local backend,degree,layers,model,options,cache,answer;
    if not Length(arg) in [4,5] then
        Error("koFull: lift solver needs backend, degree, layers, model[, options]");
    fi;
    backend:=arg[1]; degree:=arg[2]; layers:=arg[3]; model:=arg[4];
    options:=rec(); if Length(arg)=5 then options:=arg[5]; fi;
    if not IsRecord(model) or not IsBound(model.lift) or not IsFunction(model.lift) then
        Error("koFull: the complete model needs a resolution cochain lift");
    fi;
    if IsRecord(layers) then layers:=List(["D","C","B","A"],name->layers.(name)); fi;
    if not IsList(layers) then Error("koFull: lift solver needs the marked layers"); fi;
    cache:=rec();
    answer:=function(layer,index)
        local key,leading,result,marked;
        marked:=First(layers,l->IsIdenticalObj(l,layer));
        if marked=fail or not IsBound(layer.cochains) or not IsInt(index)
            or index<1 or index>Length(layer.cochains) then
            Error("koFull: flat lift requested outside its marked layer basis");
        fi;
        key:=Concatenation(layer.name,"_",String(index));
        if IsBound(cache.(key)) then
            if cache.(key).nativeLeading<>layer.cochains[index] then
                Error("koFull: a cached flat-lift generator changed its cochain basis");
            fi;
            return cache.(key).result;
        fi;
        leading:=model.lift(layer.p,layer.cochains[index],layer.name in ["A","D"]);
        result:=KOAHSS_ExtensionFlatLift(model,degree,layer.name,leading,options);
        cache.(key):=rec(nativeLeading:=ShallowCopy(layer.cochains[index]),result:=result);
        if not IsBound(layer.fullLifts) then layer.fullLifts:=[]; fi;
        layer.fullLifts[index]:=result;
        return result;
    end;
    return answer;
end);
