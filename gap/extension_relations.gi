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

# The primary-operation rows of relations whose target layer lies right below
# the generator's (doc/extension_cup_i_formulas.md) are on unless
# FERMIONAHSS_NATIVE_RELATIONS=0. Tests may set
# KOAHSS_EXTENSION_RELATION_OVERRIDE.native to true or false, which takes
# precedence, and unbind it afterwards.
BindGlobal("KOAHSS_EXTENSION_RELATION_OVERRIDE",rec());
BindGlobal("KOAHSS_NativeRelationsEnabled",function()
    if IsBound(KOAHSS_EXTENSION_RELATION_OVERRIDE.native) then
        return KOAHSS_EXTENSION_RELATION_OVERRIDE.native=true;
    fi;
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_NATIVE_RELATIONS) then
        return GAPInfo.SystemEnvironment.FERMIONAHSS_NATIVE_RELATIONS<>"0";
    fi;
    return true;
end);

# The light rows of the relations that the target-layer shortcut and the
# primary operations do not settle (doc/extensions.md, "Light rows") are on
# unless FERMIONAHSS_LIGHT_RELATIONS=0; KOAHSS_EXTENSION_RELATION_OVERRIDE.light
# takes precedence. FERMIONAHSS_LIGHT_ABSORPTION=0 keeps the pivot B rows of
# the strong absorption regime exact.
BindGlobal("KOAHSS_LightRelationsEnabled",function()
    if IsBound(KOAHSS_EXTENSION_RELATION_OVERRIDE.light) then
        return KOAHSS_EXTENSION_RELATION_OVERRIDE.light=true;
    fi;
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_RELATIONS) then
        return GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_RELATIONS<>"0";
    fi;
    return true;
end);
BindGlobal("KOAHSS_LightAbsorptionEnabled",function()
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_ABSORPTION) then
        return GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_ABSORPTION<>"0";
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
# The optional liftState(name,i) returns the state of the complete lift of a
# lower generator, or fail; by default the lifts stored in layers.X.fullLifts
# are used. A zero coordinate needs no lift.
BindGlobal("KOAHSS_ExtensionGaugeReduce",function(arg)
    local model,k,state,layers,lower,preferred,ranges,stored,i,count,choices,
        canonicals,coordinates,canonical,power,j,stage,stageIndex,stages,
        proof,attempts,used,budget,allowance,remaining,coordinateMethod,liftState,lift;
    if not Length(arg) in [5,6,7] then Error("koFull: invalid gauge reduction arguments"); fi;
    model:=arg[1]; k:=arg[2]; state:=arg[3]; layers:=arg[4]; lower:=arg[5];
    preferred:=fail; if Length(arg)>=6 then preferred:=arg[6]; fi;
    liftState:=function(name,i) return layers.(name).fullLifts[i].state; end;
    if Length(arg)=7 then liftState:=arg[7]; fi;
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
                # A lift of another prime's model (coprime order) never
                # enters a product of this model; its coordinate is irrelevant.
                if IsBound(layers.(stored.name).fullLifts) and IsBound(layers.(stored.name).fullLifts[i])
                    and IsBound(layers.(stored.name).fullLifts[i].model)
                    and layers.(stored.name).fullLifts[i].model<>"complete" then
                    Add(ranges,[0]);
                else Add(ranges,[0..stored.orders[i]-1]); fi;
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
                if coordinates[j]<>0 then
                    lift:=liftState(stored.name,i);
                    if lift=fail then
                        return rec(status:="unresolved",reason:="a lower flat lift of the gauge reduction is unavailable");
                    fi;
                    power:=KOAHSS_ExtensionPower(model,k,lift,coordinates[j]);
                    if not KOAHSS_ExtensionStateIsZero(power) then
                        if KOAHSS_ExtensionStateIsZero(canonical) then canonical:=power;
                        else canonical:=model.xtimes(k,canonical,power); fi;
                    fi;
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

# The relation measurements of one degree k. `source` is a stacking model, or
# a function returning one, called once when the first relation needs a model;
# options.layerLimited states that the model it returns takes layer-limited
# requests. A model and a flat lift are made only for a relation that is
# measured: a split relation, and a relation whose target layer shows that it
# splits, need neither. Flat lifts are solved when first needed and kept on
# the marked layer records (layers.X.fullLifts); a lower generator needs one
# only when it enters a product with a nonzero coefficient.
# engine.answer(layer,index,order,lower) is the relation oracle of
# koAHSSExtensionFromLayers.
BindGlobal("KOAHSS_ExtensionRelationEngine",function(backend,k,layers,source,options)
    local model,built,refusal,kernelCache,liftSolver,local3,layerLimited,liftFailure,
        ensureModel,ensureLocal3,completeLift,localLift,partialLift,partialLifts,markedState,
        liftState,refused,flat,times,flatProduct,boundary,reduce,answer,
        nativeRelations,rowsCache,coboundaryRows,primaryRow,primaryFallbacks,
        context,lightModes,lightFrames,lightFallbacks,aborted,lightMode,startPart,
        ensureFrame,abortLight,lightRow,auditPart,heavyCount;
    model:=fail; built:=fail; refusal:=fail; local3:=fail; liftSolver:=fail; liftFailure:=fail;
    kernelCache:=rec(); partialLifts:=rec(); rowsCache:=rec(); primaryFallbacks:=[];
    context:=fail; if IsBound(options.context) then context:=options.context; fi;
    lightModes:=rec(); lightFrames:=rec(); lightFallbacks:=[]; aborted:=false;
    heavyCount:=0;
    nativeRelations:=KOAHSS_NativeRelationsEnabled();
    if IsBound(options.nativeRelations) then nativeRelations:=options.nativeRelations=true; fi;
    if IsFunction(source) then
        layerLimited:=IsBound(options.layerLimited) and options.layerLimited=true;
    else
        layerLimited:=IsBound(source.layerLimited) and source.layerLimited=true;
    fi;
    ensureModel:=function()
        local candidate;
        if model<>fail then return model; fi;
        if refusal<>fail then return fail; fi;
        if IsFunction(source) then
            candidate:=source(); built:=candidate;
            if not IsRecord(candidate) or not IsBound(candidate.status)
               or candidate.status<>"computed"
               or (IsBound(candidate.supports) and not candidate.supports(k)) then
                refusal:=rec(status:="unresolved",code:="model-setup",
                    reason:="native extension resource budget exceeded in this degree");
                if IsRecord(candidate) and IsBound(candidate.reason) then
                    refusal.reason:=candidate.reason;
                fi;
                return fail;
            fi;
        else candidate:=source; fi;
        model:=candidate;
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
        liftSolver:=KOAHSS_ExtensionLiftSolver(backend,k,layers,model);
        return model;
    end;
    # The two-layer model of the prime three (degrees five and six) shares the
    # worker of the complete model; see doc/extensions.md, "Localization at
    # the primes".
    ensureLocal3:=function()
        if local3<>fail then return local3; fi;
        if ensureModel()=fail then return fail; fi;
        if not IsBound(model.primeLocal) then
            Error("koFull: the stacking model has no three-local model");
        fi;
        local3:=model.primeLocal(3);
        local3.gaugeKernelRepresentatives:=model.gaugeKernelRepresentatives;
        return local3;
    end;
    refused:=function()
        return rec(status:="unresolved",code:="model-setup",reason:=refusal.reason);
    end;
    # The complete flat lift of generator index of layer name, on the marked
    # layer record itself.
    completeLift:=function(name,index)
        local value;
        if ensureModel()=fail then return refused(); fi;
        value:=liftSolver(layers.(name),index);
        if value.status="obstructed" and not IsBound(model.act) then
            Error("koFull: an E6 survivor has no full flat lift: ",name," ",index," ",value);
        fi;
        return value;
    end;
    # The three-local flat lift (A,0,0,D) of an A generator: D solves
    # delta_s D = -J(A), zero in degree five.
    localLift:=function(index)
        local value,leading;
        if IsBound(layers.A.fullLifts) and IsBound(layers.A.fullLifts[index]) then
            return layers.A.fullLifts[index];
        fi;
        if ensureLocal3()=fail then return refused(); fi;
        leading:=local3.lift(layers.A.p,layers.A.cochains[index],true);
        value:=KOAHSS_ExtensionFlatLift(local3,k,"A",leading);
        if value.status="computed" then
            value:=ShallowCopy(value); value.model:="three-local"; value.prime:=3;
            MakeImmutable(value);
        fi;
        if not IsBound(layers.A.fullLifts) then layers.A.fullLifts:=[]; fi;
        layers.A.fullLifts[index]:=value;
        return value;
    end;
    # A lift solved only through the layer index upto (layer-limited
    # differential), for a relation measured through that layer when it is
    # the layer right below the generator's: the measured power does not
    # depend on the choice in that layer, since the generator's order is
    # even and the layer is binary. A complete lift already solved is used
    # as it is.
    partialLift:=function(name,index,upto)
        local key;
        if IsBound(layers.(name).fullLifts) and IsBound(layers.(name).fullLifts[index])
           and layers.(name).fullLifts[index].status="computed" then
            return layers.(name).fullLifts[index];
        fi;
        key:=Concatenation(name,"_",String(index),"_",String(upto));
        if not IsBound(partialLifts.(key)) then
            partialLifts.(key):=KOAHSS_ExtensionFlatLift(model,k,name,
                model.lift(layers.(name).p,layers.(name).cochains[index],name in ["A","D"]),
                rec(upto:=upto));
        fi;
        return partialLifts.(key);
    end;
    # A generator of the last layer that a product is measured through
    # contributes its marked cocycle only: the state with that one layer. It
    # agrees there with its flat lift, whose leading layer is that cocycle,
    # and the products and their curvature through that layer read nothing
    # else.
    markedState:=function(mdl,name,j)
        local state,signed;
        signed:=name in ["A","D"];
        state:=StructuralCopy(mdl.zero(k));
        state.(name):=mdl.lift(layers.(name).p,layers.(name).cochains[j],signed);
        if ForAny(mdl.coboundary(layers.(name).p,state.(name),signed),
                x->(signed and x<>0) or (not signed and x mod 2<>0)) then
            Error("koFull: the marked E6 representative is not closed: ",name," ",j);
        fi;
        return state;
    end;
    # The state of the complete lift of a lower generator, or fail (the
    # reason is kept in liftFailure).
    liftState:=function(name,index)
        local value;
        value:=completeLift(name,index);
        if value.status<>"computed" then
            liftFailure:=Concatenation("the flat lift of generator ",String(index)," of layer ",
                name," is unavailable: ",String(value.reason));
            return fail;
        fi;
        return value.state;
    end;
    flat:=function(mdl,state,upto)
        if upto<3 then
            if not KOAHSS_ExtensionStateIsZeroThrough(mdl.d(k,state,upto),upto) then
                Error("koFull: an extension reduction used a nonflat state");
            fi;
        elif not KOAHSS_ExtensionStateIsZero(mdl.d(k,state)) then
            Error("koFull: an extension reduction used a nonflat state");
        fi;
    end;
    times:=function(mdl,x,y,upto)
        if upto<3 then return mdl.xtimes(k,x,y,upto); fi;
        return mdl.xtimes(k,x,y);
    end;
    # The ordered product of the powers of the complete lifts of layer name
    # with the given coefficients; a zero coefficient needs no lift.
    flatProduct:=function(mdl,name,coefficients,upto)
        local product,j,power,state,last,field;
        product:=mdl.zero(k);
        last:=Position(["A","B","C","D"],name)-1=upto;
        for j in [1..Length(coefficients)] do
            if coefficients[j]<>0 then
                if last then state:=markedState(mdl,name,j);
                else
                    state:=liftState(name,j);
                    if state=fail then return fail; fi;
                    # A layer-limited product reads the layers through upto
                    # only, and returns zeros below them: so does a single
                    # factor, which reaches the comparison unmultiplied.
                    if upto<3 then
                        state:=ShallowCopy(state);
                        for field in ["A","B","C","D"]{[upto+2..4]} do
                            state.(field):=List(state.(field),x->0);
                        od;
                    fi;
                fi;
                power:=KOAHSS_ExtensionPower(mdl,k,state,coefficients[j],upto);
                if not KOAHSS_ExtensionStateIsZeroThrough(power,upto) then
                    if KOAHSS_ExtensionStateIsZeroThrough(product,upto) then product:=power;
                    else product:=times(mdl,product,power,upto); fi;
                fi;
            fi;
        od;
        return product;
    end;
    boundary:=function(mdl,name,primitive,canonical,upto)
        local gauge,output;
        gauge:=StructuralCopy(mdl.zero(k-1)); gauge.(name):=primitive;
        if IsBound(mdl.act) then
            if upto<3 then output:=mdl.act(k,gauge,canonical,upto);
            else output:=mdl.act(k,gauge,canonical); fi;
            flat(mdl,output,upto);
            return rec(gauge:=gauge,state:=output,canonical:=canonical,
                certificateLevel:="transfer-R",equation:="state = act(gauge, canonical)");
        fi;
        if upto<3 then Error("koFull: a layer-limited reduction needs the transferred model"); fi;
        output:=mdl.d(k-1,gauge); flat(mdl,output,3);
        return rec(gauge:=gauge,state:=output);
    end;
    # Reduce the layers A.. of a stacked power through the layer index upto
    # (D=3 is the complete reduction). Below the target layer of a relation
    # the state is never read, and every operation stays layer-limited. In
    # the three-local model the binary layers are absent; its D layer is
    # reduced by the marked D generators.
    reduce:=function(mdl,state,lower,upto)
        local current,coefficients,steps,name,n,signed,gens,rows,solution,
            count,values,primitive,gauge,chosen,left,next,stored,j,coordinates,
            canonical,comparison,preferred,data,class,projected,fallback,fields,
            localMode,active;
        fields:=["A","B","C","D"]; localMode:=IsBound(mdl.localPrime);
        current:=state; coefficients:=rec(A:=[],B:=[],C:=[],D:=[]); steps:=[];
        for name in fields{[1..upto+1]} do
            n:=k+rec(A:=-3,B:=-2,C:=-1,D:=1).(name); signed:=name in ["A","D"];
            stored:=First(lower.layers,l->l.name=name);
            if stored<>fail then coefficients.(name):=List(stored.orders,o->0); fi;
            if localMode and name in ["B","C"] then
                if ForAny(current.(name),x->x<>0) then
                    Error("koFull: a three-local state acquired a binary layer");
                fi;
                Add(steps,rec(layer:=name,coordinates:=coefficients.(name),absentLayer:=true));
                continue;
            fi;
            # The rows are the marked leading cochains: the leading layer of a
            # flat lift is the marked cochain itself (KOAHSS_ExtensionLiftSolver).
            gens:=[]; active:=[];
            if stored<>fail then
                for j in [1..Length(stored.orders)] do
                    Add(gens,mdl.lift(layers.(name).p,layers.(name).cochains[j],signed));
                    Add(active,j);
                od;
            fi;
            count:=Length(gens); rows:=Concatenation(gens,mdl.matrix(n-1,signed));
            if signed then solution:=koAHSSSolveIntegerSystem(rows,current.(name));
            else solution:=koAHSSSolveMod2System(rows,current.(name)); fi;
            if solution=fail then
                if upto<3 or localMode then
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
                    class:=layers.D.cell.project(data.class(mdl.project(n,current.D,true)));
                    projected:=KOAHSS_ExtensionGroupCoordinates(layers.D.group,layers.D.generators,class);
                    coefficients.D:=projected;
                    preferred:=List([1..lower.generatorCount],j->0);
                    for stored in lower.layers do
                        for j in [1..Length(stored.orders)] do
                            preferred[stored.startColumn+j-1]:=coefficients.(stored.name)[j];
                        od;
                    od;
                fi;
                fallback:=KOAHSS_ExtensionGaugeReduce(mdl,k,state,layers,lower,preferred,
                    function(name,i)
                        if name="D" then return markedState(mdl,name,i); fi;
                        return liftState(name,i);
                    end);
                fallback.reductionSteps:=steps;
                fallback.ordinaryReductionFailure:=rec(layer:=name,residualState:=current);
                return fallback;
            fi;
            values:=solution.particular{[1..count]};
            for j in [1..count] do coefficients.(name)[active[j]]:=values[j]; od;
            primitive:=solution.particular{[count+1..Length(solution.particular)]};
            chosen:=flatProduct(mdl,name,coefficients.(name),upto);
            if chosen=fail then
                return rec(status:="unresolved",reason:=liftFailure,residualState:=current,reductionSteps:=steps);
            fi;
            gauge:=boundary(mdl,name,primitive,chosen,upto);
            if IsBound(mdl.act) then left:=gauge.state;
            elif KOAHSS_ExtensionStateIsZero(gauge.state) then left:=chosen;
            elif KOAHSS_ExtensionStateIsZero(chosen) then left:=gauge.state;
            else left:=mdl.xtimes(k,gauge.state,chosen); fi;
            if KOAHSS_ExtensionStateIsZeroThrough(left,upto) then next:=current;
            elif upto<3 then next:=KOAHSS_ExtensionDivideLeft(mdl,k,left,current,upto);
            else next:=KOAHSS_ExtensionDivideLeft(mdl,k,left,current); fi;
            if ForAny(next.(name),x->x<>0) then
                Error("koFull: leading component survived its measured reduction");
            fi;
            Add(steps,rec(layer:=name,coordinates:=coefficients.(name),boundary:=gauge,
                chosenProduct:=chosen,before:=current,after:=next,equationVerified:=true));
            current:=next;
        od;
        if not KOAHSS_ExtensionStateIsZeroThrough(current,upto) then
            Error("koFull: nonzero final extension residual");
        fi;
        coordinates:=List([1..lower.generatorCount],j->0);
        canonical:=mdl.zero(k);
        for stored in lower.layers do
            values:=coefficients.(stored.name);
            for j in [1..Length(values)] do coordinates[stored.startColumn+j-1]:=values[j]; od;
            chosen:=flatProduct(mdl,stored.name,values,upto);
            if chosen=fail then
                return rec(status:="unresolved",reason:=liftFailure,residualState:=current,reductionSteps:=steps);
            fi;
            if not KOAHSS_ExtensionStateIsZeroThrough(chosen,upto) then
                if KOAHSS_ExtensionStateIsZeroThrough(canonical,upto) then canonical:=chosen;
                else canonical:=times(mdl,canonical,chosen,upto); fi;
            fi;
        od;
        # The preceding equations use an explicit ordered word. Preserve that
        # word and its boundary witnesses; do not assume strict associativity
        # or silently replace it by componentwise addition of the chosen lifts.
        if upto<3 then comparison:=KOAHSS_ExtensionGaugeCompare(mdl,k,state,canonical,rec(upto:=upto));
        else comparison:=KOAHSS_ExtensionGaugeCompare(mdl,k,state,canonical); fi;
        if comparison.status<>"computed" then
            return rec(status:="unresolved",reason:="the measured power lacks an exact gauge comparison with the fixed lower product",
                comparison:=comparison,reductionSteps:=steps,residualState:=current);
        fi;
        return rec(status:="computed",lowerCoordinates:=coordinates,
            reductionSteps:=steps,canonicalLowerProduct:=canonical,
            canonicalComparison:=comparison,measuredThrough:=fields[upto+1],
            orderedReductionVerified:=true,residualState:=current);
    end;
    # The coboundary matrix of R in degree n as rows, from the backend.
    coboundaryRows:=function(n,signed)
        local key,d;
        key:=Concatenation(String(n),"_",String(signed));
        if not IsBound(rowsCache.(key)) then
            d:=backend.dimension(n);
            rowsCache.(key):=List([1..d],function(j)
                local e;
                e:=List([1..d],i->0); e[j]:=1;
                return backend.coboundary(n,e,signed);
            end);
        fi;
        return rowsCache.(key);
    end;
    # The row of a relation whose target layer lies right below the
    # generator's, read from the class of a primary operation of the marked
    # cocycle on R (doc/extension_cup_i_formulas.md): D(rho u) with
    # delta_s u = m A for A over B, the reduced twisted Bockstein for B over C,
    # and an integral lift of D(c) for C over D. No stacking model, flat lift
    # or bar transport is used. The target-layer coordinates come from the E6
    # cell projection, which quotients by the incoming images. A C-over-D row
    # is determined modulo 2H, the relation of a lift shifted within D, and
    # its witness says so (shiftedLift). Returns fail when an ingredient is
    # unavailable; the relation is then measured in the model.
    primaryRow:=function(layer,index,order,lower,target)
        local fields,name,solution,z,value,twice,degree,q,element,coordinates,
            stored,row,j,witness,method;
        fields:=["A","B","C","D"]; name:=target.layer;
        if not IsBound(backend.nativePrimary) or not IsBound(backend.hasCupMod2)
           or backend.hasCupMod2<>true or not IsBound(layers.(name).cell)
           or not IsBound(layer.cochains) or not IsBound(layer.cochains[index]) then
            return fail;
        fi;
        if layer.name="A" then
            if k<4 or not IsEvenInt(order) then return fail; fi;
            solution:=koAHSSSolveIntegerSystem(coboundaryRows(k-4,true),
                order*layer.cochains[index]);
            if solution=fail then return fail; fi;
            z:=List(solution.particular,x->x mod 2);
            value:=backend.nativePrimary("D",k-4,z); degree:=k-2; q:=-1;
            method:="D(rho u), delta_s u = m A";
        elif layer.name="B" then
            twice:=backend.coboundary(k-2,List(layer.cochains[index],x->x mod 2),true);
            if ForAny(twice,IsOddInt) then return fail; fi;
            value:=List(twice,x->(x/2) mod 2); degree:=k-1; q:=-2;
            method:="rho beta_s(b)";
        else
            value:=backend.nativePrimary("D",k-1,List(layer.cochains[index],x->x mod 2));
            twice:=backend.coboundary(k+1,value,true);
            if ForAny(twice,IsOddInt) then return fail; fi;
            solution:=koAHSSSolveIntegerSystem(coboundaryRows(k+1,true),List(twice,x->x/2));
            if solution=fail then return fail; fi;
            value:=value-2*solution.particular; degree:=k+1; q:=-4;
            if ForAny(backend.coboundary(k+1,value,true),x->x<>0) then return fail; fi;
            method:="integral lift of D(c)";
        fi;
        element:=layers.(name).cell.project(backend.cohomologyData(degree,q).class(value));
        coordinates:=KOAHSS_ExtensionGroupCoordinates(layers.(name).group,
            layers.(name).generators,element);
        stored:=First(lower.layers,l->l.name=name);
        if stored=fail or Length(stored.orders)<>Length(coordinates) then return fail; fi;
        row:=List([1..lower.generatorCount],j->0);
        for j in [1..Length(coordinates)] do row[stored.startColumn+j-1]:=coordinates[j]; od;
        witness:=rec(operation:="xtimes",power:=order,prime:=2,model:="primary-R",
            method:=method,measuredLayers:=[name],sufficiency:=target);
        if name="D" then
            witness.truncatedBelow:=fail;
            witness.shiftedLift:="D";
        else
            witness.truncatedBelow:=fields[Position(fields,name)+1];
        fi;
        return rec(status:="computed",lowerPresentationId:=lower.presentationId,
            lowerCoordinates:=row,witness:=witness);
    end;
    # Light rows (gap/extension_light.gi). KOAHSS_ExtensionPrimeRows starts
    # every prime part; the two-primary part, and the three-primary part in
    # degrees five and six, are light when the switches and the page context
    # allow it. A light part never reaches the model measurement: a light row
    # that fails for a reason other than a resource limit aborts the part,
    # which is recorded in lightFallbacks and computed again in the model.
    lightMode:=prime->prime<>fail and IsBound(lightModes.(String(prime)))
        and lightModes.(String(prime))="light";
    startPart:=function(part)
        local allowed;
        aborted:=false;
        allowed:=KOAHSS_LightRelationsEnabled() and context<>fail and layerLimited and nativeRelations
            and KOAHSS_LayeredRelationsEnabled() and IsFunction(source)
            and IsBound(backend.nativePrimary) and IsBound(backend.hasCupMod2) and backend.hasCupMod2=true
            and IsBound(context.getCell) and IsBound(context.getMap)
            and (part.prime=2 or (part.prime=3 and k in [5,6] and KOAHSS_PrimeLocalizationEnabled()));
        if allowed then
            lightModes.(String(part.prime)):="light";
            lightFrames.(String(part.prime)):=rec(part:=part,frame:=fail);
        else
            lightModes.(String(part.prime)):="heavy";
        fi;
    end;
    ensureFrame:=function(prime)
        local entry;
        entry:=lightFrames.(String(prime));
        if entry.frame=fail then
            entry.frame:=CallFuncList(ValueGlobal("KOAHSS_LightFrame"),[rec(backend:=backend,k:=k,
                layers:=layers,context:=context,part:=entry.part,coboundaryRows:=coboundaryRows,
                call:=function(task,data)
                    if ensureModel()=fail then return refused(); fi;
                    if not IsBound(model.light) then
                        return rec(status:="unresolved",code:="model-setup",
                            reason:="the stacking model has no light evaluator");
                    fi;
                    return model.light(task,data);
                end)]);
        fi;
        return entry.frame;
    end;
    abortLight:=function(prime,code,reason,name,index)
        Add(lightFallbacks,rec(prime:=prime,layer:=name,index:=index,code:=code,reason:=reason));
        lightModes.(String(prime)):="heavy"; Unbind(lightFrames.(String(prime)));
        aborted:=true;
        # A heavy measurement restarts the worker with a fresh setup.
        if built<>fail and IsRecord(built) and IsBound(built.close) then built.close(); fi;
        return rec(status:="unresolved",code:="light-aborted",reason:=Concatenation(
            "the light rows of the prime ",String(prime)," are computed again in the model: ",reason));
    end;
    lightRow:=function(layer,index,order,lower,target,options)
        local prime,frame,before,attempt,oldBreak,oldSilent,refusalAnswer;
        prime:=KOAHSS_RelationPrime(order);
        frame:=ensureFrame(prime);
        before:=fail;
        if model<>fail and IsBound(model.lastFailure) then before:=model.lastFailure; fi;
        if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_DEBUG)
           and GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_DEBUG="raise" then
            return frame.row(layer.name,index,order,lower,target,options);
        fi;
        oldBreak:=BreakOnError; oldSilent:=SilentNonInteractiveErrors;
        BreakOnError:=false; SilentNonInteractiveErrors:=true;
        attempt:=CALL_WITH_CATCH(frame.row,[layer.name,index,order,lower,target,options]);
        BreakOnError:=oldBreak; SilentNonInteractiveErrors:=oldSilent;
        if attempt[1] then return attempt[2]; fi;
        if frame.refusal<>fail then
            refusalAnswer:=frame.refusal; frame.refusal:=fail;
            if IsBound(refusalAnswer.code) and refusalAnswer.code="model-setup" then return refusalAnswer; fi;
            if built<>fail and IsRecord(built) and IsBound(built.close) then built.close(); fi;
            return rec(status:="unresolved",code:="resource-limit",
                reason:=Concatenation("light row: ",String(refusalAnswer.reason)));
        fi;
        if model<>fail and IsBound(model.lastFailure) and not IsIdenticalObj(model.lastFailure,before)
           and IsBound(model.lastFailure.status) and model.lastFailure.status="unresolved" then
            return rec(status:="unresolved",code:="resource-limit",reason:=model.lastFailure.reason);
        fi;
        return abortLight(prime,"light-error",frame.lastStep,layer.name,index);
    end;
    auditPart:=function(part,rows,lowerOf)
        local entry,attempt,oldBreak,oldSilent,reason;
        if not lightMode(part.prime) then return; fi;
        entry:=lightFrames.(String(part.prime));
        if entry.frame=fail then return; fi;
        oldBreak:=BreakOnError; oldSilent:=SilentNonInteractiveErrors;
        BreakOnError:=false; SilentNonInteractiveErrors:=true;
        attempt:=CALL_WITH_CATCH(entry.frame.audit,[rows,lowerOf]);
        BreakOnError:=oldBreak; SilentNonInteractiveErrors:=oldSilent;
        if attempt[1] and attempt[2]=fail then return; fi;
        reason:="the final frame audit raised an error";
        if attempt[1] then reason:=attempt[2]; fi;
        abortLight(part.prime,"frame-audit",reason,fail,fail);
    end;
    # options.complete measures the relation through D, without the
    # target-layer shortcut (KOAHSS_ExtensionPrimeRows asks for it when a
    # later relation depends on the D components this relation dropped).
    answer:=function(arg)
        local layer,index,order,lower,complete,power,result,fields,target,upto,partial,position,
            below,lowerNames,kind,entry,witness,reason,oldBreak,oldSilent,attempt;
        layer:=arg[1]; index:=arg[2]; order:=arg[3]; lower:=arg[4];
        complete:=Length(arg)>=5 and IsBound(arg[5].complete) and arg[5].complete=true;
        fields:=["A","B","C","D"];
        kind:=KOAHSS_ExtensionRelationModel(k,layer.name,order);
        if kind.model="split" then return KOAHSS_ExtensionSplitResponse(order,kind.prime,lower); fi;
        if complete and lightMode(kind.prime) then
            return lightRow(layer,index,order,lower,fail,arg[5]);
        fi;
        position:=Position(fields,layer.name);
        lowerNames:=fields{[position+1..4]};
        # The target-layer shortcut: only the layers down to the lowest free
        # lower generator decide the group (KOAHSS_ExtensionTargetLayer). It
        # reads the lower presentation only, so a relation with no free lower
        # generator needs neither a model nor a lift.
        target:=fail;
        if not complete and KOAHSS_LayeredRelationsEnabled() and layerLimited and lower.generatorCount>0 then
            target:=KOAHSS_ExtensionTargetLayer(lower,order);
            if target.layer=fail then
                witness:=rec(operation:="xtimes",power:=order,prime:=kind.prime,model:=kind.model,
                    measuredLayers:=[],truncatedBelow:=lowerNames[1],sufficiency:=target);
                if model<>fail then witness.modelId:=model.modelId; fi;
                return rec(status:="computed",lowerPresentationId:=lower.presentationId,
                    lowerCoordinates:=List([1..lower.generatorCount],j->0),witness:=witness);
            fi;
            # A target layer right below the generator's: the primary-operation
            # row. An error in its evaluation is caught without a message, and
            # the relation is measured in the model and listed in primaryFallbacks.
            if nativeRelations and kind.model="complete" and IsEvenInt(order)
               and target.index=position then
                oldBreak:=BreakOnError; oldSilent:=SilentNonInteractiveErrors;
                BreakOnError:=false; SilentNonInteractiveErrors:=true;
                attempt:=CALL_WITH_CATCH(primaryRow,[layer,index,order,lower,target]);
                BreakOnError:=oldBreak; SilentNonInteractiveErrors:=oldSilent;
                if attempt[1] and attempt[2]<>fail then return attempt[2]; fi;
                if lightMode(kind.prime) then
                    if layer.name="C" then return lightRow(layer,index,order,lower,target,rec()); fi;
                    return abortLight(kind.prime,"primary-unavailable",
                        "a primary-operation row is unavailable",layer.name,index);
                fi;
                Add(primaryFallbacks,rec(layer:=layer.name,index:=index,
                    reason:=function() if attempt[1] then return "an ingredient was unavailable"; fi;
                        return "the evaluation raised an error"; end()));
            fi;
            if lightMode(kind.prime) and kind.model in ["complete","three-local"]
               and target.index>position then
                return lightRow(layer,index,order,lower,target,rec());
            fi;
        fi;
        if lightMode(kind.prime) then
            return abortLight(kind.prime,"heavy-entry","a relation would be measured in the model",
                layer.name,index);
        fi;
        # Count entry into the old measurement path, including unsuccessful
        # attempts. A worker used only to pair light residues does not count.
        heavyCount:=heavyCount+1;
        if kind.model="three-local" then
            # The two-layer measurement of the prime three: the A layer is
            # reduced by coboundaries, the D layer by the marked D generators;
            # the binary layers are absent.
            if ensureLocal3()=fail then return refused(); fi;
            entry:=localLift(index);
            if entry.status<>"computed" then
                return rec(status:="unresolved",reason:=Concatenation("three-local flat lift: ",String(entry.reason)));
            fi;
            power:=KOAHSS_ExtensionPower(local3,k,entry.state,order);
            partial:=reduce(local3,power,lower,3);
            if partial.status="computed" then
                return rec(status:="computed",lowerPresentationId:=lower.presentationId,
                    lowerCoordinates:=partial.lowerCoordinates,witness:=rec(operation:="xtimes",
                        power:=order,prime:=kind.prime,model:="three-local",flatLift:=entry,
                        stackedState:=power,reduction:=partial,
                        measuredLayers:=Filtered(lowerNames,n->n in ["A","D"]),
                        binaryLayersAbsent:=true,truncatedBelow:=fail,sufficiency:=target,
                        modelId:=local3.modelId));
            fi;
            # The complete four-layer measurement would give the relation in
            # the lifts of the binary generators; its three-local row needs
            # their complete relations, which are not measured here.
            if IsBound(partial.reason) then reason:=partial.reason;
            else reason:=Concatenation("the three-local reduction does not close in layer ",partial.layer); fi;
            return rec(status:="unresolved",code:="three-local-reduction",
                reason:=Concatenation("three-local relation: ",reason),threeLocalReduction:=partial);
        fi;
        if ensureModel()=fail then return refused(); fi;
        if target<>fail then
            upto:=target.index;
            if upto<3 then
                # Through the layer right below the generator's, a lift solved
                # through that layer suffices; a lower target needs a complete
                # lift, whose choices in the layers between are those of an
                # element of the group.
                if upto=position then entry:=partialLift(layer.name,index,upto);
                else entry:=completeLift(layer.name,index); fi;
                if entry.status<>"computed" then
                    return rec(status:="unresolved",reason:=entry.reason);
                fi;
                power:=KOAHSS_ExtensionPower(model,k,entry.state,order,upto);
                partial:=reduce(model,power,lower,upto);
                if partial.status="computed" then
                    below:=fail; if upto+2<=4 then below:=fields[upto+2]; fi;
                    return rec(status:="computed",lowerPresentationId:=lower.presentationId,
                        lowerCoordinates:=partial.lowerCoordinates,witness:=rec(operation:="xtimes",
                            power:=order,prime:=kind.prime,model:="complete",flatLift:=entry,stackedState:=power,
                            reduction:=partial,measuredLayers:=fields{[position+1..upto+1]},
                            truncatedBelow:=below,sufficiency:=target,modelId:=model.modelId));
                fi;
                # Otherwise the complete measurement below takes over.
            fi;
        fi;
        entry:=completeLift(layer.name,index);
        if entry.status<>"computed" then
            return rec(status:="unresolved",reason:=entry.reason);
        fi;
        power:=KOAHSS_ExtensionPower(model,k,entry.state,order);
        result:=reduce(model,power,lower,3);
        if result.status<>"computed" then return result; fi;
        return rec(status:="computed",lowerPresentationId:=lower.presentationId,
            lowerCoordinates:=result.lowerCoordinates,witness:=rec(operation:="xtimes",
                power:=order,prime:=kind.prime,model:="complete",flatLift:=entry,stackedState:=power,
                reduction:=result,measuredLayers:=lowerNames,truncatedBelow:=fail,
                modelId:=model.modelId));
    end;
    if not IsFunction(source) then ensureModel(); fi;
    return rec(answer:=answer,
        # the model in use, or fail when no relation has needed one
        model:=function() return model; end,
        # relations whose primary-operation row fell back to the model
        primaryFallbacks:=function() return primaryFallbacks; end,
        # light rows: the part decision, the abort flag of the current part,
        # and the light parts computed again in the model
        startPart:=startPart,
        auditPart:=auditPart,
        heavyMeasurements:=function() return heavyCount; end,
        lightAborted:=function() return aborted; end,
        resetPart:=function(prime) aborted:=false; end,
        lightFallbacks:=function() return lightFallbacks; end,
        # what the model source returned (a refusal included), or fail
        builtModel:=function() return built; end,
        close:=function()
            if built<>fail and IsRecord(built) and IsBound(built.close) then built.close(); fi;
        end);
end);

# The relation oracle of a fixed model (koAHSSExtensionFromLayers).
BindGlobal("KOAHSS_ExtensionHigherOracle",function(backend,k,layers,model)
    return KOAHSS_ExtensionRelationEngine(backend,k,layers,model,rec()).answer;
end);
