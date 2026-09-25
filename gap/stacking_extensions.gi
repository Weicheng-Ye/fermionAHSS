# Low-degree extension relations from the selected four-cochain xtimes.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.
#
# The degree-two local phase is the exact finite restriction of
# fermionAHSS_stacking/coherent_low_commutative.py to closed B0,C1, omega=0.
# Its 256 triangle values and source hashes are bundled in the data file.
# The runtime has no dependency on that sibling workspace. This module never
# infers an extension from an order or from an expected classification table.
CallFuncList(function()
    local name, slashes, directory, payload, encoded, digest;
    name := INPUT_FILENAME();
    if not IsEmpty(name) and name[1]<>'/' then
        name := Filename(DirectoryCurrent(),name);
    fi;
    slashes := Positions(name,'/');
    directory := Directory(name{[1..Last(slashes)]});
    payload := CallFuncList(ValueGlobal("JsonStringToGap"),
        [StringFile(Filename(directory,"../data/stacking-low-phase.json"))]);
    if not IsRecord(payload) or not IsBound(payload.schema) or payload.schema<>1
       or not IsBound(payload.calibration)
       or payload.calibration<>"coherent-low-commutative/closed-B-C/omega-zero/v1"
       or not IsBound(payload.values) or Length(payload.values)<>256
       or not ForAll(payload.values,v->IsList(v) and Length(v)=2
            and ForAll(v,IsInt) and v[2]>0) then
        Error("koFull: malformed low-degree stacking phase data");
    fi;
    encoded := Concatenation("[",JoinStringsWithSeparator(List(payload.values,
        v->Concatenation("[",String(v[1]),",",String(v[2]),"]")),","),"]");
    digest := "4ae828fdfa4c5ebcc93b3919c4f7ee9615d356f432dda908256130f1b456b11b";
    if not IsBound(payload.values_sha256) or payload.values_sha256<>digest
       or HexSHA256(encoded)<>digest then
        Error("koFull: low-degree stacking calibration checksum mismatch");
    fi;
    MakeImmutable(payload);
    BindGlobal("KOAHSS_STACKING_LOW_PHASE",payload);
end,[]);

# All arguments here are binary cochain evaluators in local bar coordinates.
# B,Bprime are constant zero-cochains and C,Cprime,s are one-cocycles.
BindGlobal("KOAHSS_StackingLowPhaseValue",function(b,c,bp,cp,s,sigma)
    local bits, index, value;
    bits := [b(sigma{[1]}),bp(sigma{[1]}),c(sigma{[1,2]}),
        c(sigma{[2,3]}),cp(sigma{[1,2]}),cp(sigma{[2,3]}),
        s(sigma{[1,2]}),s(sigma{[2,3]})];
    if not ForAll(bits,x->x in [0,1]) then
        Error("koFull: low stacking phase received nonbinary data");
    fi;
    index := 0;
    for value in bits do index := 2*index+value; od;
    value := KOAHSS_STACKING_LOW_PHASE.values[index+1];
    return value[1]/value[2];
end);

BindGlobal("KOAHSS_StackingLowGamma",function(ctx,b,c,bp,cp)
    local phase;
    phase := sigma -> KOAHSS_StackingLowPhaseValue(
        b.value,c.value,bp.value,cp.value,ctx.s.value,sigma);
    # On this exact domain Omega_2=0 and the successor inputs are zero.
    # Thus gamma_2 is precisely the signed coboundary of the selected phase.
    return ctx.make(3,function(sigma)
        local value,j,face;
        value := 0;
        for j in [1..4] do
            face := sigma{Filtered([1..4],i->i<>j)};
            if j=1 then
                value := value+(-1)^ctx.s.value(sigma{[1,2]})*phase(face);
            else value := value+(-1)^(j-1)*phase(face); fi;
        od;
        if not IsInt(value) then
            Error("koFull: the selected closed-input stacking carry is not integral");
        fi;
        return value;
    end);
end);

BindGlobal("KOAHSS_StackingIntegralPrimitive",function(backend,n,rhs)
    local matrix, j, basis;
    matrix := [];
    for j in [1..backend.dimension(n)] do
        basis := List([1..backend.dimension(n)],i->0); basis[j] := 1;
        Add(matrix,backend.coboundary(n,basis,true));
    od;
    return koAHSSSolveIntegerSystem(matrix,rhs);
end);

BindGlobal("KOAHSS_StackingExtensionOracle",function(backend,degree,layers)
    local ctx, layerD, layerC, cache, coordinates, embed, unavailable,
          projectD, binaryLift, closed, answer;
    cache := rec(); ctx := fail;
    if IsRecord(layers) then layers := List(["D","C","B","A"],name->layers.(name)); fi;
    unavailable := reason -> rec(status:="unresolved",reason:=reason);
    if degree in [1,2] and (IsBound(backend.naturalTransport)
        or IsBound(backend.naturalBar)) then
        ctx := KOAHSS_NaturalCochains(backend);
    fi;
    layerD := First(layers,l->l.name="D");
    layerC := First(layers,l->l.name="C");
    coordinates := function(layer,element)
        return CallFuncList(ValueGlobal("KOAHSS_ExtensionGroupCoordinates"),
            [layer.group,layer.generators,element]);
    end;
    embed := function(lower,layer,values)
        local result,stored,j;
        result := List([1..lower.generatorCount],i->0);
        stored := First(lower.layers,l->l.name=layer.name);
        if stored=fail then
            if not IsEmpty(values) then
                Error("koFull: lower presentation omitted a nonzero layer");
            fi;
            return result;
        fi;
        if Length(values)<>Length(layer.generators) then
            Error("koFull: stacking coordinate length does not match layer basis");
        fi;
        for j in [1..Length(values)] do
            result[stored.startColumn+j-1] := values[j];
        od;
        return result;
    end;
    projectD := function(cochain)
        local vector,base;
        vector := ctx.project(cochain,backend.twists.s);
        if ForAny(backend.coboundary(degree+1,vector,true),x->x<>0) then
            Error("koFull: measured stacking relation is not a closed D cochain");
        fi;
        base := backend.cohomologyData(degree+1,-4).class(vector);
        return rec(cochain:=vector,
            coordinates:=coordinates(layerD,layerD.cell.project(base)));
    end;
    binaryLift := function(n,vector)
        return ctx.reduce(ctx.lift(n,vector,0));
    end;
    closed := function(n,vector)
        return ForAll(backend.coboundary(n,vector,false),x->x mod 2=0);
    end;
    answer := function(layer,generatorIndex,quotientOrder,lower)
        local key,b,c,z0,z1,gamma,beta,betaVector,cClass,cCoordinates,
              cSum,dSum,ci,j,t,measured,result,wc,obstruction,obstructionVector,
              primitive,dLift,dFull,relation,homotopy,chosen;
        if ctx=fail then
            return unavailable("stacking lift and lower reduction are implemented only in degrees 1 and 2 on the normalized group-bar backend");
        fi;
        if quotientOrder<>2 then
            return unavailable("the production low-degree stacking adapter supports order-two quotient generators");
        fi;
        if degree=2 and ForAny(backend.twists.omega,x->x<>0) then
            return unavailable("degree-two stacking transport with nonzero omega requires the general defining-system and gauge reducer");
        fi;
        if not ((degree=1 and layer.name="C")
            or (degree=2 and layer.name in ["C","B"])) then
            return unavailable("this layer has no implemented production stacking reducer");
        fi;
        key := Concatenation(layer.name,"_",String(generatorIndex),"_",
            String(lower.presentationId));
        if IsBound(cache.(key)) then return cache.(key); fi;
        chosen := layer.cochains[generatorIndex];
        if not closed(layer.p,chosen) then
            Error("koFull: a leading layer representative is not a binary cocycle");
        fi;
        if degree=1 then
            c := binaryLift(0,chosen);
            wc := ctx.product(ctx.omega,c);
            obstruction := ctx.divide(ctx.d(wc,true),2);
            obstructionVector := ctx.project(obstruction,backend.twists.s);
            primitive := KOAHSS_StackingIntegralPrimitive(backend,2,-obstructionVector);
            if primitive=fail then
                return unavailable("the surviving C generator requires an integral upper lift not supplied by the retained defining system");
            fi;
            dLift := ctx.lift(2,primitive.particular,backend.twists.s);
            homotopy := ctx.make(2,sigma->ctx.transport.homotopyValue(
                3,obstruction.value,sigma));
            dFull := ctx.add([dLift,ctx.scale(-1,homotopy)]);
            relation := ctx.add([ctx.scale(2,dFull),wc]);
            measured := projectD(relation);
            result := rec(status:="computed",lowerPresentationId:=lower.presentationId,
                lowerCoordinates:=embed(lower,layerD,measured.coordinates),
                witness:=rec(operation:="xtimes",power:=2,
                    calibration:="coherent-low-commutative/degree-one",
                    leadingCochain:=ShallowCopy(chosen),
                    upperPrimitive:=primitive.particular,
                    upperObstruction:=obstructionVector,
                    upperHomotopyCorrection:=ctx.project(homotopy,backend.twists.s),
                    lowerD:=measured.cochain,xtimesCalls:=1,
                    lowerReductionProducts:=0));
        else
            z0 := ctx.zero(0); z1 := ctx.zero(1);
            if layer.name="C" then
                c := binaryLift(1,chosen);
                gamma := KOAHSS_StackingLowGamma(ctx,z0,c,z0,c);
                measured := projectD(gamma);
                result := rec(status:="computed",lowerPresentationId:=lower.presentationId,
                    lowerCoordinates:=embed(lower,layerD,measured.coordinates),
                    witness:=rec(operation:="xtimes",power:=2,
                        calibration:=KOAHSS_STACKING_LOW_PHASE.calibration,
                        leadingCochain:=ShallowCopy(chosen),
                        lowerD:=measured.cochain,xtimesCalls:=1,
                        lowerReductionProducts:=0));
            else
                b := binaryLift(0,chosen);
                gamma := KOAHSS_StackingLowGamma(ctx,b,z1,b,z1);
                beta := ctx.product(ctx.s,ctx.product(b,b));
                betaVector := List(ctx.project(beta,0),x->x mod 2);
                cClass := backend.cohomologyData(1,-2).class(betaVector);
                cCoordinates := coordinates(layerC,layerC.cell.project(cClass));
                if not ForAll(cCoordinates,x->x in [0,1]) then
                    Error("koFull: degree-one binary layer coordinates are not binary");
                fi;
                cSum := z1; dSum := ctx.zero(3); t := 0;
                for j in [1..Length(cCoordinates)] do
                    if cCoordinates[j]=1 then
                        ci := binaryLift(1,layerC.cochains[j]);
                        dSum := ctx.add([dSum,
                            KOAHSS_StackingLowGamma(ctx,z0,cSum,z0,ci)]);
                        cSum := ctx.reduce(ctx.add([cSum,ci])); t := t+1;
                    fi;
                od;
                # In normalized group cochains B^1(F2)=0. Equality of H^1
                # classes therefore identifies these complete binary cochains,
                # and no unrecorded C-gauge choice enters this reduction.
                if ctx.project(ctx.reduce(ctx.add([beta,cSum])),0)
                    <>List([1..backend.dimension(1)],i->0) then
                    return unavailable("the supplied degree-one representatives require a nontrivial C-gauge transport");
                fi;
                measured := projectD(ctx.add([gamma,ctx.scale(-1,dSum)]));
                result := rec(status:="computed",lowerPresentationId:=lower.presentationId,
                    lowerCoordinates:=embed(lower,layerD,measured.coordinates)
                        +embed(lower,layerC,cCoordinates),
                    witness:=rec(operation:="xtimes",power:=2,
                        calibration:=KOAHSS_STACKING_LOW_PHASE.calibration,
                        leadingCochain:=ShallowCopy(chosen),lowerC:=betaVector,
                        lowerCCoordinates:=cCoordinates,
                        lowerD:=measured.cochain,xtimesCalls:=1,
                        lowerReductionProducts:=t,
                        reductionCarry:=ctx.project(dSum,backend.twists.s)));
            fi;
        fi;
        cache.(key) := result;
        return result;
    end;
    return answer;
end);
