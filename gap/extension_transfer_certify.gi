# Independently certify a native-resolution presentation in the reference model.
# This deliberately does not infer gauge-class completeness from a retraction.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.
BindGlobal("KOAHSS_ExtensionTransferCertify",function(backend,model,k,layers,result)
    local bar,run,attempt,oldBreak,answer;
    if result.status<>"computed" or not IsBound(model.barState) then
        return rec(status:="unresolved",
            reason:="the native candidate has no complete-bar representative export");
    fi;
    bar:=KOAHSS_ExtensionBarModel(backend,k);
    if bar.status<>"computed" then
        return rec(status:="unresolved",reason:=bar.reason,
            certificateLevel:="complete-bar");
    fi;
    run:=function()
        local barLayers,name,j,state,generators,ids,relation,index,power,
            canonical,coefficient,chosen,proof,proofs,audit,source,earlier,
            expected,fields,indexInFields;
        barLayers:=rec(); generators:=[]; ids:=[];
        fields:=["A","B","C","D"];
        for name in ["D","C","B","A"] do
            barLayers.(name):=ShallowCopy(layers.(name));
            barLayers.(name).fullLifts:=[];
            for j in [1..Length(layers.(name).orders)] do
                source:=layers.(name).fullLifts[j];
                if source.status<>"computed" then
                    return rec(status:="unresolved",
                        reason:="a native generator has no fixed full representative");
                fi;
                state:=model.barState(k,source.state,bar);
                indexInFields:=Position(fields,name);
                if ForAny(fields{[1..indexInFields-1]},earlier->
                    ForAny(state.(earlier),x->x<>0)) then
                    Error("koFull: transfer export introduced a layer above its marked generator");
                fi;
                expected:=bar.lift(layers.(name).p,layers.(name).cochains[j],name in ["A","D"]);
                if state.(name)<>expected then
                    Error("koFull: transfer export changed its marked E6 leading representative");
                fi;
                if not KOAHSS_ExtensionStateIsZero(bar.d(k,state)) then
                    Error("koFull: a transferred generator is not flat on the complete bar");
                fi;
                MakeImmutable(state);
                Add(barLayers.(name).fullLifts,rec(status:="computed",state:=state,
                    witness:=rec(flatnessVerified:=true,modelId:=bar.modelId,
                        nativeState:=source.state,sourceModelId:=model.modelId)));
                Add(generators,state); Add(ids,Concatenation(name,":",String(j)));
            od;
        od;
        if ids<>result.lowerModel.generatorIds then
            Error("koFull: transfer certification changed the marked generator basis");
        fi;
        proofs:=[];
        for relation in result.extensionVectors do
            index:=Position(ids,relation.generatorId);
            if index=fail then Error("koFull: an unmarked transfer relation was supplied"); fi;
            power:=KOAHSS_ExtensionPower(bar,k,generators[index],relation.order);
            canonical:=bar.zero(k);
            for j in [1..Length(relation.lowerGeneratorIds)] do
                index:=Position(ids,relation.lowerGeneratorIds[j]);
                if index=fail then Error("koFull: an unmarked lower relation coordinate was supplied"); fi;
                coefficient:=relation.result.lowerCoordinates[j];
                chosen:=KOAHSS_ExtensionPower(bar,k,generators[index],coefficient);
                if not KOAHSS_ExtensionStateIsZero(chosen) then
                    if KOAHSS_ExtensionStateIsZero(canonical) then canonical:=chosen;
                    else canonical:=bar.xtimes(k,canonical,chosen); fi;
                fi;
            od;
            # Keep the original bar equation and its full gauge. A projection
            # of the native certificate is not accepted as this witness.
            proof:=KOAHSS_ExtensionGaugeCompare(bar,k,power,canonical);
            if proof.status<>"computed" then
                return rec(status:="unresolved",
                    reason:="a native relation lacks an independent complete-bar comparison",
                    generatorId:=relation.generatorId,comparison:=proof,
                    completedRelations:=proofs);
            fi;
            Add(proofs,rec(generatorId:=relation.generatorId,order:=relation.order,
                lowerGeneratorIds:=relation.lowerGeneratorIds,
                lowerCoordinates:=relation.result.lowerCoordinates,
                stackedState:=power,canonicalLowerProduct:=canonical,
                canonicalComparison:=proof));
        od;
        audit:=KOAHSS_ExtensionFiniteAudit(bar,k,barLayers,result);
        if audit.status<>"computed" then
            return rec(status:="unresolved",reason:=audit.reason,
                completedRelations:=proofs,algebraAudit:=audit);
        fi;
        return rec(status:="computed",certificateLevel:="complete-bar",
            method:="native lifts with independent reference relations and quotient audit",
            modelId:=bar.modelId,sourceModelId:=model.modelId,
            barLayers:=barLayers,relationWitnesses:=proofs,algebraAudit:=audit,
            gaugeCompletenessAssumed:=false,certified_ko:=false);
    end;
    oldBreak:=BreakOnError; BreakOnError:=false;
    attempt:=CALL_WITH_CATCH(run,[]); BreakOnError:=oldBreak; bar.close();
    if attempt[1] then return attempt[2]; fi;
    if IsBound(bar.lastFailure) and bar.lastFailure.status="unresolved" then
        answer:=ShallowCopy(bar.lastFailure); answer.certificateLevel:="complete-bar";
        return answer;
    fi;
    # The caller handles a native transport resource failure. Invalid exact
    # equations remain errors, rather than being converted into a split group.
    Error("koFull: independent bar certification failed; the original error is reported above");
end);
