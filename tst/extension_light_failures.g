# Fault injection at the worker and final-audit boundaries.
CallFuncList(function()
    local originalModel,faultKind,calls,ahss,fallback,refused,originalFrame,audited;
    originalModel := KOAHSS_ExtensionTransferredModel;;
    faultKind := "error";;
    calls := 0;;
    MakeReadWriteGlobal("KOAHSS_ExtensionTransferredModel");
    KOAHSS_ExtensionTransferredModel := function(backend,k)
        local model,originalLight;
        model := originalModel(backend,k);
        if model.status <> "computed" then return model; fi;
        originalLight := model.light;
        model.light := function(task,data)
            calls := calls+1;
            if calls=2 then
                if faultKind="error" then
                    return rec(status:="error",reason:="injected error",traceback:="");
                fi;
                return rec(status:="unresolved",code:="resource-limit",reason:="injected resource limit");
            fi;
            return originalLight(task,data);
        end;
        return model;
    end;;
    ahss := koAHSS_batch(CyclicGroup(4),[1],0,6,rec(details:=true));;
    fallback := koFull(ahss);;
    Assert(0,fallback.status="computed" and fallback.invariants=[4]);
    Assert(0,Length(fallback.degreeResult.lightFallbacks)=1);
    Assert(0,fallback.degreeResult.lightFallbacks[1].code="light-error");
    Assert(0,fallback.degreeResult.certificateLevel="transfer-R");
    Assert(0,fallback.degreeResult.heavyMeasurements>0);

    faultKind := "resource";; calls := 0;;
    refused := koFull(ahss);;
    Assert(0,refused.status="unresolved");
    Assert(0,refused.degreeResult.pendingLayer="resource-limit");
    Assert(0,not IsBound(refused.degreeResult.lightFallbacks));
    Assert(0,refused.degreeResult.heavyMeasurements=0);

    KOAHSS_ExtensionTransferredModel:=originalModel;;
    originalFrame:=KOAHSS_LightFrame;;
    MakeReadWriteGlobal("KOAHSS_LightFrame");
    KOAHSS_LightFrame:=function(env)
        local frame; frame:=originalFrame(env);
        frame.audit:=function(arg) return "injected stale frame"; end;
        return frame;
    end;;
    audited:=koFull(ahss);;
    Assert(0,audited.status="computed" and audited.invariants=[4]);
    Assert(0,audited.degreeResult.lightFallbacks[1].code="frame-audit");
    Assert(0,audited.degreeResult.heavyMeasurements>0);

    KOAHSS_ExtensionTransferredModel:=originalModel;
    KOAHSS_LightFrame:=originalFrame;
    MakeReadOnlyGlobal("KOAHSS_ExtensionTransferredModel");
    MakeReadOnlyGlobal("KOAHSS_LightFrame");
end,[]);
