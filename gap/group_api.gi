# Finite groups or supplied integral HAP resolutions use fixed bar transport.
# Nonzero twist vectors are coordinates in this resolution, not abstract
# cohomology-class labels; use the explicit-resolution API to control a basis.
InstallGlobalFunction(koAHSS, function(arg)
    local group, k, count, depth, resolution, space, pageArgs, options,
        hasCount, context, tables, result, completed, modelId, calibrationId,
        suppliedResolution,available;
    if not Length(arg) in [4..6] then
        Error("usage: koAHSS(group or HAP resolution, s, omega, k[, n][, options])");
    fi;
    suppliedResolution:=IsBoundGlobal("IsHapResolution") and
        CallFuncList(ValueGlobal("IsHapResolution"),[arg[1]]);
    if suppliedResolution then
        resolution:=arg[1]; group:=resolution!.group;
    else
        group:=arg[1];
        if not IsGroup(group) then
            Error("koAHSS: the first argument must be a GAP group or an integral HAP resolution");
        fi;
    fi;
    k := arg[4];
    if not IsInt(k) or not k in [-1..6] then
        Error("koAHSS: k = max(p+q+3) must be an integer in [-1..6]");
    fi;
    count := 5; hasCount := false; options := rec();
    if Length(arg)=5 and IsRecord(arg[5]) then
        options := arg[5];
    elif Length(arg)>=5 then
        hasCount := true;
        count := arg[5];
        if not IsInt(count) or not count in [1..5] then
            Error("koAHSS: n counts pages beginning with E2 and must be in [1..5]");
        fi;
    fi;
    if Length(arg)=6 then
        if not IsRecord(arg[6]) then
            Error("koAHSS: the final options argument must be a record");
        fi;
        options := arg[6];
    fi;
    if not ForAll(RecNames(options), name -> name="details") then
        Error("koAHSS: unknown option; supported options are details");
    fi;
    if IsBound(options.details) and not IsBool(options.details) then
        Error("koAHSS: details must be true or false");
    fi;
    if not suppliedResolution and not IsFinite(group) then
        Error("koAHSS: automatic resolution construction requires a finite group; use koAHSSHAPSpace with an explicit integral resolution");
    fi;
    if LoadPackage("hap")=fail then Error("fermionAHSS requires the GAP package hap"); fi;
    if count<=2 then depth:=Maximum(3,k+2); else depth:=Maximum(3,k+3); fi;
    if suppliedResolution then
        available:=ValueGlobal("EvaluateProperty")(resolution,"length");
        if not IsInt(available) or available<depth then
            Error("koAHSS: supplied resolution needs length at least ",depth);
        fi;
    else
        resolution := CallFuncList(ValueGlobal("ResolutionFiniteGroup"),[group,depth]);
        if resolution=fail then Error("koAHSS: HAP could not construct the finite-group resolution"); fi;
    fi;
    space := koAHSSHAPSpace(resolution,koAHSSNaturalOperations());
    pageArgs := [space,arg[2],arg[3],k];
    if hasCount then Add(pageArgs,count); fi;
    if not IsBound(options.details) or not options.details then
        return CallFuncList(koAHSSpages,pageArgs);
    fi;
    context := KOAHSS_ComputePageContext(pageArgs);
    tables := List(context.pageData, KOAHSS_ClassifyPage);
    completed := ForAll(context.pageData,
        table -> ForAll(table, row -> ForAll(row, cell -> not KOAHSS_IsUnresolved(cell))));
    modelId := "normalized-group-bar-five-row";
    calibrationId := "chi7_tail/Danus/zeta010/R2-cubic/muR0/prime3-2";
    context.resolution := resolution;
    context.suppliedResolution:=suppliedResolution;
    context.group := group;
    context.modelId := modelId;
    context.calibrationId := calibrationId;
    context.twists := context.backend.twists;
    result := rec(kind := "koAHSSResult", maxDegree := k,
        computedThrough := context.computedThrough,
        pages := rec(kind := "koAHSSPages",
            pageNumbers := ShallowCopy(context.pageNumbers), tables := tables),
        pageData := context.pageData, status := "unresolved",
        modelId := modelId, twists := context.twists,
        calibrationId := calibrationId, _context := context,
        scope := "five-row-associated-graded", certified_ko := false);
    if completed then result.status := "computed"; fi;
    return result;
end);
