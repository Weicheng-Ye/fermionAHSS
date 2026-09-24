# Finite groups use an integral HAP resolution with fixed bar transport.
# Nonzero twist vectors are coordinates in this resolution, not abstract
# cohomology-class labels; use the explicit-resolution API to control a basis.
InstallGlobalFunction(koAHSS, function(arg)
    local group, k, count, depth, resolution, space, pageArgs;
    if not Length(arg) in [4,5] then
        Error("usage: koAHSS(group, s, omega, k[, n])");
    fi;
    group := arg[1];
    if not IsGroup(group) then Error("koAHSS: the first argument must be a GAP group"); fi;
    k := arg[4];
    if not IsInt(k) or not k in [-1..6] then
        Error("koAHSS: k = max(p+q+3) must be an integer in [-1..6]");
    fi;
    count := 5;
    if Length(arg)=5 then
        count := arg[5];
        if not IsInt(count) or not count in [1..5] then
            Error("koAHSS: n counts pages beginning with E2 and must be in [1..5]");
        fi;
    fi;
    if not IsFinite(group) then
        Error("koAHSS: automatic resolution construction requires a finite group; use koAHSSHAPSpace with an explicit integral resolution");
    fi;
    if LoadPackage("hap")=fail then Error("fermionAHSS requires the GAP package hap"); fi;
    if count<=2 then depth:=Maximum(3,k+2); else depth:=Maximum(3,k+3); fi;
    resolution := CallFuncList(ValueGlobal("ResolutionFiniteGroup"),[group,depth]);
    if resolution=fail then Error("koAHSS: HAP could not construct the finite-group resolution"); fi;
    space := koAHSSHAPSpace(resolution,koAHSSNaturalOperations());
    pageArgs := [space,arg[2],arg[3],k];
    if Length(arg)=5 then Add(pageArgs,count); fi;
    return CallFuncList(koAHSSpages,pageArgs);
end);
