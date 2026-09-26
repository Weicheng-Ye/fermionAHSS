# Run from the package root in a fresh process:
# gap -q --quitonbreak examples/resolution_extensions.g
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file,slashes,prefix;
        file:=INPUT_FILENAME(); slashes:=Positions(file,'/'); prefix:="";
        if not IsEmpty(slashes) then prefix:=file{[1..Last(slashes)]}; fi;
        Read(Concatenation(prefix,"../load.g"));
    end,[]);
fi;
CallFuncList(function()
    local order,R,sign,ahss,native,reference,degree;
    for order in [2,4] do
        R:=ResolutionFiniteGroup(CyclicGroup(order),6);
        sign:=0; if order=4 then sign:=[1]; fi;
        ahss:=koAHSS(R,sign,0,3,rec(details:=true));
        native:=koFull(ahss,rec(extensionModel:="transfer"));
        reference:=koFull(ahss,rec(extensionModel:="bar"));
        Assert(0,IsIdenticalObj(native.ahss._context.resolution,R));
        Assert(0,native.status="computed" and reference.status="computed");
        Assert(0,native.invariants=reference.invariants);
        degree:=native.degreeResults[5];
        Assert(0,degree.modelSelection.selected="transfer" and not degree.modelSelection.fallback);
        if order=4 then
            Assert(0,degree.barCertification.status="computed");
            Assert(0,not degree.barCertification.gaugeCompletenessAssumed);
        fi;
        Print("C",order,", package degrees -1..3: ",native.invariants,
            "; native/reference agreement verified\n");
    od;
end,[]);
