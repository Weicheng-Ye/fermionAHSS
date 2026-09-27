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
    local order,R,sign,ahss,native,degree,expected;
    for order in [2,4,8] do
        R:=ResolutionFiniteGroup(CyclicGroup(order),6);
        sign:=0; if order=4 then sign:=[1]; fi;
        ahss:=koAHSS_batch(R,sign,0,3,rec(details:=true));
        native:=koFull_batch(ahss);
        Assert(0,IsIdenticalObj(native.ahss._context.resolution,R));
        Assert(0,native.status="computed");
        degree:=native.degreeResults[5];
        Assert(0,degree.modelId="transferred-normalized-bar" and degree.gaugeCompletenessAssumed);
        expected:=[0,8];
        if order=4 then expected:=[2]; elif order=8 then expected:=[0,2,16]; fi;
        Assert(0,degree.invariants=expected);
        Print("C",order,", package degrees -1..3: ",native.invariants,
            "; native-resolution computation verified\n");
    od;
end,[]);
