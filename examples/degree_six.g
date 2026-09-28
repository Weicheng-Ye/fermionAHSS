# Run from the package root in a fresh process (about ten minutes):
# gap -q --quitonbreak examples/degree_six.g
# Degree six for C4 with s=[1]: the E6 line has B=Z/2 at (4,-1) and D=Z/2 at
# (7,-4), and the measured relation 2B=D assembles them to Z/4.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file,slashes,prefix;
        file:=INPUT_FILENAME(); slashes:=Positions(file,'/'); prefix:="";
        if not IsEmpty(slashes) then prefix:=file{[1..Last(slashes)]}; fi;
        Read(Concatenation(prefix,"../load.g"));
    end,[]);
fi;
CallFuncList(function()
    local R,full,relation;
    R:=ResolutionFiniteGroup(CyclicGroup(4),9);
    full:=koFull(R,[1],0,6);
    Assert(0,full.status="computed" and full.invariants=[4]);
    Assert(0,full.degreeResult.modelId="transferred-normalized-bar");
    relation:=First(full.degreeResult.extensionVectors,v->v.layer="B");
    Assert(0,relation.order=2 and relation.result.lowerCoordinates=[1]);
    koAHSSDisplay(full);
    Print("C4, s=[1], package degree 6: ",full.invariants,"; native-resolution computation verified\n");
end,[]);
