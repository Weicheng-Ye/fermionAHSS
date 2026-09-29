# Run from the package root in a fresh process (about three minutes):
# gap -q --quitonbreak examples/pin_minus.g
#
# Pin- bordism in five spatial dimensions: C2 with the sign twist s=[1] at
# package degree 6. The E6 line has one Z/2 in each of the layers A (3,0),
# B (4,-1), C (5,-2) and D (7,-4), and the measured relations 2D=0, 2C=D,
# 2B=C, 2A=B assemble them to Z/16. Each relation is read in the layer just
# below its generator, the target layer, so the degree-six stacking
# correction of two A layers is never evaluated.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local full, vector, witness;
    full := koFull(CyclicGroup(2), [1], 0, 6);
    Assert(0, full.status = "computed" and full.invariants = [16]);
    koAHSSDisplay(full);
    for vector in full.degreeResult.extensionVectors do
        witness := vector.result.witness;
        Print(vector.generatorId, " of order ", vector.order, ": ", vector.order, "*g = ",
            vector.result.lowerCoordinates, " in ", vector.lowerGeneratorIds);
        if IsBound(witness.measuredLayers) then
            Print(", measured in the layers ", witness.measuredLayers);
        fi;
        Print("\n");
    od;
    Print("C2, s=[1], package degree 6: ", full.invariants, "\n");
end, []);
