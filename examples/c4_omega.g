# Run from the package root in a fresh process (about 20 minutes; 35 with
# FERMIONAHSS_LIGHT_RELATIONS=0):
# gap -q --quitonbreak examples/c4_omega.g
#
# C4 with the degree-two twist omega=[1] at package degree 5: the largest
# two-primary extension among the small cyclic examples. The E6 line has
# A=Z/4 at (2,0), B=Z/2 at (3,-1), C=Z/2 at (4,-2) and D=Z/4 at (6,-4), and
# the measured relations assemble them to Z/2 + Z/32. The relation of the A
# generator has its target in D: its row is the light A-over-D residue, whose
# degree-five diagonal phases evaluate the universal pair sources (in the
# transferred model with FERMIONAHSS_LIGHT_RELATIONS=0).
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local R, full, vector, witness;
    R := ResolutionFiniteGroup(CyclicGroup(4), 8);
    full := koFull(R, 0, [1], 5);
    Assert(0, full.status = "computed" and full.invariants = [2, 32]);
    koAHSSDisplay(full);
    Print("relation matrix (rows n*g - lower coordinates):\n", full.degreeResult.relationMatrix, "\n");
    for vector in full.degreeResult.extensionVectors do
        witness := vector.result.witness;
        Print(vector.generatorId, " of order ", vector.order, ": ", vector.order, "*g = ",
            vector.result.lowerCoordinates, " in ", vector.lowerGeneratorIds);
        if IsBound(witness.measuredLayers) then
            Print(", measured in the layers ", witness.measuredLayers);
        fi;
        Print("\n");
    od;
    Print("C4, omega=[1], package degree 5: ", full.invariants, "\n");
end, []);
