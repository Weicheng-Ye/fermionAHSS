# Run from the package root in a fresh process (about a minute):
# gap -q --quitonbreak examples/suspension.g
#
# An infinite group on a supplied resolution, and degree six at the prime
# three. The resolution of Z/3 x Z is the HAP direct product of the cyclic
# resolution with the resolution of Z. The window of package degree 6 of
# G x Z is the suspension of the window of degree 5 of G, so the E6 line has
# A=Z/3 at (3,0) and D=Z/3 at (7,-4), and the three-local relation of degree
# six gives Z/9, the degree-five value of Z/3. Replacing 3 by 9 gives
# Z/3 + Z/27 (the pages then take about an hour).
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local R, full, vector;
    R := ResolutionDirectProduct(ResolutionFiniteGroup(CyclicGroup(3), 9),
        ResolutionAbelianGroup([0], 9));
    Print("resolution of Z/3 x Z: ranks ", List([0..9], R!.dimension), "\n");
    full := koFull(R, 0, 0, 6);
    Assert(0, full.status = "computed" and full.invariants = [9]);
    koAHSSDisplay(full);
    vector := First(full.degreeResult.extensionVectors, v -> v.layer = "A");
    Assert(0, vector.result.witness.model = "three-local");
    Print("3*a = ", vector.result.lowerCoordinates[1], "*d, measured in the ",
        vector.result.witness.model, " model\n");
    Print("Z/3 x Z, untwisted, package degree 6: ", full.invariants, "\n");
end, []);
