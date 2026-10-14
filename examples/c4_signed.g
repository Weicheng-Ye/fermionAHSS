# Run from the package root in a fresh process (under a minute):
# gap -q --quitonbreak examples/c4_signed.g
#
# C4 with the sign twist s=[1] in every package degree from -1 to 6, from one
# E6 calculation on a supplied resolution: the groups are 0, Z/2, Z/2,
# Z/2 + Z/4, Z/2, Z/4, 0 and Z/4. At degree 6 the E6 line has B=Z/2 at
# (4,-1) and D=Z/2 at (7,-4), and the relation 2B=D assembles them to Z/4.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local R, batch, degree, result, relation;
    R := ResolutionFiniteGroup(CyclicGroup(4), 9);
    batch := koFull_batch(R, [1], 0, 6);
    Assert(0, batch.status = "computed");
    koAHSSDisplay(batch);
    for degree in batch.degrees do
        result := batch.degreeResults[degree + 2];
        Print("degree ", degree, ": ", result.invariants, "\n");
    od;
    Assert(0, batch.invariants = [[], [2], [2], [2, 4], [2], [4], [], [4]]);
    relation := First(batch.degreeResults[8].extensionVectors, v -> v.layer = "B");
    Assert(0, relation.order = 2 and relation.result.lowerCoordinates = [1]);
    Print("C4, s=[1], package degrees -1..6: ", batch.invariants, "\n");
end, []);
