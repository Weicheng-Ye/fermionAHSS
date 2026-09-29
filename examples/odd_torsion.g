# Run from the package root in a fresh process (under a minute):
# gap -q --quitonbreak examples/odd_torsion.g
#
# Odd torsion at package degree 5, where the layers A (2,0) and D (6,-4) carry
# the p+ip and bosonic classes. Relations are measured at the prime of their
# generator: the three-primary ones in the two-layer three-local model, the
# ones at primes five and above split without a measurement.
#   Z/3        Z/9            (3a = 2d modulo three)
#   Z/9        Z/3 + Z/27     (9a = 6d modulo nine)
#   Z/3 x Z/3  Z/3^2 + Z/9^2
#   Z/5        Z/5 + Z/5      (split: the rows lie in different Adams summands)
#   Z/6        Z/9            (the two-primary part of the E6 line is zero)
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local cases, entry, full, vector, witness;
    cases := [
        rec(name := "Z/3", group := CyclicGroup(3), expected := [9]),
        rec(name := "Z/9", group := CyclicGroup(9), expected := [3, 27]),
        rec(name := "Z/3 x Z/3", group := AbelianGroup([3, 3]), expected := [3, 3, 9, 9]),
        rec(name := "Z/5", group := CyclicGroup(5), expected := [5, 5]),
        rec(name := "Z/6", group := CyclicGroup(6), expected := [9])];
    for entry in cases do
        full := koFull(entry.group, 0, 0, 5);
        Assert(0, full.status = "computed" and full.invariants = entry.expected);
        Print(entry.name, ", untwisted, package degree 5: ", full.invariants, "\n");
        for vector in full.degreeResult.extensionVectors do
            witness := vector.result.witness;
            if IsBound(witness.model) then
                Print("  ", vector.generatorId, " of order ", vector.order, ": ", vector.order,
                    "*g = ", vector.result.lowerCoordinates, " in ", vector.lowerGeneratorIds,
                    " (", witness.model, " model)\n");
            fi;
        od;
        Print("  measured relations by prime: ", full.degreeResult.primes, "\n");
    od;
end, []);
