# Run from the package root in a fresh process (about three minutes):
# gap -q --quitonbreak examples/twisted_s3.g
#
# The symmetric group S3 with the sign twist at package degree 5. A named
# class is not a coordinate vector: the twist s is the cocycle representing
# the generator of H^1(S3;Z/2) in the basis of the chosen resolution, read
# from the cohomology data of the untwisted backend. With this twist the E6
# line has A=Z/3 at (2,0) and D=Z/3 at (6,-4), the anti-invariant classes of
# the 3-Sylow subgroup, and the three-local relation gives Z/9.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

CallFuncList(function()
    local R, space, backend, H1, generators, s, full, line;
    R := ResolutionFiniteGroup(SymmetricGroup(3), 8);
    space := koAHSSHAPSpace(R, koAHSSNaturalOperations());
    backend := space.koAHSS(0, 0, 8);
    H1 := backend.cohomologyData(1, -1);
    generators := IndependentGeneratorsOfAbelianGroup(H1.group);
    Assert(0, Length(generators) = 1);
    s := H1.represent(generators[1]);
    Print("sign twist of S3 in the resolution basis: ", s, "\n");
    line := koAHSS(R, s, 0, 5).lines[1];
    Print("E6 line of degree 5 (q=-4..0): ", line, "\n");
    full := koFull(R, s, 0, 5);
    Assert(0, full.status = "computed" and full.invariants = [9]);
    koAHSSDisplay(full);
    Print("S3 with the sign twist, package degree 5: ", full.invariants, "\n");
end, []);
