# Run from the package directory with gap -q --quitonbreak examples/c2.g.
# Also works via ReadPackage("fermionAHSS", "examples/c2.g").
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

# Untwisted C2; physical dimension at most 1; all five pages E2--E6.
c2Pages := koAHSS(CyclicGroup(2), 0, 0, 1, 5);;
Assert(0, Length(c2Pages) = 5);
Assert(0, ForAll(c2Pages, table -> List(table, Length) = [3,2,1,0,0]));
Assert(0, c2Pages[5] = koAHSS(CyclicGroup(2), 0, 0, 1));
Print("C2, untwisted, E2--E6: ", "\n");
koAHSSDisplay(c2Pages);
