# Run from the package directory with gap -q --quitonbreak examples/twisted_c2.g.
# Also works via ReadPackage("fermionAHSS", "examples/twisted_c2.g").
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash;
        file := INPUT_FILENAME(); slash := Positions(file, '/');
        Read(Concatenation(file{[1..Last(slash)]}, "../load.g"));
    end, []);
fi;

# HAP's standard cyclic C2 resolution has one basis element in each degree.
# s=[1] is its sign character; omega=[1] is its degree-two binary generator.
twistedC2Pages := koAHSS(CyclicGroup(2), [1], [1], 1, 5);;
Assert(0, Length(twistedC2Pages) = 5);
Print("C2, s=[1], omega=[1], E2--E6: ", twistedC2Pages, "\n");

# Reuse a chosen resolution when the coordinates of the twists matter.
c2Resolution := ResolutionFiniteGroup(CyclicGroup(2), 4);;
c2Space := koAHSSHAPSpace(c2Resolution, koAHSSNaturalOperations());;
Assert(0, twistedC2Pages = koAHSSpages(c2Space, [1], [1], 1, 5));
