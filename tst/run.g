# Run every test of the package in this session: gap -q --quitonbreak -T tst/run.g
CallFuncList(function()
    local file, slash, directory, name, ok;
    file := INPUT_FILENAME(); slash := Positions(file, '/');
    directory := file{[1..Last(slash)]};
    Read(Concatenation(directory, "../load.g"));
    ok := true;
    for name in ["smoke.tst"] do
        if not Test(Concatenation(directory, name), rec(compareFunction := "uptowhitespace")) then
            ok := false; Print("FAILED: ", name, "\n");
        fi;
    od;
    if ok then Print("all tests passed\n"); fi;
end, []);
QUIT_GAP(0);
