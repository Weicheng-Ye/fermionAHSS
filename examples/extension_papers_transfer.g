# Historical entry-point alias: extension_papers.g now uses native-resolution
# searches by default. Run from the package root in a fresh GAP process.
CallFuncList(function()
    local file,slashes,prefix;
    file:=INPUT_FILENAME(); slashes:=Positions(file,'/'); prefix:="";
    if not IsEmpty(slashes) then prefix:=file{[1..Last(slashes)]}; fi;
    Read(Concatenation(prefix,"extension_papers.g"));
end,[]);
