# The same literature fixtures using native-resolution searches and exact
# completion checks. Run from the package root in a fresh GAP process.
BindGlobal("FERMIONAHSS_EXTENSION_MODEL","transfer");
CallFuncList(function()
    local file,slashes,prefix;
    file:=INPUT_FILENAME(); slashes:=Positions(file,'/'); prefix:="";
    if not IsEmpty(slashes) then prefix:=file{[1..Last(slashes)]}; fi;
    Read(Concatenation(prefix,"extension_papers.g"));
end,[]);
