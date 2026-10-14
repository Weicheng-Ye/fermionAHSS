# A persistent cache of computed degrees: the JSON of a koFull degree result
# under the hash of everything that determines it (the package sources, the
# resolution through the degree it needs, the twists and the degree). A hit
# returns the stored invariants, presentation and witnesses as data; live
# worker objects are not stored. The cache lives in the package's cache/
# directory, or in FERMIONAHSS_CACHE_DIR when that is set (an empty value
# disables it). Delete the directory to recompute everything.

BindGlobal("KOAHSS_PackageRoot",CallFuncList(function()
    local path,slash;
    path:=INPUT_FILENAME();
    if path[1]<>'/' then path:=Filename(DirectoryCurrent(),path); fi;
    slash:=Last(Positions(path,'/'));
    return Concatenation(path{[1..slash]},"../");
end,[]));

BindGlobal("KOAHSS_CacheDirectory",function()
    local configured,path;
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_CACHE_DIR) then
        configured:=GAPInfo.SystemEnvironment.FERMIONAHSS_CACHE_DIR;
        if IsEmpty(configured) then return fail; fi;
        return Directory(configured);
    fi;
    path:=Concatenation(KOAHSS_PackageRoot,"cache");
    if not IsDirectoryPath(path) then
        if not IsBool(CreateDir(path)) then return fail; fi;
    fi;
    if not IsWritableFile(path) then return fail; fi;
    return Directory(path);
end);

# The hash of the package sources, computed once per session.
BindGlobal("KOAHSS_SourceHash",CallFuncList(function()
    local value;
    value:=fail;
    return function()
        local names,text,directory,entries,entry,sub;
        if value<>fail then return value; fi;
        names:=[];
        for sub in ["gap","python","python/stacking_model","data"] do
            directory:=Concatenation(KOAHSS_PackageRoot,sub);
            entries:=DirectoryContents(directory);
            for entry in entries do
                if EndsWith(entry,".gi") or EndsWith(entry,".gd") or EndsWith(entry,".py")
                   or EndsWith(entry,".json") or EndsWith(entry,".g") then
                    Add(names,Concatenation(sub,"/",entry));
                fi;
            od;
        od;
        Sort(names);
        text:="";
        for entry in names do
            Append(text,entry); Append(text,":");
            Append(text,HexSHA256(StringFile(Concatenation(KOAHSS_PackageRoot,entry))));
            Append(text,"\n");
        od;
        value:=HexSHA256(text);
        return value;
    end;
end,[]));

# The fingerprint of a resolution through degree depth: its ranks and the
# boundary words with their group elements printed.
BindGlobal("KOAHSS_ResolutionFingerprint",function(R,depth)
    local text,n,j,term;
    text:=Concatenation("group:",String(R!.group),"\n");
    for n in [0..depth] do
        Append(text,Concatenation("dim",String(n),":",String(R!.dimension(n)),"\n"));
        if n>0 then
            for j in [1..R!.dimension(n)] do
                for term in R!.boundary(n,j) do
                    Append(text,Concatenation(String(term[1]),"*",String(R!.elts[term[2]]),";"));
                od;
                Append(text,"\n");
            od;
        fi;
    od;
    return HexSHA256(text);
end);

BindGlobal("KOAHSS_ResultCacheKey",function(context,degree)
    local backend,R,depth,text;
    backend:=context.backend;
    if not IsBound(context.resolution) then return fail; fi;
    R:=context.resolution;
    depth:=Minimum(ValueGlobal("EvaluateProperty")(R,"length"),context.cohomologyCapacity+1);
    text:=Concatenation("v2|",KOAHSS_SourceHash(),"|",KOAHSS_ResolutionFingerprint(R,depth),
        "|s=",String(backend.twists.s),"|omega=",String(backend.twists.omega),
        "|k=",String(degree),"|capacity=",String(context.cohomologyCapacity));
    return HexSHA256(text);
end);

# Data only: live groups, maps and worker objects are omitted.
BindGlobal("KOAHSS_JSONData",function(value)
    local result,entry,name,converted;
    if value=fail or IsFunction(value) or IsGroup(value) then return fail; fi;
    if IsInt(value) or IsBool(value) or IsStringRep(value) then return value; fi;
    if IsRat(value) then return String(value); fi;
    if IsList(value) then
        result:=[];
        for entry in value do
            converted:=KOAHSS_JSONData(entry);
            if converted=fail then Add(result,rec(kind:="omitted-non-json-value"));
            else Add(result,converted); fi;
        od;
        return result;
    fi;
    if IsRecord(value) then
        result:=rec();
        for name in RecNames(value) do
            if name in ["group","cell","backend","_context","resolution","generators",
                "worker","process","cochains","fullLifts","pageData","getCell","getMap"] then
                continue;
            fi;
            converted:=KOAHSS_JSONData(value.(name));
            if converted<>fail then result.(name):=converted; fi;
        od;
        return result;
    fi;
    return fail;
end);

BindGlobal("KOAHSS_ResultCacheRead",function(key)
    local directory,path,text;
    directory:=KOAHSS_CacheDirectory();
    if directory=fail or key=fail then return fail; fi;
    path:=Filename(directory,Concatenation("result-",key,".json"));
    if not IsExistingFile(path) then return fail; fi;
    text:=StringFile(path);
    if text=fail then return fail; fi;
    return JsonStringToGap(text);
end);

BindGlobal("KOAHSS_ResultCacheWrite",function(key,answer)
    local directory,path,stream;
    directory:=KOAHSS_CacheDirectory();
    if directory=fail or key=fail then return; fi;
    path:=Filename(directory,Concatenation("result-",key,".json"));
    stream:=OutputTextFile(path,false);
    if stream=fail then return; fi;
    SetPrintFormattingStatus(stream,false);
    WriteAll(stream,GapToJsonString(KOAHSS_JSONData(answer)));
    CloseStream(stream);
end);
