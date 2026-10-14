# Full group constructors for the koAHSS batch runner.  In particular, a
# space group is never replaced by its finite point/holonomy group here.
BindGlobal("KOAHSSBatchSGCSource",function()
    local path,info;
    info:=PackageInfo("SpaceGroupCohomology");
    if IsEmpty(info) then
        Error("koAHSS batch: install the local SpaceGroupCohomology package");
    fi;
    info:=info[1];
    path:=Filename(Directory(info.InstallationPath),"gap/homotopy.gi");
    if not IsExistingFile(path) then
        Error("koAHSS batch: SpaceGroupCohomology is missing gap/homotopy.gi");
    fi;
    return rec(package:="SpaceGroupCohomology",version:=info.Version,
        module:="gap/homotopy.gi",module_path:=path,
        module_sha256:=HexSHA256(StringFile(path)));
end);

BindGlobal("KOAHSSBatchLoadSGC",function()
    local dependency,source;
    for dependency in ["hap","hapcryst","polymaking","crystcat"] do
        if LoadPackage(dependency)=fail then
            Error("koAHSS batch: space-group resolutions require ",dependency);
        fi;
    od;
    source:=KOAHSSBatchSGCSource();
    # This self-contained module is the package's actual resolution engine.
    # Loading it directly avoids the unrelated ring tables in read.g (the
    # local 2.3.0 tables currently refer to an undefined global in225).
    # Read the exact file being fingerprinted even if an identically named
    # function was already loaded from elsewhere in this GAP session.
    Read(source.module_path);
    if not IsBoundGlobal("SGC_ResolutionSpaceGroup")
       or not IsFunction(ValueGlobal("SGC_ResolutionSpaceGroup")) then
        Error("koAHSS batch: SGC_ResolutionSpaceGroup is not a function");
    fi;
    return source;
end);

BindGlobal("KOAHSSBatchMakeGroup",function(spec,length)
    local group,resolution,result,affine,catalog,number,symbol,
          orientation,property,isResolution,initialGroup,source,cyclicGroups,
          coordinateGenerators,coordinatePowers,projections,coordinateOrders;
    if not IsRecord(spec) or not IsBound(spec.kind)
       or not IsString(spec.kind) then
        Error("koAHSS batch: group specification needs a string kind");
    fi;
    if not IsInt(length) or length<1 then
        Error("koAHSS batch: resolution length must be a positive integer");
    fi;
    if LoadPackage("hap")=fail then Error("koAHSS batch: HAP is required"); fi;
    if spec.kind="abelian" then
        if not IsBound(spec.factors) or not IsList(spec.factors)
           or not ForAll(spec.factors,n->IsInt(n) and n>0) then
            Error("koAHSS batch: finite abelian factors must be positive integers");
        fi;
        # Retain actual projections onto the factors named by the paper. A Pc
        # basis is not a cyclic-product coordinate system (C6 already has two
        # Pc generators), and an H2 basis does not identify an extension class.
        coordinateOrders:=ShallowCopy(spec.factors);
        cyclicGroups:=List(coordinateOrders,m->CyclicGroup(IsPcGroup,m));
        coordinateGenerators:=List([1..Length(cyclicGroups)],i->
            First(Elements(cyclicGroups[i]),g->Order(g)=coordinateOrders[i]));
        coordinatePowers:=List([1..Length(cyclicGroups)],i->
            List([0..coordinateOrders[i]-1],j->coordinateGenerators[i]^j));
        if IsEmpty(cyclicGroups) then
            initialGroup:=AbelianGroup(IsPcGroup,[]); projections:=[];
        else
            initialGroup:=DirectProduct(cyclicGroups);
            projections:=List([1..Length(cyclicGroups)],i->Projection(initialGroup,i));
        fi;
        resolution:=CallFuncList(ValueGlobal("ResolutionFiniteGroup"),[initialGroup,length]);
        result:=rec(kind:="abelian",factors:=ShallowCopy(spec.factors),
            name:=Concatenation("AbelianGroup(",String(spec.factors),")"),
            resolutionMethod:="ResolutionFiniteGroup");
        result.productCoordinates:=function(g)
            local answer,i,value;
            answer:=[];
            for i in [1..Length(projections)] do
                value:=Position(coordinatePowers[i],Image(projections[i],g));
                if value=fail then Error("koAHSS batch: invalid cyclic product coordinate"); fi;
                Add(answer,value-1);
            od;
            return answer;
        end;
        result.productCoordinateConvention:=rec(
            orders:=coordinateOrders,
            factor_generators:=List(coordinateGenerators,String),
            convention:="direct-product projections; least nonnegative powers of the displayed cyclic generators");
        orientation:=g->0;
    elif spec.kind="spacegroup" then
        if not IsBound(spec.number) or not IsInt(spec.number)
           or not spec.number in [1..230] then
            Error("koAHSS batch: a three-dimensional IT number must lie in [1..230]");
        fi;
        source:=KOAHSSBatchLoadSGC();
        number:=spec.number;
        # SGC expects the full affine group in the right-action convention.
        # Its resolution itself uses left-action matrices after standardizing.
        affine:=CallFuncList(ValueGlobal("SpaceGroupOnRightBBNWZ"),[3,number]);
        resolution:=CallFuncList(ValueGlobal("SGC_ResolutionSpaceGroup"),[affine,length]);
        catalog:=ValueGlobal("CrystGroupsCatalogue")[3];
        symbol:=catalog.HermannMauguinSymbol[number];
        result:=rec(kind:="spacegroup",itNumber:=number,
            name:=Concatenation("SG ",String(number)," ",symbol),
            hermannMauguinSymbol:=symbol,affineGroup:=affine,
            catalogueParameters:=ShallowCopy(CallFuncList(ValueGlobal("CrystCatRecord"),[affine]).parameters),
            catalogueITNumbers:=Positions(catalog.spaceGroupTypeInternatTable,
                catalog.spaceGroupTypeInternatTable[number]),
            resolutionMethod:="SGC_ResolutionSpaceGroup",resolutionSource:=source);
        orientation:=function(g)
            local determinant;
            if not g in group then
                Error("koAHSS batch: orientation input is outside the full space group");
            fi;
            # The homogeneous affine matrix has final row (0,0,0,1), so
            # its determinant is the determinant of the linear action. Use
            # the actual resolution group, including SGC's change of basis.
            determinant:=DeterminantMat(g);
            if determinant=1 then return 0;
            elif determinant=-1 then return 1;
            fi;
            Error("koAHSS batch: an affine space-group element has determinant other than +/-1");
        end;
    else
        Error("koAHSS batch: supported kinds are abelian and spacegroup");
    fi;
    isResolution:=ValueGlobal("IsHapResolution");
    property:=ValueGlobal("EvaluateProperty");
    if resolution=fail or not isResolution(resolution)
       or not IsBound(resolution!.homotopy) or not IsFunction(resolution!.homotopy)
       or property(resolution,"characteristic")<>0
       or property(resolution,"length")<length or resolution!.dimension(0)<1 then
        Error("koAHSS batch: constructor did not return an integral free resolution with a contracting homotopy");
    fi;
    group:=resolution!.group;
    if spec.kind="abelian" and group<>initialGroup then
        Error("koAHSS batch: HAP changed the finite group underlying the paper coordinates");
    fi;
    if spec.kind="spacegroup" then
        if not IsMatrixGroup(group) or Length(One(group))<>4 then
            Error("koAHSS batch: SGC did not return a three-dimensional affine matrix group");
        fi;
    fi;
    result.group:=group;
    result.resolution:=resolution;
    result.orientationCharacter:=orientation;
    result.orientationConvention:="0 for determinant +1, 1 for determinant -1";
    return result;
end);
