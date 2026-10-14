# Exact page homology uses genuine group homomorphisms.  Missing operations
# produce explicit unresolved cells; they are never silently replaced by zero.

BindGlobal("KOAHSS_IsUnresolved", function(value)
    return IsRecord(value) and IsBound(value.status)
        and value.status = "unresolved";
end);

BindGlobal("KOAHSS_OperationName", function(r, q)
    if r = 2 and q = 0 then return "Dbar";
    elif r = 2 and q = -1 then return "D";
    elif r = 3 and q = -2 then return "Dtilde";
    elif r = 3 and q = 0 then return "Tau";
    elif r = 4 and q = -1 then return "Psi";
    elif r = 5 and q = 0 then return "T";
    fi;
    return fail;
end);

BindGlobal("KOAHSS_ArrowExists", function(r, q)
    return (r = 2 and q in [0, -1])
        or (r = 3 and q in [0, -2])
        or (r = 4 and q = -1)
        or (r = 5 and q = 0);
end);

BindGlobal("KOAHSS_ZeroMap", function(source, target)
    local gens, map;
    gens := GeneratorsOfGroup(source);
    map := GroupHomomorphismByImages(source, target, gens,
        List(gens, x -> One(target)));
    if map = fail then
        Error("koAHSS: cannot construct the zero group homomorphism");
    fi;
    return map;
end);

# GAP's generic infinite fp-group kernel/quotient algorithms can attempt an
# infinite coset enumeration even for Z.  Use exact abelian coordinates before
# page homology, retaining explicit maps to the backend's original E2 group.
BindGlobal("KOAHSS_NormalizeGroup", function(group)
    local generators, free, relations, matrix, smith, transform, inverse,
        orders, active, i, current, currentGens, images, representatives,
        map, productPowers;
    if IsPcpGroup(group) then
        return rec(group := group, lift := x -> x, project := x -> x);
    fi;
    if not IsFpGroup(group) then
        map := IsomorphismPcpGroup(group);
        return rec(group := Range(map), lift := x -> PreImagesRepresentative(map, x),
            project := x -> Image(map, x));
    fi;
    generators := GeneratorsOfGroup(group);
    free := FreeGeneratorsOfFpGroup(group);
    relations := RelatorsOfFpGroup(group);
    matrix := List(relations, rel -> List(free, gen -> ExponentSumWord(rel, gen)));
    if IsEmpty(matrix) or IsEmpty(generators) then
        transform := IdentityMat(Length(generators));
        orders := List(generators, x -> 0);
    else
        smith := NormalFormIntMat(matrix, 1 + 8 + 16);
        transform := smith.coltrans;
        orders := List(generators, x -> 0);
        for i in [1..Minimum(Length(matrix), Length(generators))] do
            orders[i] := AbsInt(smith.normal[i][i]);
        od;
    fi;
    active := Filtered([1..Length(orders)], i -> orders[i] <> 1);
    current := AbelianPcpGroup(orders{active});
    currentGens := GeneratorsOfGroup(current);
    productPowers := function(gens, exponents, identity)
        local result, j;
        result := identity;
        for j in [1..Length(gens)] do result := result * gens[j]^exponents[j]; od;
        return result;
    end;
    images := List([1..Length(generators)], i ->
        productPowers(currentGens, transform[i]{active}, One(current)));
    map := GroupHomomorphismByImages(group, current, generators, images);
    if map = fail then Error("koAHSS: failed to construct abelian Smith coordinates"); fi;
    if IsEmpty(generators) then inverse := []; else inverse := transform^-1; fi;
    representatives := List(active, i -> productPowers(generators, inverse[i], One(group)));
    return rec(group := current,
        lift := x -> productPowers(representatives, Exponents(x), One(group)),
        project := x -> Image(map, x));
end);

BindGlobal("KOAHSS_ComputePageContext", function(args)
    local space, s, omega, bound, count, backend, cells, maps, zeroGroup,
        getCell, getMap, key, emptyCell, validateGroup, makeBaseCell,
        makeDerivedCell, makeUnresolvedCell, mergeUnresolved, isZero,
        lastPage, pageNumbers, result, page, row, q, p, maxDegree;

    if not Length(args) in [4, 5] then
        Error("usage: koAHSSpages(space, s, omega, k[, n])");
    fi;
    space := args[1]; s := args[2]; omega := args[3]; bound := args[4];
    if IsBoundGlobal("IsHapResolution")
       and CallFuncList(ValueGlobal("IsHapResolution"), [space]) then
        space := koAHSSHAPSpace(space);
    fi;
    if not IsInt(bound) or not bound in [-1..6] then
        Error("koAHSS: k = max(p+q+3) must be an integer in [-1..6]");
    fi;
    if Length(args) = 5 then
        count := args[5];
        if not IsInt(count) or not count in [1..5] then
            Error("koAHSS: n counts pages beginning with E2 and must be in [1..5]");
        fi;
        pageNumbers := [2..count + 1];
    else
        pageNumbers := [6];
    fi;
    lastPage := Maximum(pageNumbers);
    if not IsRecord(space) or not IsBound(space.koAHSS)
       or not IsFunction(space.koAHSS) then
        Error("koAHSS: space must provide koAHSS(s, omega, maxDegree)");
    fi;
    # The display ends at p+q+3=k. Outgoing differentials raise that total
    # by one, so their targets must remain available beyond the display.
    # The five-row dependency graph needs cohomology through k+1 for E2/E3,
    # and through k+2 from E4 onward when the q=-2 row is present.
    # This is a preparation hint, not a request for every cohomology cell.
    maxDegree := bound+1;
    if lastPage>=4 and bound>=1 then maxDegree:=bound+2; fi;
    backend := space.koAHSS(s, omega, maxDegree);
    if not IsRecord(backend) or not IsBound(backend.cohomology)
       or not IsFunction(backend.cohomology) then
        Error("koAHSS: backend must provide cohomology(degree, q)");
    fi;
    if not IsBound(backend.differential) or not IsFunction(backend.differential) then
        Error("koAHSS: backend must provide differential(r, degree, q, source, target)");
    fi;
    cells := rec(); maps := rec(); zeroGroup := AbelianPcpGroup([]);
    key := function(r, degree, line)
        return Concatenation("e", String(r), "p", String(degree), "q", String(line));
    end;
    isZero := value -> IsRecord(value) and IsBound(value.status)
        and value.status = "zero";
    mergeUnresolved := function(values)
        local result, value;
        result := rec(status := "unresolved", operations := [], reasons := [],
            dependencies := [], tertiaryCoefficient := 1);
        if IsBound(backend.tertiaryCoefficient) then
            result.tertiaryCoefficient:=backend.tertiaryCoefficient;
        fi;
        if IsBound(backend.tertiaryReference) then
            result.tertiaryReference:=backend.tertiaryReference;
        fi;
        for value in values do
            if KOAHSS_IsUnresolved(value) then
                if not IsBound(value.operations) or not IsList(value.operations)
                   or not ForAll(value.operations, IsString)
                   or not IsBound(value.reasons) or not IsList(value.reasons)
                   or not ForAll(value.reasons, IsString) then
                    Error("koAHSS: unresolved operations and reasons must be lists of strings");
                fi;
                UniteSet(result.operations, Set(value.operations));
                UniteSet(result.reasons, Set(value.reasons));
                if IsBound(value.dependencies) then
                    # Internal dependency entries are [r,p,q,operation].
                    if not IsList(value.dependencies)
                       or not ForAll(value.dependencies, d -> IsList(d)
                           and Length(d) = 4 and ForAll(d{[1..3]}, IsInt)
                           and IsString(d[4])) then
                        Error("koAHSS: unresolved dependencies must be [r,p,q,operation] entries");
                    fi;
                    UniteSet(result.dependencies, Set(value.dependencies));
                fi;
            fi;
        od;
        return result;
    end;
    makeUnresolvedCell := function(r, degree, line, previous, outgoing, incoming)
        local cell;
        cell := mergeUnresolved([previous, outgoing, incoming]);
        cell.page := r; cell.degree := degree; cell.q := line;
        cell.previous := previous;
        cell.incoming := incoming; cell.outgoing := outgoing;
        if IsBound(previous.base) then cell.base := previous.base; fi;
        if KOAHSS_IsUnresolved(previous) then
            if IsBound(previous.lastKnown) then cell.lastKnown := previous.lastKnown; fi;
            if IsBound(previous.lastKnownInvariants) then
                cell.lastKnownInvariants := ShallowCopy(previous.lastKnownInvariants);
                cell.lastKnownPage := previous.lastKnownPage;
            fi;
        else
            cell.lastKnown := previous;
            cell.lastKnownInvariants := AbelianInvariants(previous.group);
            cell.lastKnownPage := previous.page;
        fi;
        return cell;
    end;
    validateGroup := function(group, degree, line)
        if not IsGroup(group) then
            Error("koAHSS: cohomology(", degree, ",", line, ") must return a group");
        fi;
        if not IsAbelian(group) then
            Error("koAHSS: cohomology(", degree, ",", line, ") must be abelian");
        fi;
        # Accessing a finite generating list is also required by the maps.
        GeneratorsOfGroup(group);
    end;
    makeBaseCell := function(group, degree, line)
        local normalized, current;
        normalized := KOAHSS_NormalizeGroup(group);
        current := normalized.group;
        return rec(page := 2, degree := degree, q := line,
            group := current, current := current, base := group,
            lift := function(x)
                if not x in current then
                    Error("koAHSS: lift argument is not in the current cell");
                fi;
                return normalized.lift(x);
            end,
            project := function(x)
                if not x in group then
                    Error("koAHSS: projection argument is not in the E2 cell");
                fi;
                return normalized.project(x);
            end);
    end;
    emptyCell := function(r, degree, line)
        local cell;
        cell := makeBaseCell(zeroGroup, degree, line);
        cell.page := r;
        return cell;
    end;
    makeDerivedCell := function(r, degree, line, previous, outgoing, incoming)
        local cycles, boundaries, imageGens, quotient, group;
        if isZero(outgoing) then cycles := previous.group;
        else cycles := Kernel(outgoing); fi;
        if isZero(incoming) then imageGens := [];
        else imageGens := GeneratorsOfGroup(Image(incoming)); fi;
        if not isZero(outgoing)
           and not ForAll(imageGens, x -> Image(outgoing, x) = One(Range(outgoing))) then
            Error("koAHSS: d", r - 1, " squared is nonzero at (p,q)=(",
                degree, ",", line, ")");
        fi;
        boundaries := Subgroup(cycles, imageGens);
        quotient := NaturalHomomorphismByNormalSubgroup(cycles, boundaries);
        group := Image(quotient);
        return rec(page := r, degree := degree, q := line,
            group := group, current := group, base := previous.base,
            previous := previous, cycles := cycles, boundaries := boundaries,
            quotientMap := quotient, incoming := incoming, outgoing := outgoing,
            lift := function(x)
                local representative;
                if not x in group then
                    Error("koAHSS: lift argument is not in the current cell");
                fi;
                representative := PreImagesRepresentative(quotient, x);
                if representative = fail then
                    Error("koAHSS: unable to lift a quotient element");
                fi;
                return previous.lift(representative);
            end,
            project := function(x)
                local representative;
                representative := previous.project(x);
                if not representative in cycles then
                    Error("koAHSS: E2 representative does not survive to E", r,
                        " at (p,q)=(", degree, ",", line, ")");
                fi;
                return Image(quotient, representative);
            end);
    end;
    getCell := function(r, degree, line)
        local id, group, previous, outgoing, incoming;
        id := key(r, degree, line);
        if IsBound(cells.(id)) then return cells.(id); fi;
        if degree < 0 or not line in [-4..0] or line = -3 then
            cells.(id) := emptyCell(r, degree, line);
        elif r = 2 then
            group := backend.cohomology(degree, line);
            validateGroup(group, degree, line);
            cells.(id) := makeBaseCell(group, degree, line);
        else
            previous := getCell(r - 1, degree, line);
            # A zero group remains zero on all later pages, even when an
            # adjacent group or operation is not yet determined.
            if not KOAHSS_IsUnresolved(previous) and IsTrivial(previous.group) then
                cells.(id) := ShallowCopy(previous);
                cells.(id).page := r;
                cells.(id).previous := previous;
            else
                outgoing := getMap(r - 1, degree, line);
                incoming := getMap(r - 1, degree - (r - 1), line + (r - 2));
                if KOAHSS_IsUnresolved(previous) or KOAHSS_IsUnresolved(outgoing)
                   or KOAHSS_IsUnresolved(incoming) then
                    cells.(id) := makeUnresolvedCell(r, degree, line,
                        previous, outgoing, incoming);
                else
                    cells.(id) := makeDerivedCell(r, degree, line,
                        previous, outgoing, incoming);
                fi;
            fi;
        fi;
        return cells.(id);
    end;
    getMap := function(r, degree, line)
        local id, source, target, map, gens, images, name;
        id := key(r, degree, line);
        if IsBound(maps.(id)) then return maps.(id); fi;
        # Symbolic zero maps do not require an unknown endpoint's group.
        # Check the arrow before materializing any unnecessary dependencies.
        if degree < 0 or not KOAHSS_ArrowExists(r, line) then
            maps.(id) := rec(status := "zero");
            return maps.(id);
        fi;
        source := getCell(r, degree, line);
        if not KOAHSS_IsUnresolved(source) and IsTrivial(source.group) then
            maps.(id) := rec(status := "zero");
            return maps.(id);
        fi;
        target := getCell(r, degree + r, line - r + 1);
        if not KOAHSS_IsUnresolved(target) and IsTrivial(target.group) then
            map := rec(status := "zero");
        elif KOAHSS_IsUnresolved(source) or KOAHSS_IsUnresolved(target) then
            map := mergeUnresolved([source, target]);
            AddSet(map.dependencies, [r, degree, line, KOAHSS_OperationName(r,line)]);
        else
            map := backend.differential(r, degree, line, source, target);
            if map = fail then
                name := KOAHSS_OperationName(r, line);
                map := rec(status := "unresolved", operations := [name],
                    reasons := [Concatenation("Unavailable ", name, " operation for d",
                        String(r), " at (p,q)=(", String(degree), ",", String(line), ")")]);
            fi;
            if KOAHSS_IsUnresolved(map) then
                map := mergeUnresolved([map]);
                AddSet(map.dependencies, [r, degree, line, KOAHSS_OperationName(r,line)]);
                maps.(id) := map;
                return map;
            fi;
            if not IsGeneralMapping(map) or not IsGroupHomomorphism(map) then
                Error("koAHSS: d", r, " at (p,q)=(", degree, ",", line,
                    ") must be a group homomorphism");
            fi;
            if not IsIdenticalObj(Source(map), source.group)
               or not IsIdenticalObj(Range(map), target.group) then
                Error("koAHSS: d", r, " at (p,q)=(", degree, ",", line,
                    ") has the wrong source or target group");
            fi;
            # Reconstruct with the checked constructor even if a backend used
            # an unchecked GAP constructor and merely asserted homomorphism.
            gens := GeneratorsOfGroup(source.group);
            images := List(gens, x -> Image(map, x));
            map := GroupHomomorphismByImages(source.group, target.group, gens, images);
            if map = fail then
                Error("koAHSS: d", r, " at (p,q)=(", degree, ",", line,
                    ") does not respect the source group relations");
            fi;
        fi;
        maps.(id) := map;
        return map;
    end;
    result := [];
    for page in pageNumbers do
        row := [];
        for q in [-4..0] do
            if bound-q-3<0 then Add(row,[]);
            else Add(row,List([0..bound-q-3],p->getCell(page,p,q))); fi;
        od;
        Add(result, row);
    od;
    # Keep the exact backend, hidden cells and quotient maps together.  In
    # particular, rebuilding a backend would change the meaning of twist
    # vectors and of representatives already lifted through these pages.
    return rec(kind := "koAHSSContext", space := space, backend := backend,
        maxDegree := bound, cohomologyCapacity := maxDegree,
        computedThrough := lastPage, pageNumbers := pageNumbers,
        pageData := result, getCell := getCell, getMap := getMap,
        cells := cells, maps := maps);
end);

BindGlobal("KOAHSS_ComputePageData", function(args)
    local context;
    context := KOAHSS_ComputePageContext(args);
    if Length(args) = 4 then return context.pageData[1]; fi;
    return context.pageData;
end);

InstallGlobalFunction(koAHSSPageData, function(arg)
    return KOAHSS_ComputePageData(arg);
end);

BindGlobal("KOAHSS_ClassifyPage", function(data)
    local classifyCell;
    classifyCell := function(cell)
        local result;
        if not KOAHSS_IsUnresolved(cell) then return AbelianInvariants(cell.group); fi;
        result := rec(status := "unresolved", page := cell.page,
            degree := cell.degree, q := cell.q, tertiaryCoefficient := cell.tertiaryCoefficient,
            operations := ShallowCopy(cell.operations),
            reasons := ShallowCopy(cell.reasons),
            dependencies := List(cell.dependencies, ShallowCopy));
        if IsBound(cell.tertiaryReference) then
            result.tertiaryReference:=cell.tertiaryReference;
        fi;
        if IsBound(cell.lastKnownInvariants) then
            result.lastKnownInvariants := ShallowCopy(cell.lastKnownInvariants);
            result.lastKnownPage := cell.lastKnownPage;
        fi;
        return result;
    end;
    return List(data, row -> List(row, classifyCell));
end);

InstallGlobalFunction(koAHSSpages, function(arg)
    local data;
    data := KOAHSS_ComputePageData(arg);
    if Length(arg) = 4 then return KOAHSS_ClassifyPage(data); fi;
    return List(data, KOAHSS_ClassifyPage);
end);
