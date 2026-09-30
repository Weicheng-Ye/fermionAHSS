# Run from the package root: gap -q --quitonbreak examples/extension_papers.g
#
# Compares koFull_batch with the literature fixtures of
# data/extension-paper-samples.json (see doc/extension-paper-comparisons.md).
# Every case names a group or a resolution, the twist vectors in its basis,
# a cutoff, the compared quantity and the printed values by package degree.
# All cases are run by default (about five minutes for the finite groups);
# FERMIONAHSS_PAPER_CASES selects cases by a comma-separated list of id
# prefixes. The space groups are resolved by SGC_ResolutionSpaceGroup of the
# SpaceGroupCohomology package, whose resolutions carry the contracting
# homotopy that koFull needs (HAP's own space-group resolutions have none);
# the runner loads the package when it is installed. Cases whose group
# expression cannot be evaluated in this session are reported as skipped.
# The expected groups are read only after koFull_batch has completed.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash, prefix;
        file:=INPUT_FILENAME(); slash:=Positions(file,'/');
        prefix:="";
        if not IsEmpty(slash) then prefix:=file{[1..Last(slash)]}; fi;
        Read(Concatenation(prefix,"../load.g"));
    end,[]);
fi;
CallFuncList(function()
    local file, slash, prefix, fixture, selection, wanted, entry, results, counts,
        subject, evaluate, full, key, degree, expected, actual, comparison, stage,
        layer, name, table, cell, report;
    file:=INPUT_FILENAME(); slash:=Positions(file,'/'); prefix:="";
    if not IsEmpty(slash) then prefix:=file{[1..Last(slash)]}; fi;
    fixture:=JsonStringToGap(StringFile(Concatenation(prefix,"../data/extension-paper-samples.json")));
    selection:="all";
    if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_PAPER_CASES) then
        selection:=GAPInfo.SystemEnvironment.FERMIONAHSS_PAPER_CASES;
    fi;
    wanted:=SplitString(selection,",");
    # Evaluate a group expression; fail when the session cannot build it. The
    # error is caught without a break loop, so that --quitonbreak does not end
    # the session.
    evaluate:=function(expression)
        local attempt, oldBreak, oldSilent;
        oldBreak:=BreakOnError; oldSilent:=SilentNonInteractiveErrors;
        BreakOnError:=false; SilentNonInteractiveErrors:=true;
        attempt:=CALL_WITH_CATCH(EvalString,[expression]);
        BreakOnError:=oldBreak; SilentNonInteractiveErrors:=oldSilent;
        if attempt[1] then return attempt[2]; fi;
        return fail;
    end;
    # The resolutions of the space groups with their contracting homotopy; a
    # package that fails to load leaves those cases skipped.
    if not IsBoundGlobal("SGC_ResolutionSpaceGroup")
       and not IsEmpty(PackageInfo("SpaceGroupCohomology")) then
        evaluate("LoadPackage(\"SpaceGroupCohomology\",false)");
    fi;
    results:=[]; counts:=rec(match:=0,mismatch:=0,unresolved:=0,skipped:=0);
    for entry in fixture.cases do
        if selection<>"all" and not ForAny(wanted,w->StartsWith(entry.id,w)) then continue; fi;
        subject:=evaluate(entry.group_expression);
        if subject=fail then
            Add(results,rec(id:=entry.id,status:="skipped",
                reason:="the group expression cannot be evaluated in this session"));
            counts.skipped:=counts.skipped+1; continue;
        fi;
        Print("Computing ",entry.id," (",entry.fermionic_symmetry,")\n");
        full:=koFull_batch(subject,entry.s,entry.omega,entry.cutoff);
        for key in SortedList(RecNames(entry.expected_by_package_degree)) do
            degree:=Int(key); expected:=entry.expected_by_package_degree.(key);
            comparison:=rec(id:=entry.id,source:=entry.source,quantity:=entry.quantity,
                packageDegree:=degree,spatialDimension:=degree-1,expected:=expected,
                input:=rec(group_expression:=entry.group_expression,s:=entry.s,omega:=entry.omega));
            if entry.quantity="E6 layers" then
                table:=full.pages.tables[Position(full.pages.pageNumbers,6)];
                actual:=rec();
                for name in RecNames(expected) do
                    # The layer name ends with its bidegree (p,q).
                    cell:=SplitString(Filtered(name,c->c in "0123456789,-"),",");
                    actual.(name):=table[Int(cell[2])+5][Int(cell[1])+1];
                od;
                comparison.actual:=actual; comparison.status:="computed";
                comparison.match:=ForAll(RecNames(expected),function(n)
                    return expected.(n)=fail or SortedList(actual.(n))=SortedList(expected.(n)); end);
            else
                if full.degreeResults[degree+2].status<>"computed" then
                    comparison.status:="unresolved";
                    comparison.reason:=full.degreeResults[degree+2].reason;
                    counts.unresolved:=counts.unresolved+1; Add(results,comparison); continue;
                fi;
                comparison.status:="computed";
                if entry.quantity="full group" or degree>=4 then
                    actual:=full.degreeResults[degree+2].invariants;
                else
                    # The filtration stage below the p=0 layer: C at degree 1, B at degree 2, A at degree 3.
                    layer:=["D","C","B"][degree];
                    stage:=First(full.degreeResults[degree+2].filtration,f->f.layer=layer);
                    actual:=stage.invariants;
                fi;
                comparison.actual:=actual; comparison.match:=SortedList(actual)=SortedList(expected);
                comparison.relationMatrix:=full.degreeResults[degree+2].relationMatrix;
            fi;
            if comparison.match then counts.match:=counts.match+1; else counts.mismatch:=counts.mismatch+1; fi;
            Add(results,comparison);
        od;
    od;
    report:=rec(schemaVersion:=2,gapVersion:=GAPInfo.Version,selection:=selection,
        scope:="five-row-stacking-model",certified_ko:=false,counts:=counts,results:=results);
    WriteAll(OutputTextUser(),Concatenation("EXTENSION_PAPER_RESULTS ",GapToJsonString(report),"\n"));
    Print("matches ",counts.match,", mismatches ",counts.mismatch,", unresolved ",counts.unresolved,
        ", skipped ",counts.skipped,"\n");
    for comparison in Filtered(results,r->IsBound(r.match) and not r.match) do
        Print("  mismatch ",comparison.id," degree ",comparison.packageDegree,": paper ",comparison.expected,
            ", computed ",comparison.actual,"\n");
    od;
    # The Table VII groups of Wang and Gu are reproduced exactly.
    Assert(0,ForAll(Filtered(results,r->r.source="wang_gu_2020"),r->r.status="computed" and r.match));
end,[]);
