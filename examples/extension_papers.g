# Run from the package root: gap -q --quitonbreak examples/extension_papers.g
# The expected groups are used only after koFull has completed independently.
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
    local fixture, cache, rows, case, full, degree, result, row, matches,
        unresolved, mismatches, key, report, file, slash, prefix, name, j, vector,
        extensionOptions,comparison,powerWitness;
    file:=INPUT_FILENAME(); slash:=Positions(file,'/'); prefix:="";
    if not IsEmpty(slash) then prefix:=file{[1..Last(slash)]}; fi;
    fixture:=JsonStringToGap(StringFile(Concatenation(
        prefix,"../data/extension-paper-samples.json")));
    cache:=rec(); rows:=[]; matches:=0; unresolved:=0; mismatches:=0;
    extensionOptions:=rec();
    if IsBoundGlobal("FERMIONAHSS_EXTENSION_MODEL") then
        extensionOptions.extensionModel:=ValueGlobal("FERMIONAHSS_EXTENSION_MODEL");
    fi;
    for case in fixture.cases do
        if IsBound(case.reuse_case) then full:=cache.(case.reuse_case);
        else
            full:=koFull(CyclicGroup(2),case.s,case.omega,case.cutoff,extensionOptions);
            cache.(case.id):=full;
            Print("Completed ",case.id,": ",full.status,"\n");
        fi;
        for key in SortedList(RecNames(case.expected_by_package_degree)) do
            degree:=Int(key); result:=full.degreeResults[degree+2];
            row:=rec(caseId:=case.id,source:=case.source,packageDegree:=degree,
                spatialDimension:=degree-1,expected:=case.expected_by_package_degree.(key),
                input:=rec(group:="CyclicGroup(2)",s:=case.s,omega:=case.omega,cutoff:=case.cutoff),
                status:=result.status);
            if IsBound(result.modelSelection) then row.modelSelection:=result.modelSelection; fi;
            if IsBound(result.certificateLevel) then row.certificateLevel:=result.certificateLevel; fi;
            if result.status="computed" then
                row.actual:=result.invariants; row.relationMatrix:=result.relationMatrix;
                row.match:=row.actual=row.expected;
                if IsBound(result.algebraAudit) then
                    row.finiteAudit:=rec(normalFormCount:=Length(result.algebraAudit.table),
                        table:=result.algebraAudit.table,
                        associativityVerified:=result.algebraAudit.associativityVerified,
                        commutativityVerified:=result.algebraAudit.commutativityVerified,
                        scope:=result.algebraAudit.scope);
                    row.fixedFlatLifts:=[]; row.powerWitnesses:=[];
                    for name in ["D","C","B","A"] do
                        for j in [1..Length(result.layers.(name).fullLifts)] do
                            Add(row.fixedFlatLifts,rec(layer:=name,generator:=j,
                                state:=result.layers.(name).fullLifts[j].state));
                        od;
                    od;
                    for vector in result.extensionVectors do
                        if IsBound(vector.result.witness.stackedState) then
                            comparison:=vector.result.witness.reduction.canonicalComparison;
                            powerWitness:=rec(layer:=vector.layer,order:=vector.order,
                                lowerCoordinates:=vector.result.lowerCoordinates,
                                stackedState:=vector.result.witness.stackedState,
                                canonicalLowerProduct:=vector.result.witness.reduction.canonicalLowerProduct,
                                gauge:=comparison.gauge,equalityVerified:=comparison.equalityVerified);
                            if IsBound(comparison.boundary) then powerWitness.boundary:=comparison.boundary; fi;
                            if IsBound(comparison.certificateLevel) then
                                powerWitness.certificateLevel:=comparison.certificateLevel;
                            fi;
                            Add(row.powerWitnesses,powerWitness);
                        fi;
                    od;
                fi;
                if row.match then matches:=matches+1; else mismatches:=mismatches+1; fi;
            else
                row.reason:=result.reason; row.pendingLayer:=result.pendingLayer;
                unresolved:=unresolved+1;
            fi;
            Add(rows,row);
        od;
    od;
    report:=rec(schemaVersion:=1,gapVersion:=GAPInfo.Version,
        scope:="five-row-stacking-model",certified_ko:=false,
        extensionModel:=full.extensionModel,
        distinctCalculations:=Length(RecNames(cache)),matches:=matches,
        unresolved:=unresolved,mismatches:=mismatches,results:=rows);
    WriteAll(OutputTextUser(),Concatenation("EXTENSION_PAPER_RESULTS ",GapToJsonString(report),"\n"));
    Assert(0,mismatches=0);
    if IsBound(extensionOptions.extensionModel) and extensionOptions.extensionModel="transfer" then
        Assert(0,unresolved=0 and matches=Length(rows));
        Assert(0,ForAll(Filtered(rows,row->row.packageDegree in [3..5]),row->
            IsBound(row.modelSelection) and row.modelSelection.selected="transfer" and
            not row.modelSelection.fallback));
    fi;
end,[]);
