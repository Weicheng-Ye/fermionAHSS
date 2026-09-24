# One resumable group job. Invoked by ../run_koahss_groups.py, not by hand.
if LoadPackage("json")=fail then Error("the GAP json package is required"); fi;
KOAHSSBatchJob:=JsonStringToGap(StringFile(GAPInfo.SystemEnvironment.KOAHSS_JOB_FILE));;
Read(Concatenation(KOAHSSBatchJob.root,"/load.g"));;
if LoadPackage("hap")=fail then Error("the GAP hap package is required"); fi;
Read(Concatenation(KOAHSSBatchJob.root,"/batch/koahss_groups.g"));;
Read(Concatenation(KOAHSSBatchJob.root,"/batch/koahss_spacegroup_twists.g"));;
Read(Concatenation(KOAHSSBatchJob.root,"/batch/koahss_paper_twists.g"));;
BreakOnError:=false;;
Reset(GlobalMersenneTwister,KOAHSSBatchJob.seed);;
Reset(GlobalRandomSource,KOAHSSBatchJob.seed);;

KOAHSSBatchEmit:=function(event)
    local stream;
    event.group_id:=KOAHSSBatchJob.spec.id;
    stream:=OutputTextFile(KOAHSSBatchJob.results,true);
    SetPrintFormattingStatus(stream,false);
    WriteAll(stream,Concatenation(GapToJsonString(event),"\n"));
    CloseStream(stream);
end;;

KOAHSSBatchCoordinateIndex:=function(coordinates)
    return Sum([1..Length(coordinates)],i->coordinates[i]*2^(i-1));
end;;

# Recover the character on this resolution from its values on group generators.
# Contracting g*e0 to e0 gives a 1-chain on which a 1-cocycle evaluates as s(g).
KOAHSSBatchOrientation:=function(model,backend,H)
    local R,generators,matrix,rhs,g,index,chain,i,solution,basis;
    R:=model.resolution; generators:=GeneratorsOfGroup(R!.group);
    matrix:=List(IdentityMat(backend.dimension(1)),v->backend.coboundary(1,v,false));
    rhs:=List([1..backend.dimension(2)],i->0);
    for g in generators do
        index:=Position(R!.elts,g);
        if index=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(g); else Add(R!.elts,g); fi;
            index:=Position(R!.elts,g);
        fi;
        if index=fail then Error("cannot index an orientation-character generator"); fi;
        chain:=R!.homotopy(0,[1,index]);
        for i in [1..Length(matrix)] do
            Add(matrix[i],Sum(Filtered(chain,t->AbsInt(t[1])=i),t->SignInt(t[1])) mod 2);
        od;
        Add(rhs,model.orientationCharacter(g));
    od;
    solution:=koAHSSSolveMod2System(matrix,rhs);
    if solution=fail then Error("orientation character has no cocycle representative"); fi;
    return rec(cochain:=solution.particular,coordinates:=Exponents(H.class(solution.particular)));
end;;

KOAHSSBatchSecondaryChecks:=[];;
KOAHSSBatchSecondary:=function(backend,n,a,options)
    local result,check;
    result:=koAHSSNaturalSecondary(backend,n,a,options);
    if result.status="candidate" then
        check:=rec(input_degree:=n,operation:=result.operation,
            reference_convention:=result.referenceConvention,
            normalization_coefficients:=result.normalizationCoefficients,
            identities_verified:=ForAll(result.equations,e->e.satisfied));
        if not check in KOAHSSBatchSecondaryChecks then Add(KOAHSSBatchSecondaryChecks,check); fi;
    fi;
    return result;
end;;

KOAHSSBatchSecondaryCallback:=function(name)
    return function(ctx)
        local options,result;
        options:=rec(inputType:="mod2");
        if name="Tau" then options.inputType:="integral"; fi;
        result:=KOAHSSBatchSecondary(ctx.backend,ctx.degree,ctx.cochain,options);
        if result.status="candidate" then return result.cochain; fi;
        return rec(status:="unresolved",operations:=[name],
            reasons:=[Concatenation("natural secondary formula: ",result.reason)]);
    end;
end;;

KOAHSSBatchTertiarySummary:=function(evaluation)
    local result;
    result:=rec(status:=evaluation.status);
    if IsBound(evaluation.reason) then result.reason:=evaluation.reason; fi;
    if evaluation.status="computed" then
        result.degree:=evaluation.inputDegree;
        result.modulus:=evaluation.modulus;
        result.phase_numerator:=evaluation.phaseNumerator;
        result.T:=evaluation.cochain;
        result.reference_convention:=evaluation.referenceConvention;
        result.local_R_choice:=false;
        result.old_correction_applied:=evaluation.oldCorrectionApplied;
        result.mu_R:=evaluation.muR;
        result.prime_three_coefficient:=evaluation.oddPrimaryCoefficient;
        result.bar_samples:=evaluation.barSamples;
        result.model_samples:=evaluation.modelSamples;
        result.transport_model:=evaluation.transportModel;
        result.input_group_bar_used:=evaluation.inputGroupBarUsed;
        result.carry_identity_verified:=evaluation.carryAudit.identityVerified;
        result.kernel_audit:=evaluation.kernelAudit;
    fi;
    return result;
end;;

KOAHSSBatchMain:=function()
    local job,model,R,base,H1,H2,r1,r2,basis1,basis2,orientation,header,
        sc,wc,s,omega,total,selected,attempted,skipped,errors,ops,evaluations,
        space,result,event,key,completed,start,depth,partial,coordinates,
        paperTwists,selectedTwists,twist,recovery,pageCount,recovered,distinctSelected;
    job:=KOAHSSBatchJob;
    # p+q+3 <= k bounds every displayed source. Outgoing differentials and
    # the matched secondary/tertiary checks need only this extra closure degree.
    depth:=Maximum(3,job.k+2);
    if job.pages>=3 then depth:=Maximum(3,job.k+3); fi;
    KOAHSSBatchEmit(rec(type:="group_start",resolution_length:=depth));
    model:=KOAHSSBatchMakeGroup(job.spec,depth); R:=model.resolution;
    base:=koAHSSHAPSpace(R).koAHSS(0,0,depth);
    H1:=base.data(1,-1); H2:=base.data(2,-1);
    basis1:=List(GeneratorsOfGroup(H1.group),H1.represent);
    basis2:=List(GeneratorsOfGroup(H2.group),H2.represent);
    if ForAny(Concatenation(GeneratorsOfGroup(H1.group),GeneratorsOfGroup(H2.group)),g->Order(g)<>2) then
        Error("twist cohomology must be elementary abelian over F2");
    fi;
    r1:=Length(basis1); r2:=Length(basis2); orientation:=fail;
    if model.kind="spacegroup" and IsBound(model.orientationCharacter) then
        orientation:=KOAHSSBatchOrientation(model,base,H1);
    fi;
    # Fix paper cocycles before taking the boundary snapshot: bar transport
    # may append group elements to the resolution's indexing table.
    paperTwists:=KOAHSSBatchPaperTwists(model,base,H1,H2,job.spec,orientation);
    for twist in paperTwists do
        twist.class_key:=Concatenation(String(KOAHSSBatchCoordinateIndex(twist.s_coordinates)),":",
            String(KOAHSSBatchCoordinateIndex(twist.omega_coordinates)));
        twist.key:=twist.class_key;
        if IsBound(twist.spacegroup_twist) then
            twist.key:=Concatenation(twist.key,":",twist.spacegroup_twist);
        fi;
    od;
    Sort(paperTwists,function(a,b)
        if a.class_key=b.class_key then
            return IsBound(a.spacegroup_twist) and IsBound(b.spacegroup_twist)
                and a.spacegroup_twist="zero" and b.spacegroup_twist<>"zero";
        fi;
        return [KOAHSSBatchCoordinateIndex(a.s_coordinates),KOAHSSBatchCoordinateIndex(a.omega_coordinates)]
             < [KOAHSSBatchCoordinateIndex(b.s_coordinates),KOAHSSBatchCoordinateIndex(b.omega_coordinates)];
    end);
    if IsBound(job.spacegroup_twist) then
        paperTwists:=Filtered(paperTwists,t->t.spacegroup_twist=job.spacegroup_twist);
    fi;
    if job.twists="paper" then selectedTwists:=paperTwists;
    elif job.twists="untwisted" then
        selectedTwists:=Filtered(paperTwists,t->ForAll(t.s_coordinates,x->x=0)
            and ForAll(t.omega_coordinates,x->x=0));
    elif job.twists="electronic" then
        if orientation=fail then Error("electronic twist selection requires a space group"); fi;
        selectedTwists:=Filtered(paperTwists,t->t.spacegroup_twist="zero");
    else Error("unknown paper twist filter"); fi;
    header:=rec(type:="basis",gap_version:=GAPInfo.Version,
        package_version:=KOAHSS_PACKAGE_VERSION,
        hap_version:=GAPInfo.PackagesInfo.hap[1].Version,
        group_name:=model.name,resolution_method:=model.resolutionMethod,
        dimensions:=List([0..depth],R!.dimension),
        h1_basis:=basis1,h2_basis:=basis2,h1_rank:=r1,h2_rank:=r2,
        total_twist_pairs:=2^(r1+r2),twist_policy:="paper-table-rows-and-spacegroup-pair-v1",
        paper_twists:=paperTwists,paper_twist_pairs:=Length(paperTwists),
        paper_distinct_twist_pairs:=Length(Set(List(paperTwists,t->t.class_key))),
        selected_twist_pairs:=Length(selectedTwists),
        selected_distinct_twist_pairs:=Length(Set(List(selectedTwists,t->t.class_key))),
        selected_twist_keys:=List(selectedTwists,t->t.key));
    # GAP's JSON writer does not serialize fail, despite reading null as fail.
    if orientation<>fail then header.orientation:=orientation; fi;
    if IsBound(model.productCoordinateConvention) then
        header.product_coordinate_convention:=model.productCoordinateConvention;
    fi;
    if IsBound(model.resolutionSource) then header.resolution_source:=model.resolutionSource; fi;
    if job.previous_basis<>fail and header<>job.previous_basis then
        Error("resolution/cohomology basis changed: use a fresh output directory");
    fi;
    # The complete boundary data (not just dimensions) prevent unsafe resume
    # across a change of high-degree resolution bases with the same H1/H2.
    coordinates:=rec(terms:=List([1..depth],n->List([1..R!.dimension(n)],i->R!.boundary(n,i))),
        elements:=List(R!.elts,String));
    if job.previous_boundary<>fail and coordinates<>job.previous_boundary then
        Error("resolution boundary changed: use a fresh output directory");
    fi;
    if job.previous_basis=fail then
        KOAHSSBatchEmit(header);
        KOAHSSBatchEmit(rec(type:="resolution_boundary",boundary:=coordinates));
    fi;
    total:=2^(r1+r2); completed:=Set(job.completed);
    attempted:=0; skipped:=0; errors:=0;
    selected:=Length(selectedTwists);
    if job.max_twists>0 then selected:=Minimum(selected,job.max_twists); fi;
    distinctSelected:=Length(Set(List(selectedTwists{[1..selected]},t->t.class_key)));
    for twist in selectedTwists{[1..selected]} do
        key:=twist.key;
        if key in completed then skipped:=skipped+1; continue; fi;
        sc:=twist.s_coordinates; wc:=twist.omega_coordinates;
        s:=twist.s_cochain; omega:=twist.omega_cochain;
        evaluations:=[]; ops:=rec(); KOAHSSBatchSecondaryChecks:=[];
        if job.operations="chosen" then
            ops.useNaturalPrimary:=true;
            ops.usesNaturalT:=true;
            ops.tertiaryReference:="Danus-degree0to3-chi7_tail-epsilon100-eta101-mu0";
            ops.Tau:=KOAHSSBatchSecondaryCallback("Tau");
            ops.Psi:=KOAHSSBatchSecondaryCallback("Psi");
            ops.T:=koAHSSNaturalTCallback(rec(evaluations:=evaluations));
        fi;
        # New adapter per twist keeps twist-specific caches separate.
        space:=koAHSSHAPSpace(R,ops);
        event:=rec(type:="case",key:=key,s_coordinates:=sc,omega_coordinates:=wc,
            s:=s,omega:=omega,operation_convention:=job.operations,
            paper_twist:=true,paper_source_rows:=ShallowCopy(twist.paper_source_rows),
            certified_ko:=false,k:=job.k,display_condition:="p+q+3<=k",
            page_numbers:=[2..job.pages+1],q_rows:=[-4..0]);
        if orientation<>fail then
            event.spacegroup_twist:=twist.spacegroup_twist;
            event.twist_class_key:=twist.class_key;
            # Paper comparison targets the zero-omega family, even if the
            # second family's pullback happens to give the same class.
            event.paper_electronic_case:=twist.spacegroup_twist="zero";
            if IsBound(twist.twist_audit) then event.twist_audit:=twist.twist_audit; fi;
        fi;
        start:=Runtime();
        KOAHSSBatchEmit(rec(type:="case_start",key:=key));
        # Persist every completed page before starting the next. Later failure
        # cannot turn missing E6 entries into zero or erase verified early pages.
        recovered:=[];
        for pageCount in [1..job.pages] do
            KOAHSSBatchEmit(rec(type:="page_start",key:=key,page:=pageCount+1));
            result:=CALL_WITH_CATCH(function()
                return koAHSSpages(space,s,omega,job.k,pageCount);
            end,[]);
            if not result[1] then break; fi;
            recovered:=result[2];
            KOAHSSBatchEmit(rec(type:="page_complete",key:=key,
                page:=pageCount+1,table:=recovered[pageCount],
                runtime_ms:=Runtime()-start));
        od;
        attempted:=attempted+1;
        event.tables:=recovered;
        event.page_numbers:=[2..Length(recovered)+1];
        event.status:="computed";
        if not result[1] then
            errors:=errors+1; event.status:="error";
            event.reason:="GAP rejected this calculation; see the group log";
            event.requested_page_numbers:=[2..job.pages+1];
            event.first_failed_page:=Length(recovered)+2;
        elif ForAny(recovered,t->ForAny(t,row->ForAny(row,IsRecord))) then
            event.status:="partial";
        fi;
        event.runtime_ms:=Runtime()-start;
        event.tertiary_evaluations:=List(evaluations,KOAHSSBatchTertiarySummary);
        if job.operations="chosen" then
            event.tertiary_convention:="Danus-degree0to3-chi7_tail-epsilon100-eta101-mu0";
            event.tertiary_coefficient:=0;
            event.tertiary_mu_R:=0;
            event.tertiary_prime_three_coefficient:=2;
        fi;
        event.secondary_formula_checks:=ShallowCopy(KOAHSSBatchSecondaryChecks);
        KOAHSSBatchEmit(event);
    od;
    KOAHSSBatchEmit(rec(type:="group_complete",selected_pairs:=selected,
        all_pairs:=total,attempted:=attempted,skipped:=skipped,errors:=errors,
        paper_pairs:=Length(paperTwists),selected_distinct_pairs:=distinctSelected,
        exhaustive_twists:=distinctSelected=total,
        exhaustive_paper_twists:=selected=Length(paperTwists)));
    return true;
end;;

KOAHSSBatchResult:=CALL_WITH_CATCH(KOAHSSBatchMain,[]);;
if not KOAHSSBatchResult[1] then
    KOAHSSBatchEmit(rec(type:="group_error",reason:="group setup or worker failed; see log"));
    FORCE_QUIT_GAP(2);
fi;
FORCE_QUIT_GAP(0);
