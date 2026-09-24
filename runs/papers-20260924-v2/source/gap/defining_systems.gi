# Construct an actual pair of nullhomotopies.  Vanishing in a secondary
# quotient may require changing b before the equation for c is soluble.
InstallGlobalFunction(koAHSSDefiningSystem,function(arg)
    local backend,n,A,options,secondaryOptions,evaluate,secondary,csol,
          H,K,generators,matrix,rhs,adjustment,h,i,again,opts,answer,
          suppliedC,suppliedBoundary,chooseC;
    if not Length(arg) in [3,4] then
        Error("koAHSSDefiningSystem(backend,n,A[,options])");
    fi;
    backend:=arg[1]; n:=arg[2]; A:=arg[3]; options:=rec();
    KOAHSS_EquationBackend(backend,n);
    KOAHSS_CC_CheckVector(A,backend.dimension(n),"defining-system input A");
    if ForAny(backend.coboundary(n,A,true),x->x<>0) then
        Error("koAHSS: defining-system A must be a sign-integral cocycle");
    fi;
    if Length(arg)=4 then options:=arg[4]; fi;
    if not IsRecord(options) then Error("koAHSS: defining-system options must be a record"); fi;
    secondaryOptions:=rec(inputType:="integral");
    if IsBound(options.secondaryOptions) then
        if not IsRecord(options.secondaryOptions) then
            Error("koAHSS: secondaryOptions must be a record");
        fi;
        secondaryOptions:=ShallowCopy(options.secondaryOptions);
        secondaryOptions.inputType:="integral";
    fi;
    if IsBound(options.b) then secondaryOptions.b:=options.b; fi;
    evaluate:=koAHSSSecondaryReference;
    if IsBound(backend.usesNaturalPrimary) and backend.usesNaturalPrimary then
        evaluate:=koAHSSNaturalSecondary;
    fi;
    if IsBound(options.secondaryEvaluator) then
        if not IsFunction(options.secondaryEvaluator) then
            Error("koAHSS: secondaryEvaluator must be a function");
        fi;
        evaluate:=options.secondaryEvaluator;
    fi;
    secondary:=evaluate(backend,n,A,secondaryOptions);
    if not IsRecord(secondary) or not IsBound(secondary.status) then
        Error("koAHSS: secondary evaluator must return a reference-result record");
    fi;
    if secondary.status<>"candidate" then return secondary; fi;
    if not IsBound(secondary.cochains) or not IsBound(secondary.cochains.psi)
       or not IsBound(secondary.cochains.b) then
        Error("koAHSS: integral secondary result needs cochains.psi and cochains.b");
    fi;
    h:=List([1..backend.dimension(n+1)],i->0);
    suppliedC:=fail; suppliedBoundary:=fail;
    if IsBound(options.c) then
        KOAHSS_CC_CheckVector(options.c,backend.dimension(n+2),"defining cochain c");
        suppliedC:=List(options.c,x->x mod 2);
        suppliedBoundary:=List(backend.coboundary(n+2,options.c,false),x->x mod 2);
    fi;
    chooseC:=function(rhs)
        if suppliedC<>fail and suppliedBoundary=rhs then
            # A checked supplied c proves existence directly. Do not invent
            # a homogeneous solution family when no elimination was needed.
            return rec(status:="solved",primitive:=ShallowCopy(suppliedC),
                particular:=ShallowCopy(suppliedC),primitiveDegree:=n+2,
                rhs:=ShallowCopy(rhs),modulus:=2,constraintsApplied:=0,
                choice:="validated supplied cochain",solutionFamilyComputed:=false);
        fi;
        return koAHSSSolveCochainEquation(backend,n+2,rhs);
    end;
    csol:=chooseC(secondary.cochains.psi);
    adjustment:=fail;
    if csol=fail then
        if not IsBound(backend.data) or not IsBound(backend.primary) then
            return rec(status:="unavailable",reason:="adjusting b requires cohomology representatives and D",
                       secondary:=secondary);
        fi;
        H:=backend.data(n+1,-1); K:=backend.data(n+3,-1);
        generators:=GeneratorsOfGroup(H.group);
        matrix:=List(generators,g->Exponents(K.class(backend.primary("D",n+1,H.represent(g)))));
        rhs:=Exponents(K.class(secondary.cochains.psi));
        adjustment:=koAHSSSolveMod2System(matrix,rhs);
        if adjustment=fail then
            return rec(status:="obstructed",equation:="[D h]=[psi(A,b)]",secondary:=secondary,
                       reason:="the input does not survive the secondary quotient");
        fi;
        for i in [1..Length(generators)] do
            if adjustment.particular[i]=1 then h:=h+H.represent(generators[i]); fi;
        od;
        h:=List(h,x->x mod 2);
        opts:=ShallowCopy(secondaryOptions);
        opts.b:=List(secondary.cochains.b+h,x->x mod 2);
        # Reuse the actual helpers; changing their choices here would change
        # the secondary reference while trying to vary only its nullhomotopy.
        if IsBound(secondary.helpers) then opts.helpers:=secondary.helpers; fi;
        again:=evaluate(backend,n,A,opts);
        if not IsRecord(again) or not IsBound(again.status) then
            Error("koAHSS: secondary evaluator returned an invalid adjusted result");
        fi;
        if again.status<>"candidate" then return again; fi;
        csol:=chooseC(again.cochains.psi);
        if csol=fail then
            return rec(status:="unsatisfied-identity",equation:="[psi(A,b+h)] = [psi(A,b)] + D[h]",
                       reason:="the chosen secondary helpers fail the required variation test",
                       secondary:=secondary,adjustedSecondary:=again,bAdjustment:=h);
        fi;
        secondary:=again;
    fi;
    if IsBound(options.c) then
        if suppliedBoundary<>secondary.cochains.psi then
            Error("koAHSS: supplied c does not satisfy dc=psi for the selected b");
        fi;
        csol.primitive:=ShallowCopy(suppliedC);
        csol.particular:=ShallowCopy(csol.primitive);
    fi;
    answer:=rec(status:="candidate",backend:=backend,n:=n,A:=ShallowCopy(A),
        b:=ShallowCopy(secondary.cochains.b),c:=ShallowCopy(csol.primitive),
        psi:=ShallowCopy(secondary.cochains.psi),secondary:=secondary,
        bAdjustment:=h,bAdjustmentSolutions:=adjustment,cSolutions:=csol,
        isNormalizedOperation:=false);
    return answer;
end);
