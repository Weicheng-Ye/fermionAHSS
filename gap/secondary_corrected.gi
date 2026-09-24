# The revised Psi note calibrates a particular degree-nine universal helper
# and its stable continuation.  Declaring that convention is required here;
# solving a primitive on one resolution does not establish it.
BindGlobal("koAHSSPsiNoteConvention",function()
    return rec(id:="chi9_C1_stable",referenceDegree:=9,
        helperSimplexCode:=87112285934295547847080502420714144399368,
        helperCoordinateIndicesZeroBased:=[83,118,162,216],helperValue:=1,
        coefficientLabels:=["epsilon3","epsilon4","epsilon5"],
        coefficients:=[0,1,1],operation:="Tau",noteOperation:="Psi",
        requirement:="the chi9(C)=1 universal reference and its stable continuation");
end);

# This primary cochain can be computed independently of any assertion about
# the secondary reference.  Use the note's actual cup-one representative.
BindGlobal("koAHSSPsiNoteCorrection",function(backend,n,input)
    local ctx,s,omega,a,somega,sq1omega,kappa;
    if not IsInt(n) or n<0 then
        Error("koAHSS: Psi-note correction degree must be nonnegative");
    fi;
    ctx:=KOAHSS_SecondaryContext(backend);
    if not IsBound(backend.twists) or not IsRecord(backend.twists)
       or not IsBound(backend.twists.s) or not IsBound(backend.twists.omega) then
        Error("koAHSS: Psi-note correction requires backend twists");
    fi;
    s:=ctx.cocycle(1,backend.twists.s);
    omega:=ctx.cocycle(2,backend.twists.omega);
    a:=ctx.cocycle(n,input);
    somega:=ctx.cup(0,1,s,2,omega);
    sq1omega:=ctx.cup(1,2,omega,2,omega);
    kappa:=ctx.reduce(ctx.cup(0,3,somega,n,a)
                    +ctx.cup(0,3,sq1omega,n,a));
    if ForAny(ctx.d(n+3,kappa,false),x->x mod 2<>0) then
        Error("koAHSS: the Psi-note correction is not a mod-two cocycle");
    fi;
    return kappa;
end);

BindGlobal("koAHSSSecondaryCorrected",function(arg)
    local backend,n,input,options,inputType,convention,reference,rawOptions,
          replay,ctx,required,name,result,kappa,numerator,deltaTau,
          cochains,checkIdentity,invalid;
    if not Length(arg) in [3,4] then
        Error("koAHSSSecondaryCorrected(backend,n,input[,options])");
    fi;
    backend:=arg[1]; n:=arg[2]; input:=arg[3]; options:=rec();
    if Length(arg)=4 then options:=arg[4]; fi;
    if not IsRecord(options) then Error("koAHSS: corrected secondary options must be a record"); fi;
    if not IsInt(n) or n<0 then Error("koAHSS: corrected secondary degree must be nonnegative"); fi;
    inputType:="integral";
    if IsBound(options.inputType) then inputType:=options.inputType; fi;
    if not inputType in ["integral","mod2"] then
        Error("koAHSS: inputType must be integral or mod2");
    fi;
    convention:=fail; reference:=fail;
    if IsBound(options.reference) then reference:=options.reference; fi;
    if IsRecord(reference) and IsBound(reference.referenceConvention) then
        convention:=reference.referenceConvention;
    fi;
    if IsBound(options.referenceConvention) then
        convention:=options.referenceConvention;
    fi;
    required:=koAHSSPsiNoteConvention();
    if convention<>required.id then
        return rec(status:="unavailable",failedStage:="reference convention",
            reason:="the Psi-note correction requires an explicit chi9_C1_stable reference declaration",
            requiredReferenceConvention:=required,isNormalizedOperation:=false);
    fi;
    ctx:=KOAHSS_SecondaryContext(backend);
    ctx.check(n,input);
    if inputType="integral" and ForAny(ctx.d(n,input,true),x->x<>0) then
        Error("koAHSS: corrected integral secondary input must be a sign-integral cocycle");
    fi;
    rawOptions:=ShallowCopy(options);
    rawOptions.inputType:=inputType;
    if IsBound(rawOptions.reference) then Unbind(rawOptions.reference); fi;
    if IsBound(rawOptions.referenceConvention) then Unbind(rawOptions.referenceConvention); fi;
    invalid:=function(reason)
        return rec(status:="invalid-reference",failedStage:="reference formula",
            reason:=reason,isNormalizedOperation:=false);
    end;
    if reference=fail then
        reference:=koAHSSSecondaryReference(backend,n,input,rawOptions);
        if reference.status<>"candidate" then return reference; fi;
    else
        if not IsRecord(reference) or not IsBound(reference.status)
           or reference.status<>"candidate" or not IsBound(reference.inputDegree)
           or reference.inputDegree<>n or not IsBound(reference.inputType)
           or reference.inputType<>inputType or not IsBound(reference.cochains)
           or not IsRecord(reference.cochains) or not IsBound(reference.helpers)
           or not IsRecord(reference.helpers) or not IsBound(reference.equations)
           or not IsList(reference.equations) then
            return invalid("expected a raw secondary reference at the supplied degree and input type");
        fi;
        if IsBound(reference.isPsiNoteCorrected) and reference.isPsiNoteCorrected=true then
            return invalid("the Psi-note correction has already been applied");
        fi;
        if not IsBound(reference.cochains.b) then return invalid("the reference is missing b"); fi;
        if IsBound(rawOptions.b) and rawOptions.b<>reference.cochains.b then
            return invalid("the supplied b differs from the stored reference");
        fi;
        rawOptions.b:=ShallowCopy(reference.cochains.b);
        rawOptions.helpers:=reference.helpers;
        replay:=koAHSSSecondaryReference(backend,n,input,rawOptions);
        if replay.status<>"candidate" then
            return invalid("the stored helpers do not reproduce a valid raw reference");
        fi;
        for name in RecNames(replay.cochains) do
            if not IsBound(reference.cochains.(name))
               or reference.cochains.(name)<>replay.cochains.(name) then
                return invalid(Concatenation("the stored ",name," differs from the raw reference formula"));
            fi;
        od;
        if not IsBound(reference.cochain) or reference.cochain<>replay.cochain then
            return invalid("the stored output differs from the raw reference formula");
        fi;
        for name in ["operation","noteOperation","operationDegree"] do
            if not IsBound(reference.(name)) or reference.(name)<>replay.(name) then
                return invalid("the stored operation metadata differs from the raw reference");
            fi;
        od;
    fi;
    kappa:=koAHSSPsiNoteCorrection(backend,n,input);
    numerator:=ctx.d(n+3,kappa,true);
    if ForAny(numerator,x->x mod 2<>0) then
        Error("koAHSS: the corrected secondary Bockstein numerator is not even");
    fi;
    deltaTau:=List(numerator,x->QuoInt(x,2));
    result:=ShallowCopy(reference);
    cochains:=StructuralCopy(reference.cochains);
    result.cochains:=cochains;
    result.reference:=reference;
    result.equations:=ShallowCopy(reference.equations);
    result.referenceConvention:=required.id;
    result.normalizationDeclaration:=required;
    result.referenceConventionVerified:=false;
    result.isPsiNoteCorrected:=true;
    result.isNormalizedOperation:=false;
    result.isStableFamily:=false;
    result.isLocalCandidate:=true;
    result.generalPsiNormalized:=false;
    result.correction:=ShallowCopy(kappa);
    result.deltaTau:=ShallowCopy(deltaTau);
    result.scope:="exact cochain correction conditional on the caller's declared reference convention";
    cochains.kappa:=ShallowCopy(kappa);
    cochains.tau:=cochains.tau+deltaTau;
    checkIdentity:=function(label,degree,lhs,rhs,modulus)
        if modulus=2 then lhs:=ctx.reduce(lhs); rhs:=ctx.reduce(rhs); fi;
        Add(result.equations,rec(name:=label,degree:=degree,lhs:=lhs,rhs:=rhs,
            coefficientModulus:=modulus,satisfied:=lhs=rhs));
        if lhs<>rhs then Error("koAHSS: corrected secondary identity failed: ",label); fi;
    end;
    checkIdentity("d_s corrected tau = 0",n+5,
        ctx.d(n+4,cochains.tau,true),ctx.zero(n+5),0);
    if inputType="integral" then
        cochains.psi:=ctx.reduce(cochains.psi+kappa);
        cochains.M:=cochains.M+kappa;
        result.cochain:=ShallowCopy(cochains.psi);
        result.comparisonScope:="matched Bockstein comparison on sign-integral input";
        checkIdentity("rho corrected M = corrected psi",n+3,cochains.M,cochains.psi,2);
        checkIdentity("d_s corrected M = 2 corrected tau",n+4,
            ctx.d(n+3,cochains.M,true),2*cochains.tau,0);
        checkIdentity("d corrected psi = 0",n+4,
            ctx.d(n+3,cochains.psi,false),ctx.zero(n+4),2);
    else
        result.cochain:=ShallowCopy(cochains.tau);
        result.comparisonScope:="formal local beta_s(kappa) extension on mod-two input; general Theta normalization is not established";
    fi;
    return result;
end);
