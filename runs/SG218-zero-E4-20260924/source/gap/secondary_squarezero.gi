# The revised Psi note proves J Psi=0, J=beta_s(Sq2+omega). When
# J is injective on H^(n+3)(F2)/D H^(n+1)(F2), this determines Psi's
# coset uniquely. This test does not choose a primitive to hide a nonzero
# square, and makes no claim to construct the degree-nine universal helper.
BindGlobal("koAHSSSecondarySquareZeroTarget",function(backend,n)
    local source,target,outgoing,images,denominator,projection,quotient,
          generators,jImages,map,quotientMap,g;
    KOAHSS_EquationBackend(backend,n);
    if not IsBound(backend.data) or not IsBound(backend.primary) then
        return rec(status:="unavailable",reason:="cohomology and primary operations are required");
    fi;
    source:=backend.data(n+1,-1); target:=backend.data(n+3,-2);
    outgoing:=backend.data(n+6,-4);
    images:=List(GeneratorsOfGroup(source.group),g->
        target.class(backend.primary("D",n+1,source.represent(g))));
    denominator:=Subgroup(target.group,images);
    projection:=NaturalHomomorphismByNormalSubgroup(target.group,denominator);
    quotient:=Image(projection);
    generators:=GeneratorsOfGroup(target.group);
    jImages:=List(generators,g->outgoing.class(
        backend.primary("Dtilde",n+3,target.represent(g))));
    map:=GroupHomomorphismByImages(target.group,outgoing.group,generators,jImages);
    if map=fail then Error("koAHSS: outgoing J is not a cohomology homomorphism"); fi;
    if not ForAll(GeneratorsOfGroup(denominator),g->Image(map,g)=One(outgoing.group)) then
        return rec(status:="unsatisfied-identity",reason:="J D is nonzero on cohomology");
    fi;
    generators:=GeneratorsOfGroup(quotient);
    quotientMap:=GroupHomomorphismByImages(quotient,outgoing.group,generators,
        List(generators,g->Image(map,PreImagesRepresentative(projection,g))));
    if quotientMap=fail then Error("koAHSS: J does not descend to the secondary quotient"); fi;
    return rec(status:="computed",inputDegree:=n,
        quotient:=quotient,projection:=projection,outgoingMap:=quotientMap,
        quotientInvariants:=AbelianInvariants(quotient),
        kernelInvariants:=AbelianInvariants(Kernel(quotientMap)),
        determinesZero:=IsTrivial(Kernel(quotientMap)),
        reason:="the revised Psi note proves J Psi=0; an injective J determines the zero coset");
end);

# A matched, local cochain realization of a theorem-determined zero coset.
# The canonical integral lift and b are selected from a, before evaluating
# the requested A or b. Thus later b changes reuse the same chi correction.
# This is not transport of the universal degree-nine primitive. If the
# theorem does not determine a unique class, keep the raw local convention
# and say so explicitly instead of forcing an arbitrary class into ker J.
BindGlobal("koAHSSSecondaryWithSquareZero",function(arg)
    local backend,n,input,options,inputType,ctx,a,proof,matrix,i,v,rhs,
          lift,canonicalA,canonicalOptions,canonical,helpers,chiShift,
          selectedOptions,raw,result,kappa,deltaTau,cochains,checkIdentity;
    if not Length(arg) in [3,4] then
        Error("koAHSSSecondaryWithSquareZero(backend,n,input[,options])");
    fi;
    backend:=arg[1]; n:=arg[2]; input:=arg[3]; options:=rec();
    if Length(arg)=4 then options:=arg[4]; fi;
    if not IsRecord(options) then Error("koAHSS: secondary options must be a record"); fi;
    if IsBound(options.referenceConvention) then
        return koAHSSSecondaryCorrected(backend,n,input,options);
    fi;
    inputType:="integral";
    if IsBound(options.inputType) then inputType:=options.inputType; fi;
    if not inputType in ["integral","mod2"] then Error("koAHSS: invalid secondary inputType"); fi;
    ctx:=KOAHSS_SecondaryContext(backend); a:=ctx.cocycle(n,input);
    if inputType="integral" and ForAny(ctx.d(n,input,true),x->x<>0) then
        Error("koAHSS: integral secondary input must be a sign-integral cocycle");
    fi;
    # Cache only the input-independent cohomology proof, within this backend
    # and hence these particular twists and resolution conventions.
    if not IsBound(backend.psiNoteSquareZeroTargets) then backend.psiNoteSquareZeroTargets:=[]; fi;
    if not IsBound(backend.psiNoteSquareZeroTargets[n+1]) then
        backend.psiNoteSquareZeroTargets[n+1]:=koAHSSSecondarySquareZeroTarget(backend,n);
    fi;
    proof:=backend.psiNoteSquareZeroTargets[n+1];
    if proof.status<>"computed" then return proof; fi;
    if not proof.determinesZero then
        result:=koAHSSSecondaryReference(backend,n,input,options);
        result.squareZeroTarget:=proof;
        result.noteNormalizationStatus:="not determined by injectivity; raw local convention retained";
        return result;
    fi;
    # A mod-two input need not have a sign-integral lift. The new note's
    # Bockstein comparison does not normalize Theta on those inputs.
    rhs:=ctx.d(n,a,true);
    if ForAny(rhs,x->x mod 2<>0) then Error("koAHSS: cocycle has odd Bockstein numerator"); fi;
    rhs:=List(rhs,x->-QuoInt(x,2)); matrix:=[];
    for i in [1..ctx.dimension(n)] do
        v:=ctx.zero(n); v[i]:=1; Add(matrix,ctx.d(n,v,true));
    od;
    lift:=koAHSSSolveIntegerSystem(matrix,rhs);
    if lift=fail then
        if inputType="integral" then Error("koAHSS: supplied integral lift was not recovered"); fi;
        result:=koAHSSSecondaryReference(backend,n,input,options);
        result.noteNormalizationStatus:="mod-two input has no sign-integral lift; raw local convention retained";
        return result;
    fi;
    canonicalA:=a+2*lift.particular;
    canonicalOptions:=ShallowCopy(options); canonicalOptions.inputType:="integral";
    if IsBound(canonicalOptions.b) then Unbind(canonicalOptions.b); fi;
    if IsBound(canonicalOptions.helpers) then
        canonicalOptions.helpers:=ShallowCopy(canonicalOptions.helpers);
        if IsBound(canonicalOptions.helpers.b) then Unbind(canonicalOptions.helpers.b); fi;
    fi;
    canonical:=koAHSSSecondaryReference(backend,n,canonicalA,canonicalOptions);
    if canonical.status<>"candidate" then return canonical; fi;
    kappa:=koAHSSPsiNoteCorrection(backend,n,a);
    chiShift:=ctx.reduce(canonical.cochains.psi+kappa);
    helpers:=ShallowCopy(canonical.helpers);
    helpers.chi:=ShallowCopy(helpers.chi);
    helpers.chi.selectedPrimitive:=ctx.reduce(helpers.chi.selectedPrimitive+chiShift);
    if IsBound(options.helpers) and IsBound(options.helpers.b) then
        helpers.b:=options.helpers.b;
    fi;
    selectedOptions:=ShallowCopy(options); selectedOptions.helpers:=helpers;
    raw:=koAHSSSecondaryReference(backend,n,input,selectedOptions);
    if raw.status<>"candidate" then return raw; fi;
    deltaTau:=ctx.d(n+3,kappa,true);
    if ForAny(deltaTau,x->x mod 2<>0) then Error("koAHSS: kappa Bockstein is not integral"); fi;
    deltaTau:=List(deltaTau,x->QuoInt(x,2));
    result:=ShallowCopy(raw); result.reference:=raw;
    result.cochains:=StructuralCopy(raw.cochains); cochains:=result.cochains;
    result.equations:=ShallowCopy(raw.equations);
    cochains.kappa:=kappa; cochains.tau:=cochains.tau+deltaTau;
    result.squareZeroTarget:=proof; result.canonicalIntegralLift:=canonicalA;
    result.chiAdjustment:=chiShift; result.correction:=kappa;
    result.noteNormalizationStatus:="zero Tau coset determined by injective outgoing J";
    result.isUniversalHelper:=false; result.isNormalizedOperation:=false;
    result.scope:="theorem-determined Tau coset with a matched local cochain realization";
    checkIdentity:=function(label,degree,lhs,rhs,modulus)
        if modulus=2 then lhs:=ctx.reduce(lhs); rhs:=ctx.reduce(rhs); fi;
        Add(result.equations,rec(name:=label,degree:=degree,lhs:=lhs,rhs:=rhs,
            coefficientModulus:=modulus,satisfied:=lhs=rhs));
        if lhs<>rhs then Error("koAHSS: square-zero secondary realization failed: ",label); fi;
    end;
    checkIdentity("d_s matched tau = 0",n+5,ctx.d(n+4,cochains.tau,true),ctx.zero(n+5),0);
    if inputType="integral" then
        cochains.psi:=ctx.reduce(cochains.psi+kappa); cochains.M:=cochains.M+kappa;
        result.cochain:=cochains.psi;
        checkIdentity("rho matched M = psi",n+3,cochains.M,cochains.psi,2);
        checkIdentity("d_s matched M = 2 tau",n+4,ctx.d(n+3,cochains.M,true),2*cochains.tau,0);
        if Image(proof.projection,backend.data(n+3,-2).class(cochains.psi))<>One(proof.quotient) then
            return rec(status:="unsatisfied-identity",reason:="the local cochain realization does not represent the theorem-determined zero coset");
        fi;
    else result.cochain:=cochains.tau; fi;
    return result;
end);
