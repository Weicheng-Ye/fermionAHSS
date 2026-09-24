# Transport the exact lower compatibility M through supplied secondary
# coefficient changes. This is a cochain equation at one input, not a proof
# that the supplied coefficients or choices define a normalized stable family.
InstallGlobalFunction(koAHSSNormalizeSecondaryPair,
function(backend,n,A,reference,tauCoefficients,psiCoefficients)
    local ctx, result, checkIdentity, invalid, a, stored, name, replay,
          tauBasis, psiBasis, deltaPsi, deltaTau, i, numerator, betaDelta,
          rhs, matrix, basis, solution, v, generator, deltaM, M, psi, tau;
    if not IsInt(n) or n<0 then
        Error("koAHSS: secondary transport degree must be nonnegative");
    fi;
    for basis in [tauCoefficients,psiCoefficients] do
        if not IsList(basis) or Length(basis)<>3
           or not ForAll(basis,x -> IsInt(x) and x in [0,1]) then
            Error("koAHSS: secondary transport coefficients must be three binary integers");
        fi;
    od;
    ctx:=KOAHSS_SecondaryContext(backend);
    ctx.check(n,A);
    if ForAny(ctx.d(n,A,true),x -> x<>0) then
        Error("koAHSS: secondary transport input must be a sign-integral cocycle");
    fi;
    a:=ctx.cocycle(n,ctx.reduce(A));
    result:=rec(status:="candidate",inputDegree:=n,inputType:="integral",A:=ShallowCopy(A),
        reference:=reference,equations:=[],
        tauCoefficients:=ShallowCopy(tauCoefficients),
        psiCoefficients:=ShallowCopy(psiCoefficients),
        coefficientLabels:=rec(Tau:=["epsilon3","epsilon4","epsilon5"],
                               Psi:=["epsilon2","epsilon3","epsilon4"]),
        coefficientsSupplied:=true,isCoherentAtInput:=false,
        isNormalizedOperation:=false,isStableFamily:=false,
        scope:="coherent transport at this input and supplied coefficient pair");
    invalid:=function(stage,reason)
        result.status:="invalid-reference";
        result.failedStage:=stage; result.reason:=reason;
        return result;
    end;
    checkIdentity:=function(label,degree,lhs,rhs,modulus)
        if modulus=2 then lhs:=ctx.reduce(lhs); rhs:=ctx.reduce(rhs); fi;
        Add(result.equations,rec(name:=label,degree:=degree,lhs:=lhs,rhs:=rhs,
            coefficientModulus:=modulus,satisfied:=lhs=rhs));
        return lhs=rhs;
    end;
    if not IsRecord(reference) or not IsBound(reference.status)
       or reference.status<>"candidate" or not IsBound(reference.inputType)
       or reference.inputType<>"integral" or not IsBound(reference.inputDegree)
       or reference.inputDegree<>n or not IsBound(reference.cochains)
       or not IsRecord(reference.cochains) or not IsBound(reference.helpers)
       or not IsRecord(reference.helpers) then
        return invalid("reference metadata",
            "reference must be an integral koAHSSSecondaryReference candidate at this degree");
    fi;
    stored:=reference.cochains;
    if not ForAll(["A","a","s","omega","b","psi","tau","M"],
                  key -> IsBound(stored.(key))) then
        return invalid("reference cochains","the reference is missing its defining cochains");
    fi;
    if not IsBound(backend.twists) or not IsRecord(backend.twists)
       or not IsBound(backend.twists.s) or not IsBound(backend.twists.omega) then
        Error("koAHSS: secondary transport requires backend twists");
    fi;
    if stored.A<>A or stored.a<>a or stored.s<>backend.twists.s
       or stored.omega<>backend.twists.omega then
        return invalid("reference input",
            "the reference input or twist representatives differ from the supplied data");
    fi;
    ctx.check(n+1,stored.b); ctx.check(n+3,stored.psi);
    ctx.check(n+4,stored.tau); ctx.check(n+3,stored.M);
    if not ForAll(Concatenation(stored.b,stored.psi),x -> x in [0,1]) then
        return invalid("binary representatives","reference b and psi must be binary cochains");
    fi;
    if not checkIdentity("reference rho M = psi",n+3,stored.M,stored.psi,2)
       or not checkIdentity("reference d_s M = 2 tau",n+4,
                            ctx.d(n+3,stored.M,true),2*stored.tau,0)
       or not checkIdentity("reference d psi = 0",n+4,
                            ctx.d(n+3,stored.psi,false),ctx.zero(n+4),2)
       or not checkIdentity("reference d_s tau = 0",n+5,
                            ctx.d(n+4,stored.tau,true),ctx.zero(n+5),0) then
        return invalid("reference compatibility","the claimed lower compatibility does not hold");
    fi;
    # Recompute the actual reference formula using its chosen defining data.
    # Compatibility alone would allow a different secondary pair to be
    # relabelled as the reference operation.
    replay:=koAHSSSecondaryReference(backend,n,A,
        rec(b:=ShallowCopy(stored.b),helpers:=reference.helpers));
    result.referenceVerification:=replay;
    if replay.status<>"candidate" then
        return invalid("reference formula","the stored defining data do not reproduce a valid reference");
    fi;
    for name in ["psi","tau","M"] do
        if replay.cochains.(name)<>stored.(name) then
            return invalid("reference formula",
                Concatenation("the stored ",name," does not match its defining cochain formula"));
        fi;
    od;

    tauBasis:=koAHSSSecondaryCorrections(backend,n,A,"Tau");
    psiBasis:=koAHSSSecondaryCorrections(backend,n,a,"Psi");
    deltaPsi:=ctx.zero(n+3); deltaTau:=ctx.zero(n+4);
    for i in [1..3] do
        deltaPsi:=deltaPsi+tauCoefficients[i]*tauBasis[i];
        deltaTau:=deltaTau+psiCoefficients[i]*psiBasis[i];
    od;
    deltaPsi:=ctx.reduce(deltaPsi);
    result.deltaPsi:=deltaPsi; result.deltaTau:=deltaTau;
    result.correctionBases:=rec(Tau:=tauBasis,Psi:=psiBasis);
    numerator:=ctx.d(n+3,deltaPsi,true);
    if not checkIdentity("d_s lift deltaPsi is even",n+4,numerator,ctx.zero(n+4),2) then
        result.status:="unsatisfied-identity";
        result.reason:="the mod-two secondary correction has no integral Bockstein";
        return result;
    fi;
    betaDelta:=List(numerator,x -> QuoInt(x,2));
    rhs:=deltaTau-betaDelta;
    if not checkIdentity("transport rhs is closed",n+5,
                         ctx.d(n+4,rhs,true),ctx.zero(n+5),0) then
        result.status:="unsatisfied-identity";
        result.reason:="the integral compatibility obstruction is not closed";
        return result;
    fi;
    matrix:=[];
    for i in [1..ctx.dimension(n+3)] do
        basis:=ctx.zero(n+3); basis[i]:=1;
        Add(matrix,ctx.d(n+3,basis,true));
    od;
    result.transportEquation:=rec(
        name:="d_s v = deltaTau - d_s lift(deltaPsi)/2",
        primitiveDegree:=n+3,coefficientRing:=Integers,
        matrix:=matrix,rhs:=rhs,equationConvention:="v * matrix = rhs over Z");
    solution:=koAHSSSolveIntegerSystem(matrix,rhs);
    if solution=fail then
        result.status:="obstructed"; result.failedStage:="coherent secondary transport";
        result.reason:="the supplied coefficient pair has a nonzero integral compatibility obstruction";
        result.obstructionCocycle:=rhs;
        return result;
    fi;
    v:=solution.particular;
    if not checkIdentity("d_s v equals transport rhs",n+4,ctx.d(n+3,v,true),rhs,0) then
        Error("koAHSS: secondary transport solver does not match the backend coboundary");
    fi;
    for generator in solution.homogeneousGenerators do
        if ForAny(ctx.d(n+3,generator,true),x -> x<>0) then
            Error("koAHSS: a homogeneous secondary transport witness is not closed");
        fi;
    od;
    deltaM:=deltaPsi+2*v;
    M:=stored.M+deltaM; psi:=ctx.reduce(stored.psi+deltaPsi);
    tau:=stored.tau+deltaTau;
    if not checkIdentity("transported rho M = psi",n+3,M,psi,2)
       or not checkIdentity("transported d_s M = 2 tau",n+4,ctx.d(n+3,M,true),2*tau,0)
       or not checkIdentity("transported d psi = 0",n+4,ctx.d(n+3,psi,false),ctx.zero(n+4),2)
       or not checkIdentity("transported d_s tau = 0",n+5,ctx.d(n+4,tau,true),ctx.zero(n+5),0) then
        Error("koAHSS: solved secondary transport failed exact compatibility");
    fi;
    result.solutionLattice:=solution;
    # Preserve the common helper choices when this result is used by the
    # defining-system constructor and b subsequently has to be adjusted.
    result.helpers:=reference.helpers;
    result.v:=ShallowCopy(v); result.deltaM:=deltaM;
    result.homogeneousMCorrections:=List(solution.homogeneousGenerators,z -> 2*z);
    result.cochains:=rec(A:=ShallowCopy(A),a:=a,s:=ShallowCopy(stored.s),
        omega:=ShallowCopy(stored.omega),b:=ShallowCopy(stored.b),psi:=psi,tau:=tau,M:=M);
    result.operation:="Tau"; result.operationDegree:=n+3; result.cochain:=psi;
    result.isCoherentAtInput:=true;
    return result;
end);
