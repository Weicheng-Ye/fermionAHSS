# Solve the missing tertiary phase equation on one cochain complex.
# These are LOCAL candidates, not a normalized natural cohomology operation.
InstallGlobalFunction(koAHSSTertiaryCandidates,
function(arg)
    local backend, n, psi, c, modulus, dimension, checkVector,
          differential, cup, reduce, omega, s,
          upperDegree, outputDegree, upperDimension, outputDimension,
          dc, q2c, upperBinary, obstruction, half, upperPhase,
          upperDifferential, rhs, matrix, basis, row, i, solutions,
          correctionDifferential, numerator, particularCocycle,
          homogeneousCocycles, generator, cocycle, multiply,
          correction, coefficient, unshiftedParticularCocycle;
    if not Length(arg) in [5,6] then
        Error("koAHSSTertiaryCandidates(backend, n, psi, c, modulus[, A])");
    fi;
    backend := arg[1]; n := arg[2]; psi := arg[3];
    c := arg[4]; modulus := arg[5];
    if not IsInt(n) or n < 0 then
        Error("koAHSS: tertiary input degree must be a nonnegative integer");
    fi;
    if not IsInt(modulus) or modulus < 2 or modulus mod 2 <> 0 then
        Error("koAHSS: tertiary phase modulus must be a positive even integer");
    fi;
    if not IsRecord(backend) or not IsBound(backend.dimension)
       or not IsFunction(backend.dimension)
       or not IsBound(backend.coboundary)
       or not IsFunction(backend.coboundary)
       or not IsBound(backend.cupMod2)
       or not IsFunction(backend.cupMod2)
       or not IsBound(backend.twists) or not IsRecord(backend.twists)
       or not IsBound(backend.twists.s)
       or not IsBound(backend.twists.omega) then
        Error("koAHSS: tertiary candidates require a cochain backend with twists");
    fi;
    dimension := function(degree)
        local value;
        value := backend.dimension(degree);
        if not IsInt(value) or value < 0 then
            Error("koAHSS: cochain dimensions must be nonnegative integers");
        fi;
        return value;
    end;
    checkVector := function(vector, degree, binary)
        if not IsList(vector) or Length(vector) <> dimension(degree)
           or not ForAll(vector, IsInt) then
            Error("koAHSS: tertiary cochain has incorrect dimension or noninteger entries");
        fi;
        if binary and not ForAll(vector, x -> x in [0,1]) then
            Error("koAHSS: psi, c, and twists must use binary representatives");
        fi;
    end;
    differential := function(degree, vector, sign)
        local value;
        checkVector(vector, degree, false);
        value := backend.coboundary(degree, vector, sign);
        checkVector(value, degree+1, false);
        return value;
    end;
    reduce := vector -> List(vector, x -> x mod 2);
    cup := function(index, degreeA, a, degreeB, b)
        local value;
        value := backend.cupMod2(index, degreeA, a, degreeB, b);
        checkVector(value, degreeA+degreeB-index, false);
        return reduce(value);
    end;
    s := backend.twists.s;
    omega := backend.twists.omega;
    checkVector(s, 1, true);
    checkVector(omega, 2, true);
    checkVector(c, n+2, true);
    checkVector(psi, n+3, true);
    if not ForAll(differential(1,s,false), x -> x mod 2 = 0)
       or not ForAll(differential(2,omega,false), x -> x mod 2 = 0) then
        Error("koAHSS: tertiary twists must be mod-two cocycles");
    fi;
    if not ForAll(differential(n+3,psi,false), x -> x mod 2 = 0) then
        Error("koAHSS: psi must be a mod-two cocycle");
    fi;
    dc := reduce(differential(n+2,c,false));
    if dc <> psi then
        Error("koAHSS: the defining cochain c must satisfy dc = psi modulo two");
    fi;

    upperDegree := n+4;
    outputDegree := n+5;
    upperDimension := dimension(upperDegree);
    outputDimension := dimension(outputDegree);
    half := QuoInt(modulus,2);
    q2c := reduce(cup(n,n+2,c,n+2,c) + cup(n+1,n+2,c,n+3,dc));
    upperBinary := reduce(q2c + cup(0,2,omega,n+2,c));
    obstruction := reduce(cup(n+1,n+3,psi,n+3,psi)
                           + cup(0,2,omega,n+3,psi));
    upperPhase := List(upperBinary, x -> half*x);
    upperDifferential := differential(upperDegree,upperPhase,true);
    if List(upperDifferential, x -> x mod modulus)
       <> List(obstruction, x -> half*x) then
        Error("koAHSS: cup operations fail the required upper phase identity");
    fi;
    rhs := List(obstruction, x -> (-half*x) mod modulus);

    # Row-vector convention: a cochain x has coboundary x * matrix.
    matrix := [];
    for i in [1..upperDimension] do
        basis := List([1..upperDimension], j -> 0);
        basis[i] := 1;
        row := differential(upperDegree,basis,true);
        if not ForAll(differential(outputDegree,row,true), x -> x = 0) then
            Error("koAHSS: the twisted coboundary does not square to zero");
        fi;
        Add(matrix,row);
    od;
    multiply := function(vector)
        if upperDimension = 0 or outputDimension = 0 then
            return List([1..outputDimension], j -> 0);
        fi;
        return vector * matrix;
    end;
    if upperDifferential <> multiply(upperPhase) then
        Error("koAHSS: coboundary disagrees with its integer matrix");
    fi;
    solutions := koAHSSSolvePhaseSystem(matrix,rhs,modulus);
    if solutions = fail then
        # Given the upper identity, -upperPhase is already a solution.
        Error("koAHSS: phase solver failed despite the explicit solution -upperPhase");
    fi;
    correctionDifferential := differential(upperDegree,solutions.particular,true);
    if correctionDifferential <> multiply(solutions.particular)
       or List(correctionDifferential, x -> x mod modulus) <> rhs then
        Error("koAHSS: particular phase correction fails the residual equation");
    fi;
    numerator := upperDifferential + correctionDifferential;
    if not ForAll(numerator, x -> x mod modulus = 0) then
        Error("koAHSS: tertiary numerator is not divisible by the phase modulus");
    fi;
    particularCocycle := List(numerator, x -> QuoInt(x,modulus));
    if not ForAll(differential(outputDegree,particularCocycle,true), x -> x = 0) then
        Error("koAHSS: particular tertiary candidate is not closed");
    fi;
    homogeneousCocycles := [];
    for generator in solutions.homogeneousGenerators do
        numerator := differential(upperDegree,generator,true);
        if numerator <> multiply(generator)
           or not ForAll(numerator, x -> x mod modulus = 0) then
            Error("koAHSS: homogeneous phase correction fails its residual equation");
        fi;
        cocycle := List(numerator, x -> QuoInt(x,modulus));
        if not ForAll(differential(outputDegree,cocycle,true), x -> x = 0) then
            Error("koAHSS: homogeneous tertiary candidate is not closed");
        fi;
        Add(homogeneousCocycles,cocycle);
    od;
    unshiftedParticularCocycle := ShallowCopy(particularCocycle);
    coefficient := 0;
    correction := List([1..outputDimension], j -> 0);
    if Length(arg) = 6 then
        correction := koAHSSTertiaryCorrection(backend,n,arg[6]);
        if correction = fail then return fail; fi;
        particularCocycle := particularCocycle + correction;
        coefficient := 1;
    fi;
    return rec(
        phaseSolutions := solutions,
        upperPhase := upperPhase,
        particularCorrection := ShallowCopy(solutions.particular),
        unshiftedParticularCocycle := unshiftedParticularCocycle,
        particularCocycle := particularCocycle,
        homogeneousCocycles := homogeneousCocycles,
        homogeneousOrders := ShallowCopy(solutions.homogeneousOrders),
        modulus := modulus,
        degree := outputDegree,
        tertiaryCoefficient := coefficient,
        correction := correction,
        correctionApplied := coefficient = 1,
        isLocalCandidateFamily := true,
        isNormalizedOperation := false
    );
end);
