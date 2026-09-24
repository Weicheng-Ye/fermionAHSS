# The primary differences between secondary operations in the supplied notes.
# Coefficients are determined only from supplied comparisons, never defaulted.
InstallGlobalFunction(koAHSSSecondaryNormalizationEquations, function(operation)
    local detection, labels, corrections, outputDegree;
    if operation = "Tau" then
        detection := [[1,0,0],[1,1,0],[0,0,1]];
        labels := ["epsilon3","epsilon4","epsilon5"];
        corrections := ["s^3 rho2", "s omega rho2", "(Sq1 omega) rho2"];
        outputDegree := 12;
    elif operation = "Psi" then
        detection := [[1,1,0],[0,1,1],[0,0,1]];
        labels := ["epsilon2","epsilon3","epsilon4"];
        corrections := ["beta_s(s Sq2)", "beta_s(s^3)", "beta_s(omega Sq1)"];
        outputDegree := 13;
    else
        Error("koAHSS: secondary operation must be Tau or Psi");
    fi;
    return rec(operation := operation, coefficientLabels := labels,
        correctionLabels := corrections, discrepancyLabels := ["d1","d2","d3"],
        comparisonLabels := ["L + R^3", "3L + R", "gamma3 + R"],
        detectionMatrix := detection, matrix := TransposedMat(detection),
        equationConvention := "epsilon * matrix = d modulo two",
        comparisonConvention := "prescribed bundle value minus reference value",
        inputDegree := 9, outputDegree := outputDegree);
end);

InstallGlobalFunction(koAHSSSecondaryCorrections, function(backend,n,input,operation)
    local equations, dimension, differential, cup, closed, s, omega, a,
        ss, sss, somega, sq1omega, sq1a, sq2a, values, degree,
        vector, numerator, result;
    equations := koAHSSSecondaryNormalizationEquations(operation);
    KOAHSS_EquationBackend(backend,n);
    if not IsBound(backend.cupMod2) or not IsFunction(backend.cupMod2)
       or not IsBound(backend.twists) or not IsRecord(backend.twists)
       or not IsBound(backend.twists.s) or not IsBound(backend.twists.omega) then
        Error("koAHSS: secondary corrections require mod-two cups and twists");
    fi;
    dimension := degree -> KOAHSS_EquationDimension(backend,degree);
    differential := function(degree,vector,sign)
        local value;
        KOAHSS_CC_CheckVector(vector,dimension(degree),"correction input");
        value := backend.coboundary(degree,vector,sign);
        KOAHSS_CC_CheckVector(value,dimension(degree+1),"correction coboundary");
        return value;
    end;
    closed := function(degree,vector,sign,modTwo)
        local value;
        value := differential(degree,vector,sign);
        if modTwo then return ForAll(value,x -> x mod 2 = 0); fi;
        return ForAll(value,x -> x = 0);
    end;
    cup := function(p,a,q,b)
        local value;
        value := backend.cupMod2(0,p,a,q,b);
        KOAHSS_CC_CheckVector(value,dimension(p+q),"secondary correction cup");
        return List(value,x -> x mod 2);
    end;
    s := backend.twists.s; omega := backend.twists.omega;
    KOAHSS_CC_CheckVector(s,dimension(1),"s");
    KOAHSS_CC_CheckVector(omega,dimension(2),"omega");
    if not ForAll(s,x -> x in [0,1]) or not ForAll(omega,x -> x in [0,1])
       or not closed(1,s,false,true) or not closed(2,omega,false,true) then
        Error("koAHSS: secondary correction twists must be binary cocycles");
    fi;
    KOAHSS_CC_CheckVector(input,dimension(n),"secondary correction input");
    if operation = "Tau" and not closed(n,input,true,false) then
        Error("koAHSS: Tau correction input must be a sign-integral cocycle");
    fi;
    a := List(input,x -> x mod 2);
    if not closed(n,a,false,true) then
        Error("koAHSS: secondary correction input must reduce to a mod-two cocycle");
    fi;
    ss := cup(1,s,1,s);
    sss := cup(2,ss,1,s);
    if operation = "Tau" then
        somega := cup(1,s,2,omega);
        numerator := differential(2,omega,false);
        sq1omega := List(numerator,x -> QuoInt(x,2) mod 2);
        values := [cup(3,sss,n,a),cup(3,somega,n,a),cup(3,sq1omega,n,a)];
        degree := n+3;
        for vector in values do
            if not closed(degree,vector,false,true) then
                Error("koAHSS: Tau correction is not a mod-two cocycle");
            fi;
        od;
        return values;
    fi;
    sq2a := koAHSSSq(backend,n,a,2);
    numerator := differential(n,a,false);
    sq1a := List(numerator,x -> QuoInt(x,2) mod 2);
    values := [cup(1,s,n+2,sq2a),cup(3,sss,n,a),cup(2,omega,n+1,sq1a)];
    result := [];
    for vector in values do
        numerator := differential(n+3,vector,true);
        if not ForAll(numerator,x -> x mod 2 = 0) then
            Error("koAHSS: Psi correction Bockstein numerator is not even");
        fi;
        vector := List(numerator,x -> QuoInt(x,2));
        if not closed(n+4,vector,true,false) then
            Error("koAHSS: Psi correction is not a sign-integral cocycle");
        fi;
        Add(result,vector);
    od;
    return result;
end);

InstallGlobalFunction(koAHSSSolveSecondaryNormalization, function(operation,comparisons)
    local equations, solutions, mode, records, comparison, group,
        coefficient, value, j, valid, result, rhs;
    equations := koAHSSSecondaryNormalizationEquations(operation);
    if not IsList(comparisons) then
        Error("koAHSS: normalization comparisons must be a list");
    fi;
    solutions := [];
    if Length(comparisons)=3 and ForAll(comparisons,x -> IsInt(x) and x in [0,1]) then
        mode := "supplied bundle discrepancy numbers";
        rhs := ShallowCopy(comparisons);
        for coefficient in Tuples([0,1],3) do
            if List(coefficient*equations.matrix,x -> x mod 2)=rhs then
                Add(solutions,coefficient);
            fi;
        od;
    else
        mode := "supplied quotient-group comparisons";
        records := [];
        for comparison in comparisons do
            if not IsRecord(comparison) or not IsBound(comparison.target)
               or not IsGroup(comparison.target) or not IsAbelian(comparison.target)
               or not IsBound(comparison.reference) or not IsBound(comparison.prescribed)
               or not IsBound(comparison.corrections) or not IsList(comparison.corrections)
               or Length(comparison.corrections)<>3 then
                Error("koAHSS: each comparison needs an abelian target group, reference, prescribed, and three corrections");
            fi;
            group := comparison.target;
            if not comparison.reference in group or not comparison.prescribed in group
               or not ForAll(comparison.corrections,x -> x in group) then
                Error("koAHSS: comparison values must belong to its target quotient group");
            fi;
            if not ForAll(comparison.corrections,x -> x^2=One(group)) then
                Error("koAHSS: normalization corrections must have order dividing two");
            fi;
            Add(records,comparison);
        od;
        for coefficient in Tuples([0,1],3) do
            valid := true;
            for comparison in records do
                value := comparison.reference;
                for j in [1..3] do
                    value := value * comparison.corrections[j]^coefficient[j];
                od;
                if value <> comparison.prescribed then valid := false; break; fi;
            od;
            if valid then Add(solutions,coefficient); fi;
        od;
    fi;
    result := rec(operation := operation, coefficientLabels := equations.coefficientLabels,
        solutions := solutions, comparisonCount := Length(comparisons),
        scope := mode, equations := equations);
    if Length(solutions)=0 then result.status := "inconsistent";
    elif Length(solutions)=1 then
        result.status := "solved";
        result.coefficients := ShallowCopy(solutions[1]);
    else result.status := "ambiguous"; fi;
    return result;
end);
