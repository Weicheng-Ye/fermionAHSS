# The explicitly specified k=1 primary correction from the tertiary note.
# This is a correction to an existing tertiary operation, not a reference T.
InstallGlobalFunction(koAHSSTertiaryCorrection, function(backend, n, A)
    local dimension, checkVector, differential, cupIntegral, reduce, s, omega,
          a, square, numerator, bockstein, W, dW, Romega, Qomega, first,
          second, pontryagin, product, correction;
    if not IsInt(n) or n < 0 then
        Error("koAHSS: tertiary correction degree must be a nonnegative integer");
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
        Error("koAHSS: tertiary correction requires a cochain backend with twists");
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
            Error("koAHSS: correction cochain has incorrect dimension or noninteger entries");
        fi;
        if binary and not ForAll(vector, x -> x in [0,1]) then
            Error("koAHSS: correction twists must use binary representatives");
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
    s := backend.twists.s;
    omega := backend.twists.omega;
    checkVector(s, 1, true);
    checkVector(omega, 2, true);
    checkVector(A, n, false);
    if not ForAll(differential(1,s,false), x -> x mod 2 = 0)
       or not ForAll(differential(2,omega,false), x -> x mod 2 = 0) then
        Error("koAHSS: correction twists must be mod-two cocycles");
    fi;
    if not ForAll(differential(n,A,true), x -> x = 0) then
        Error("koAHSS: tertiary correction input must be a sign-integral cocycle");
    fi;
    if ForAll(A, x -> x = 0) then
        return List([1..dimension(n+5)], j -> 0);
    fi;

    a := reduce(A);
    if not ForAll(differential(n,a,false), x -> x mod 2 = 0) then
        Error("koAHSS: sign-integral input does not reduce to a mod-two cocycle");
    fi;
    if n < 4 or ForAll(a, x -> x = 0) then
        square := List([1..dimension(n+4)], j -> 0);
    else
        square := backend.cupMod2(n-4,n,a,n,a);
        checkVector(square,n+4,false);
        square := reduce(square);
    fi;
    numerator := differential(n+4,square,true);
    if not ForAll(numerator, x -> x mod 2 = 0) then
        Error("koAHSS: Sq4 Bockstein numerator is not even");
    fi;
    bockstein := List(numerator, x -> QuoInt(x,2));
    if not ForAll(differential(n+5,bockstein,true), x -> x = 0) then
        Error("koAHSS: Sq4 Bockstein correction is not closed");
    fi;

    W := ShallowCopy(omega);
    dW := differential(2,W,false);
    if not ForAll(dW, x -> x mod 2 = 0) then
        Error("koAHSS: omega Bockstein numerator is not even");
    fi;
    Romega := List(dW, x -> QuoInt(x,2));
    if not ForAll(differential(3,Romega,false), x -> x = 0) then
        Error("koAHSS: omega Bockstein is not closed");
    fi;
    # Qomega vanishes literally when the chosen binary lift is closed.
    # This includes omega=0 and needs no integral cup implementation.
    if ForAll(Romega, x -> x = 0) then return bockstein; fi;
    if (IsBound(backend.hasIntegralCups) and backend.hasIntegralCups = false)
       or not IsBound(backend.cupIntegral)
       or not IsFunction(backend.cupIntegral) then
        return fail;
    fi;
    cupIntegral := function(index, p, left, q, right, signLeft, signRight)
        local value;
        value := backend.cupIntegral(index,p,left,q,right,signLeft,signRight);
        if value = fail then return fail; fi;
        checkVector(value,p+q-index,false);
        return value;
    end;
    first := cupIntegral(0,2,W,3,Romega,false,false);
    if first = fail then return fail; fi;
    second := cupIntegral(1,3,Romega,3,Romega,false,false);
    if second = fail then return fail; fi;
    Qomega := first + second;
    if not ForAll(differential(5,Qomega,false), x -> x = 0) then
        Error("koAHSS: the Pontryagin-square Bockstein cochain is not closed");
    fi;
    first := cupIntegral(0,2,W,2,W,false,false);
    if first = fail then return fail; fi;
    second := cupIntegral(1,2,W,3,dW,false,false);
    if second = fail then return fail; fi;
    pontryagin := first + second;
    if differential(4,pontryagin,false) <> List(Qomega, x -> 4*x) then
        Error("koAHSS: integral cups fail d(PontryaginSquare) = 4 Qomega");
    fi;
    if ForAll(Qomega, x -> x = 0) then return bockstein; fi;
    # The sign-valued input must be the first factor in this representative.
    product := cupIntegral(0,n,A,5,Qomega,true,false);
    if product = fail then return fail; fi;
    correction := bockstein + List(product, x -> (-1)^n*x);
    if not ForAll(differential(n+5,correction,true), x -> x = 0) then
        Error("koAHSS: the k=1 tertiary correction is not closed");
    fi;
    return correction;
end);
