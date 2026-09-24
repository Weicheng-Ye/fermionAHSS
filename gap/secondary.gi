# Explicit secondary defining equations and the notes' reference cochains.
# Independently solved local helpers are retained as affine families. They
# are not asserted to form a natural, Thom-normalized secondary operation.
BindGlobal("KOAHSS_SecondaryContext", function(backend)
    local context;
    KOAHSS_EquationBackend(backend,0);
    if not IsBound(backend.cupMod2) or not IsFunction(backend.cupMod2) then
        Error("koAHSS: secondary cochains require mod-two higher cups");
    fi;
    context:=rec();
    context.dimension:=n -> KOAHSS_EquationDimension(backend,n);
    context.zero:=n -> List([1..context.dimension(n)],j -> 0);
    context.reduce:=a -> List(a,x -> x mod 2);
    context.check:=function(n,a)
        KOAHSS_CC_CheckVector(a,context.dimension(n),"secondary cochain");
    end;
    context.d:=function(n,a,sign)
        local value;
        context.check(n,a);
        value:=backend.coboundary(n,a,sign);
        context.check(n+1,value);
        return value;
    end;
    context.cup:=function(i,p,a,q,b)
        local value, degree;
        degree:=p+q-i;
        if degree<0 then return []; fi;
        if i<0 or ForAll(a,x -> x=0) or ForAll(b,x -> x=0) then
            return context.zero(degree);
        fi;
        value:=backend.cupMod2(i,p,a,q,b);
        context.check(degree,value);
        return context.reduce(value);
    end;
    context.integralCup:=function(i,p,a,q,b)
        local value, degree;
        degree:=p+q-i;
        if degree<0 then return []; fi;
        if i<0 or ForAll(a,x -> x=0) or ForAll(b,x -> x=0) then
            return context.zero(degree);
        fi;
        if not IsBound(backend.cupIntegral) or not IsFunction(backend.cupIntegral) then
            return fail;
        fi;
        value:=backend.cupIntegral(i,p,a,q,b,false,false);
        if value=fail then return fail; fi;
        context.check(degree,value);
        return value;
    end;
    context.cocycle:=function(n,a)
        context.check(n,a);
        a:=context.reduce(a);
        if ForAny(context.d(n,a,false),x -> x mod 2<>0) then
            Error("koAHSS: secondary input must be a mod-two cocycle");
        fi;
        return a;
    end;
    context.bockstein:=function(n,a)
        local value;
        value:=context.d(n,a,false);
        if ForAny(value,x -> x mod 2<>0) then
            Error("koAHSS: secondary Bockstein numerator is not even");
        fi;
        return List(value,x -> QuoInt(x,2));
    end;
    context.square2:=function(n,a) return context.cup(n-2,n,a,n,a); end;
    context.q2:=function(n,a)
        return context.reduce(context.cup(n-2,n,a,n,a)
            +context.cup(n-1,n,a,n+1,context.reduce(context.d(n,a,false))));
    end;
    return context;
end);

InstallGlobalFunction(koAHSSCartanPrimitive, function(arg)
    local backend, kind, n, x, y, ctx, xy, xx, ex, ey, rhs, args, result;
    if not Length(arg) in [5,6] then
        Error("koAHSSCartanPrimitive(backend, kind, n, x, y[, constraints])");
    fi;
    backend:=arg[1]; kind:=arg[2]; n:=arg[3];
    if not kind in [1,2] or not IsInt(n) or n<0 then
        Error("koAHSS: Cartan kind must be 1 or 2 and input degree nonnegative");
    fi;
    ctx:=KOAHSS_SecondaryContext(backend);
    x:=ctx.cocycle(kind,arg[4]); y:=ctx.cocycle(n,arg[5]);
    xy:=ctx.cup(0,kind,x,n,y);
    xx:=ctx.cup(0,kind,x,kind,x);
    ey:=ctx.reduce(ctx.bockstein(n,y));
    rhs:=ctx.square2(kind+n,xy)+ctx.cup(0,kind,x,n+2,ctx.square2(n,y));
    if kind=1 then
        rhs:=rhs+ctx.cup(0,2,xx,n+1,ey);
    else
        ex:=ctx.reduce(ctx.bockstein(2,x));
        rhs:=rhs+ctx.cup(0,4,xx,n,y)+ctx.cup(0,3,ex,n+1,ey);
    fi;
    rhs:=ctx.reduce(rhs);
    if ForAny(ctx.d(n+kind+2,rhs,false),x -> x mod 2<>0) then
        return rec(status:="unsatisfied-identity",kind:=kind,inputDegree:=n,
            primitiveDegree:=n+kind+1,rhs:=rhs,
            reason:="the Cartan relation cochain is not closed in the supplied conventions");
    fi;
    args:=[backend,n+kind+1,rhs];
    if Length(arg)=6 then Add(args,arg[6]); fi;
    result:=CallFuncList(koAHSSSolveCochainEquation,args);
    if result=fail then
        return rec(status:="obstructed",kind:=kind,inputDegree:=n,
            primitiveDegree:=n+kind+1,rhs:=rhs,
            reason:="the Cartan equation has no solution with the supplied constraints");
    fi;
    result.kind:=kind; result.inputDegree:=n;
    result.equation:=rec(name:=Concatenation("Cartan zeta",String(kind)),
        coefficientModulus:=2,primitiveDegree:=n+kind+1,rhs:=rhs);
    result.scope:="this cochain complex and these input cocycles";
    result.isNormalizedOperation:=false;
    return result;
end);

InstallGlobalFunction(koAHSSSecondaryReference, function(arg)
    local backend, n, input, options, inputType, ctx, s, omega, a, A,
          B, e, CB, u, v, w, primary, b, chi, zeta2, zeta1e, zeta1a,
          K, t, x, y, F, G, q, L, psi, tau, M, dq, first, second,
          omegaS, ss, numerator, helpers, constraints, result, family,
          selectPrimitive, solveHelper, recordIdentity, stop, choice, args;
    if not Length(arg) in [3,4] then
        Error("koAHSSSecondaryReference(backend, n, input[, options])");
    fi;
    backend:=arg[1]; n:=arg[2]; input:=arg[3];
    if not IsInt(n) or n<0 then Error("koAHSS: secondary degree must be nonnegative"); fi;
    options:=rec(); if Length(arg)=4 then options:=arg[4]; fi;
    if not IsRecord(options) then Error("koAHSS: secondary options must be a record"); fi;
    inputType:="integral";
    if IsBound(options.inputType) then inputType:=options.inputType; fi;
    if not inputType in ["integral","mod2"] then
        Error("koAHSS: inputType must be integral or mod2");
    fi;
    helpers:=rec(); constraints:=rec();
    if IsBound(options.helpers) then helpers:=options.helpers; fi;
    if IsBound(options.constraints) then constraints:=options.constraints; fi;
    if not IsRecord(helpers) or not IsRecord(constraints) then
        Error("koAHSS: helpers and constraints must be records");
    fi;
    ctx:=KOAHSS_SecondaryContext(backend);
    if not IsBound(backend.twists) or not IsRecord(backend.twists)
       or not IsBound(backend.twists.s) or not IsBound(backend.twists.omega) then
        Error("koAHSS: secondary cochains require backend.twists.s and omega");
    fi;
    s:=ctx.cocycle(1,backend.twists.s);
    omega:=ctx.cocycle(2,backend.twists.omega);
    result:=rec(status:="candidate",inputDegree:=n,inputType:=inputType,
        cochains:=rec(s:=s,omega:=omega),helpers:=rec(),equations:=[],
        isLocalCandidate:=true,isNormalizedOperation:=false,isStableFamily:=false,
        scope:="one defining system on this cochain complex",
        noteStableDegree:=n>=6);
    stop:=function(status,stage,reason)
        result.status:=status; result.failedStage:=stage; result.reason:=reason;
        return result;
    end;
    recordIdentity:=function(name,degree,lhs,rhs,modulus)
        local identity;
        if modulus=2 then lhs:=ctx.reduce(lhs); rhs:=ctx.reduce(rhs); fi;
        identity:=rec(name:=name,degree:=degree,lhs:=lhs,rhs:=rhs,
            coefficientModulus:=modulus,satisfied:=lhs=rhs);
        Add(result.equations,identity);
        return identity.satisfied;
    end;
    solveHelper:=function(name,solver,arguments)
        local family;
        if IsBound(constraints.(name)) then Add(arguments,constraints.(name)); fi;
        family:=CallFuncList(solver,arguments);
        if family=fail then
            family:=rec(status:="obstructed",reason:="the defining equation has no solution");
        fi;
        result.helpers.(name):=family;
        return family;
    end;
    selectPrimitive:=function(name,family)
        local selected, supplied;
        if not IsBound(family.status) or family.status<>"solved" then return fail; fi;
        supplied:=fail;
        if name="b" and IsBound(options.b) then supplied:=options.b;
        elif IsBound(helpers.(name)) then supplied:=helpers.(name); fi;
        selected:=family.primitive;
        family.choiceSource:="local equation solution";
        if supplied<>fail then
            if IsRecord(supplied) and IsBound(supplied.selectedPrimitive) then
                supplied:=supplied.selectedPrimitive;
            elif IsRecord(supplied) and IsBound(supplied.primitive) then
                supplied:=supplied.primitive;
            fi;
            ctx.check(family.primitiveDegree,supplied);
            selected:=ctx.reduce(supplied);
            family.choiceSource:="supplied cochain";
        fi;
        if not recordIdentity(Concatenation("d ",name),family.primitiveDegree+1,
            ctx.d(family.primitiveDegree,selected,false),family.rhs,2) then
            family.status:="unsatisfied-identity";
            family.reason:="supplied primitive does not satisfy its defining equation";
            return fail;
        fi;
        if IsBound(constraints.(name)) then
            supplied:=constraints.(name);
            if Length(selected)=0 then choice:=List(supplied.rhs,j -> 0);
            elif Length(supplied.rhs)=0 then choice:=[];
            else choice:=selected*supplied.matrix; fi;
            if ctx.reduce(choice)<>ctx.reduce(supplied.rhs) then
                family.status:="unsatisfied-identity";
                family.reason:="supplied primitive violates its linear constraints";
                return fail;
            fi;
        fi;
        family.selectedPrimitive:=ShallowCopy(selected);
        return selected;
    end;

    ctx.check(n,input);
    if inputType="integral" then
        A:=ShallowCopy(input);
        if ForAny(ctx.d(n,A,true),x -> x<>0) then
            Error("koAHSS: integral secondary input must be a sign-integral cocycle");
        fi;
        a:=ctx.cocycle(n,ctx.reduce(A)); result.cochains.A:=A;
    else a:=ctx.cocycle(n,input); fi;
    B:=ctx.bockstein(n,a); e:=ctx.reduce(B);
    CB:=List(B+e,x -> QuoInt(x,2));
    u:=ctx.square2(n,a); v:=ctx.cup(0,2,omega,n,a);
    w:=ctx.cup(0,1,s,n+1,e); primary:=ctx.reduce(u+v+w);
    result.cochains.a:=a; result.cochains.B:=B; result.cochains.e:=e;
    result.cochains.CB:=CB; result.cochains.primary:=primary;
    family:=solveHelper("b",koAHSSSolveCochainEquation,[backend,n+1,primary]);
    b:=selectPrimitive("b",family);
    if b=fail then return stop(family.status,"b",family.reason); fi;
    result.cochains.b:=b;
    family:=solveHelper("chi",koAHSSAdem22Primitive,[backend,n,a]);
    chi:=selectPrimitive("chi",family);
    if chi=fail then return stop(family.status,"chi",family.reason); fi;
    family:=solveHelper("zeta2",koAHSSCartanPrimitive,[backend,2,n,omega,a]);
    zeta2:=selectPrimitive("zeta2",family);
    if zeta2=fail then return stop(family.status,"zeta2",family.reason); fi;
    family:=solveHelper("zeta1e",koAHSSCartanPrimitive,[backend,1,n+1,s,e]);
    zeta1e:=selectPrimitive("zeta1e",family);
    if zeta1e=fail then return stop(family.status,"zeta1e",family.reason); fi;

    omegaS:=ctx.cup(1,2,omega,1,s); ss:=ctx.cup(0,1,s,1,s);
    F:=ctx.reduce(ctx.q2(n+1,b)+ctx.cup(0,2,omega,n+1,b)+zeta2+chi
        +ctx.cup(n+1,n+2,u,n+2,v)+ctx.cup(n+1,n+2,u,n+2,w)
        +ctx.cup(n+1,n+2,v,n+2,w)+zeta1e
        +ctx.cup(0,2,omegaS,n+1,e)+ctx.cup(0,1,s,n+2,u)
        +ctx.cup(0,2,ss,n+1,ctx.reduce(CB)));
    result.cochains.F:=F;
    first:=ctx.integralCup(0,2,omega,n+1,B);
    second:=ctx.integralCup(n-1,n+1,B,n+1,B);
    if first=fail or second=fail then
        return stop("unavailable","q","signed integral higher cups are unavailable");
    fi;
    q:=first+second; result.cochains.q:=q;
    dq:=ctx.d(n+3,q,true);
    if not recordIdentity("d_s q is even",n+4,dq,ctx.zero(n+4),2) then
        return stop("unsatisfied-identity","q parity","the reference q has odd signed coboundary");
    fi;
    dq:=List(dq,x -> QuoInt(x,2));
    if not recordIdentity("d F = rho(d_s q / 2)",n+4,ctx.d(n+3,F,false),dq,2) then
        return stop("unsatisfied-identity","F compatibility",
            "the notes' exact F identity fails in the supplied cochain conventions");
    fi;
    numerator:=ctx.d(n+3,F,true)+dq;
    if not recordIdentity("tau numerator is even",n+4,numerator,ctx.zero(n+4),2) then
        return stop("unsatisfied-identity","tau parity","the tau numerator is not even");
    fi;
    tau:=List(numerator,x -> QuoInt(x,2)); result.cochains.tau:=tau;
    if not recordIdentity("d_s tau = 0",n+5,ctx.d(n+4,tau,true),ctx.zero(n+5),0) then
        return stop("unsatisfied-identity","tau closure","the integral secondary cochain is not closed");
    fi;
    if inputType="mod2" then
        result.operation:="Psi"; result.noteOperation:="Theta";
        result.operationDegree:=n+4; result.cochain:=tau;
        return result;
    fi;

    K:=List(A-a,x -> QuoInt(x,2)); t:=ctx.reduce(K);
    x:=ctx.reduce(ctx.d(n,t,false)); y:=ctx.cup(0,1,s,n,a);
    result.cochains.K:=K; result.cochains.t:=t;
    if not recordIdentity("d t = e + s a",n+1,x,e+y,2) then
        return stop("unsatisfied-identity","integral carry",
            "the notes' exact twisted carry identity fails in the supplied cochain conventions");
    fi;
    family:=solveHelper("zeta1a",koAHSSCartanPrimitive,[backend,1,n,s,a]);
    zeta1a:=selectPrimitive("zeta1a",family);
    if zeta1a=fail then return stop(family.status,"zeta1a",family.reason); fi;
    G:=ctx.reduce(ctx.q2(n,t)+ctx.cup(0,2,omega,n,t)+ctx.cup(n,n+1,x,n+1,y)
        +zeta1a+ctx.cup(0,2,omegaS,n,a)+ctx.cup(0,1,s,n+1,b));
    result.cochains.G:=G;
    if not recordIdentity("d G = rho q",n+3,ctx.d(n+2,G,false),q,2) then
        return stop("unsatisfied-identity","G compatibility",
            "the notes' exact G identity fails in the supplied cochain conventions");
    fi;
    numerator:=q-ctx.d(n+2,G,true);
    if not recordIdentity("L numerator is even",n+3,numerator,ctx.zero(n+3),2) then
        return stop("unsatisfied-identity","L parity","the L numerator is not even");
    fi;
    L:=List(numerator,x -> QuoInt(x,2));
    psi:=ctx.reduce(F+L); M:=F+L;
    result.cochains.L:=L; result.cochains.psi:=psi; result.cochains.M:=M;
    if not recordIdentity("d psi = 0",n+4,ctx.d(n+3,psi,false),ctx.zero(n+4),2) then
        return stop("unsatisfied-identity","psi closure","the mod-two secondary cochain is not closed");
    fi;
    if not recordIdentity("d_s M = 2 tau",n+4,ctx.d(n+3,M,true),2*tau,0) then
        return stop("unsatisfied-identity","lower compatibility","the exact common lower boundary fails");
    fi;
    if not recordIdentity("rho M = psi",n+3,M,psi,2) then
        return stop("unsatisfied-identity","lower reduction","the exact common lower reduction fails");
    fi;
    result.operation:="Tau"; result.noteOperation:="Psi";
    result.operationDegree:=n+3; result.cochain:=psi;
    return result;
end);
