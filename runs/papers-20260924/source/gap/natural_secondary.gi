# Natural simplicial operations on the group bar model, projected to the
# supplied resolution. Only the input nullhomotopy b is solved locally.
BindGlobal("KOAHSS_NaturalCochains",function(backend)
    local transport, ctx, make, s, omega, cacheLimit, checkSign, directLift;
    if IsBound(backend.naturalTransport) then transport:=backend.naturalTransport();
    elif IsBound(backend.naturalBar) then transport:=backend.naturalBar();
    else return fail; fi;
    ctx:=rec(transport:=transport);
    cacheLimit:=256;
    if IsBound(transport.cacheEntries) then cacheLimit:=transport.cacheEntries; fi;
    make:=function(n,evaluate)
        local cache,value,order,next,memoToken;
        if IsBound(transport.cochainKey) then cache:=NewDictionary("",true);
        else cache:=NewDictionary([],true); fi;
        order:=[]; next:=1;
        value:=function(simplex)
            local answer,key;
            if Length(simplex)<>n+1 then Error("natural cochain input degree mismatch"); fi;
            if IsBound(transport.isDegenerate) then
                if transport.isDegenerate(simplex) then return 0; fi;
            elif ForAny([1..n],j->simplex[j]=simplex[j+1]) then return 0; fi;
            # LOCAL values are invariant under simultaneous group translation.
            # Bound memo tables: nonlinear expressions otherwise retain all faces.
            # General simplicial test models need not supply this method.
            if IsBound(transport.normalizeSimplex) then simplex:=transport.normalizeSimplex(simplex); fi;
            key:=simplex;
            if IsBound(transport.cochainKey) then key:=transport.cochainKey(simplex); fi;
            answer:=LookupDictionary(cache,key);
            if answer=fail then
                answer:=evaluate(simplex);
                if not IsInt(answer) then Error("natural cochain must be integral"); fi;
                if cacheLimit>0 then
                    # Retain the other entries when one new value arrives.
                    # A fixed FIFO ring keeps the same strict entry bound.
                    key:=Immutable(key);
                    if Length(order)<cacheLimit then Add(order,key);
                    else
                        RemoveDictionary(cache,order[next]); order[next]:=key;
                        next:=next mod cacheLimit+1;
                    fi;
                    AddDictionary(cache,key,answer);
                fi;
            fi;
            return answer;
        end;
        if IsBound(transport.scalarMemo) then
            # Reuse this exact scalar expression on an immutable input face
            # before repeating its normalization and full memo-key work.
            # The transport owns the bounded projection scope; a distinct
            # mutable token identifies this closure without cochain equality.
            memoToken:=rec();
            return rec(degree:=n,value:=sigma->transport.scalarMemo(sigma,memoToken,value));
        fi;
        return rec(degree:=n,value:=value);
    end;
    ctx.make:=make;
    # For affine DK simplices and a cocycle sign, top-face evaluation with
    # the first-vertex frame is already invariant under group translation.
    # Keep legacy evaluation for structurally accepted non-affine simplices.
    directLift:=function(n,evaluate,sign,evaluateMod2)
        local fallback,checked;
        fallback:=fail;
        checked:=function(simplex,binary)
            local answer;
            if Length(simplex)<>n+1 then Error("natural cochain input degree mismatch"); fi;
            if transport.isDegenerate(simplex) then return 0; fi;
            if not transport.localCoordinateLiftAllowed(simplex,sign) then
                if fallback=fail then fallback:=make(n,evaluate); fi;
                answer:=fallback.value(simplex);
            elif binary then answer:=evaluateMod2(simplex);
            else answer:=evaluate(simplex); fi;
            if not IsInt(answer) then Error("natural cochain must be integral"); fi;
            if binary then return answer mod 2; fi;
            return answer;
        end;
        # All coefficient characters and the first-vertex frame are +/-1,
        # so reduction modulo two uses the unsigned top-face evaluation.
        # The eligibility and legacy fallback remain exactly the signed ones.
        return rec(degree:=n,directDKCoordinateLift:=true,
            value:=sigma->checked(sigma,false),mod2Value:=sigma->checked(sigma,true));
    end;
    checkSign:=function(sign)
        if IsInt(sign) and sign=0 then return; fi;
        KOAHSS_CC_CheckVector(sign,backend.dimension(1),"natural cochain sign");
        if not ForAll(sign,x->x in [0,1]) then
            Error("natural cochain sign must be zero or a binary degree-one vector");
        fi;
    end;
    # These flags certify literal zero on the entire envelope, never merely
    # zero after projection. In particular an even integral lift reduces to
    # zero before any large homotopy chain needs to be constructed.
    ctx.zero:=function(n)
        local a;
        a:=make(n,sigma->0); a.knownZero:=true; a.mod2Zero:=true;
        return a;
    end;
    ctx.isZero:=a->IsBound(a.knownZero) and a.knownZero;
    ctx.lift:=function(n,a,sign)
        local value;
        KOAHSS_CC_CheckVector(a,backend.dimension(n),"natural cochain lift");
        checkSign(sign);
        a:=Immutable(a);
        if IsList(sign) then sign:=Immutable(sign); fi;
        if ForAll(a,x->x=0) then return ctx.zero(n); fi;
        if IsBound(transport.localCoordinateLiftAllowed) and
           ((IsInt(sign) and sign=0) or
            ForAll(backend.coboundary(1,sign,false),x->x mod 2=0)) then
            value:=directLift(n,transport.lift(n,a,sign),sign,transport.lift(n,a,0));
        else
            # Arbitrary binary noncocycles were accepted here historically.
            # Their raw frame need not be G-invariant; preserve make's path.
            value:=make(n,transport.lift(n,a,sign));
        fi;
        value.mod2Zero:=ForAll(a,x->x mod 2=0);
        return value;
    end;
    ctx.project:=function(a,sign)
        checkSign(sign);
        if ctx.isZero(a) then return List([1..backend.dimension(a.degree)],j->0); fi;
        return transport.project(a.degree,a.value,sign);
    end;
    ctx.projectMany:=function(cochains,signs)
        local degree,active,answer,values,j;
        if not IsList(cochains) or IsEmpty(cochains) or not IsList(signs)
           or Length(cochains)<>Length(signs)
           or not ForAll(cochains,a->a.degree=cochains[1].degree) then
            Error("natural simultaneous projection needs cochains in one degree and one sign per cochain");
        fi;
        degree:=cochains[1].degree;
        for j in [1..Length(signs)] do checkSign(signs[j]); od;
        if not IsBound(transport.projectMany) then
            return List([1..Length(cochains)],j->ctx.project(cochains[j],signs[j]));
        fi;
        answer:=List(cochains,a->List([1..backend.dimension(degree)],j->0));
        active:=Filtered([1..Length(cochains)],j->not ctx.isZero(cochains[j]));
        if not IsEmpty(active) then
            values:=transport.projectMany(degree,
                List(active,j->cochains[j].value),signs{active});
            for j in [1..Length(active)] do answer[active[j]]:=values[j]; od;
        fi;
        return answer;
    end;
    ctx.reduce:=function(a)
        if ctx.isZero(a) or (IsBound(a.mod2Zero) and a.mod2Zero) then
            return ctx.zero(a.degree);
        fi;
        if IsBound(a.directDKCoordinateLift) and a.directDKCoordinateLift=true then
            # The underlying value retains all validation; applying mod two
            # adds no geometry query and requires no second simplex memo.
            if IsBound(a.mod2Value) then
                return rec(degree:=a.degree,directDKCoordinateLift:=true,
                    value:=a.mod2Value,mod2Value:=a.mod2Value);
            fi;
            return rec(degree:=a.degree,directDKCoordinateLift:=true,
                value:=sigma->a.value(sigma) mod 2);
        fi;
        return make(a.degree,sigma->a.value(sigma) mod 2);
    end;
    ctx.scale:=function(k,a)
        if k=0 or ctx.isZero(a) then return ctx.zero(a.degree); fi;
        return make(a.degree,sigma->k*a.value(sigma));
    end;
    ctx.add:=function(terms)
        local degree;
        if IsEmpty(terms) or not ForAll(terms,a->a.degree=terms[1].degree) then
            Error("natural cochain sum needs a nonempty list in one degree");
        fi;
        degree:=terms[1].degree;
        terms:=Filtered(terms,a->not ctx.isZero(a));
        if IsEmpty(terms) then return ctx.zero(degree); fi;
        if Length(terms)=1 then return terms[1]; fi;
        return make(terms[1].degree,sigma->Sum(terms,a->a.value(sigma)));
    end;
    ctx.divide:=function(a,denominator)
        if not IsInt(denominator) or denominator=0 then
            Error("natural cochain denominator must be a nonzero integer");
        fi;
        if ctx.isZero(a) then return ctx.zero(a.degree); fi;
        return make(a.degree,function(sigma)
            local value;
            value:=a.value(sigma);
            if value mod denominator<>0 then Error("natural cochain numerator is not divisible"); fi;
            return QuoInt(value,denominator);
        end);
    end;
    s:=ctx.reduce(ctx.lift(1,backend.twists.s,0));
    omega:=ctx.reduce(ctx.lift(2,backend.twists.omega,0));
    ctx.s:=s; ctx.omega:=omega;
    ctx.d:=function(a,sign)
        if ctx.isZero(a) then return ctx.zero(a.degree+1); fi;
        return make(a.degree+1,function(sigma)
            local value,j,face;
            value:=0;
            for j in [1..Length(sigma)] do
                face:=sigma{Filtered([1..Length(sigma)],k->k<>j)};
                if j=1 and sign then
                    value:=value+(-1)^s.value(sigma{[1,2]})*a.value(face);
                else value:=value+(-1)^(j-1)*a.value(face); fi;
            od;
            return value;
        end);
    end;
    ctx.cup:=function(i,a,b,integral)
        local word,n;
        n:=a.degree+b.degree-i;
        if n<0 then Error("negative natural cochain degree"); fi;
        if i<0 or ctx.isZero(a) or ctx.isZero(b) then return ctx.zero(n); fi;
        word:=List([0..i+1],j->1+j mod 2);
        if integral then
            return make(n,sigma->koAHSSNaturalWordValue(word,[a.degree,b.degree],
                [a.value,b.value],sigma,i));
        fi;
        return make(n,sigma->koAHSSNaturalWordValue(word,[a.degree,b.degree],
            [a.value,b.value],sigma));
    end;
    ctx.product:=function(a,b) return ctx.cup(0,a,b,false); end;
    ctx.square2:=a->ctx.cup(a.degree-2,a,a,false);
    ctx.q2:=a->ctx.reduce(ctx.add([ctx.square2(a),
        ctx.cup(a.degree-1,a,ctx.reduce(ctx.d(a,false)),false)]));
    ctx.zeta:=function(kind,x,a)
        if not kind in [1,2] then Error("natural Cartan helper kind must be one or two"); fi;
        if a.degree=0 or ctx.isZero(x) or ctx.isZero(a) then
            return ctx.zero(a.degree+kind+1);
        fi;
        return make(a.degree+kind+1,sigma->koAHSSNaturalZetaValue(
            kind,a.degree,x.value,a.value,sigma));
    end;
    ctx.chi:=function(a)
        local faces,cache,order,next;
        if ctx.isZero(a) then return ctx.zero(a.degree+3); fi;
        faces:=Combinations([1..a.degree+4],a.degree+1);
        cache:=NewDictionary([],true); order:=[]; next:=1;
        return make(a.degree+3,function(sigma)
            local values,value;
            values:=List(faces,face->a.value(sigma{face}) mod 2);
            if ForAll(values,x->x=0) then return 0; fi;
            value:=LookupDictionary(cache,values);
            if value=fail then
                value:=koAHSSNaturalChiValue(a.degree,a.value,sigma,"tail");
                if cacheLimit>0 then
                    MakeImmutable(values);
                    if Length(order)<cacheLimit then Add(order,values);
                    else
                        RemoveDictionary(cache,order[next]); order[next]:=values;
                        next:=next mod cacheLimit+1;
                    fi;
                    AddDictionary(cache,values,value);
                fi;
            fi;
            return value;
        end);
    end;
    return ctx;
end);

BindGlobal("KOAHSS_NaturalPrimaryData",function(backend,n,input,inputType)
    local ctx,a,A,B,e,u,v,w,c,legacy,fast,scalarE,word,transport;
    ctx:=KOAHSS_NaturalCochains(backend);
    if ctx=fail then return fail; fi;
    A:=ctx.lift(n,input,0);
    if inputType="integral" then A:=ctx.lift(n,input,backend.twists.s); fi;
    a:=ctx.reduce(A); B:=ctx.divide(ctx.d(a,false),2); e:=ctx.reduce(B);
    u:=ctx.square2(a); v:=ctx.product(ctx.omega,a); w:=ctx.product(ctx.s,e);
    c:=ctx.reduce(ctx.add([u,v,w]));
    transport:=ctx.transport;
    if n<=7 and not ctx.isZero(c) and
       IsBound(transport.localScalarExpressionAllowed) then
        legacy:=c;
        scalarE:=function(sigma)
            local value,j,face;
            if Length(sigma)<>n+2 then Error("natural cochain input degree mismatch"); fi;
            if transport.isDegenerate(sigma) then return 0; fi;
            value:=0;
            for j in [1..Length(sigma)] do
                face:=sigma{Filtered([1..Length(sigma)],k->k<>j)};
                value:=value+(-1)^(j-1)*a.value(face);
            od;
            if not IsInt(value) then Error("natural cochain must be integral"); fi;
            if value mod 2<>0 then Error("natural cochain numerator is not divisible"); fi;
            return QuoInt(value,2) mod 2;
        end;
        word:=List([0..n-1],j->1+j mod 2);
        fast:=ctx.make(n+2,function(sigma)
            local value;
            value:=0;
            # Keep the existing word evaluator's cut order and zero-product
            # short circuit; in particular s(face)=0 must not force delta a/2.
            if not ctx.isZero(u) then
                value:=value+koAHSSNaturalWordValue(word,[n,n],[a.value,a.value],sigma);
            fi;
            if not ctx.isZero(v) then
                value:=value+koAHSSNaturalWordValue([1,2],[2,n],
                    [ctx.omega.value,a.value],sigma);
            fi;
            if not ctx.isZero(w) then
                value:=value+koAHSSNaturalWordValue([1,2],[1,n+1],
                    [ctx.s.value,scalarE],sigma);
            fi;
            return value mod 2;
        end);
        # Decide on the ORIGINAL simplex, before entering the outer memo.
        # Non-affine or noncanonical inputs must not be normalized twice.
        c:=rec(degree:=n+2,primaryScalarSpecialized:=true,value:=function(sigma)
            if Length(sigma)<>n+3 or transport.isDegenerate(sigma) or
               not transport.localScalarExpressionAllowed(sigma) then
                return legacy.value(sigma);
            fi;
            return fast.value(sigma);
        end);
    fi;
    return rec(ctx:=ctx,A:=A,a:=a,B:=B,e:=e,u:=u,v:=v,w:=w,primary:=c);
end);

BindGlobal("KOAHSS_NaturalPrimary",function(backend,name,n,input)
    local data,ctx,value;
    KOAHSS_CC_CheckVector(input,backend.dimension(n),"natural primary input");
    # These primary operations send the zero cochain to zero identically.
    # Avoid constructing high-dimensional simplicial comparison chains for it.
    if ForAll(input,x->x mod 2=0) then
        if name in ["D","Dbar"] then
            return List([1..backend.dimension(n+2)],j->0);
        elif name="Dtilde" then
            return List([1..backend.dimension(n+3)],j->0);
        fi;
    fi;
    data:=KOAHSS_NaturalPrimaryData(backend,n,input,"mod2");
    if data=fail then Error("natural primary requires a simplicial transport"); fi;
    ctx:=data.ctx;
    if name in ["D","Dbar"] then
        return List(ctx.project(data.primary,0),x->x mod 2);
    elif name="Dtilde" then
        value:=ctx.reduce(ctx.add([data.u,data.v]));
        # Integral chain transfer commutes with d_s. Differentiate after
        # transfer to avoid evaluating a redundant higher simplicial dimension.
        value:=backend.coboundary(n+2,ctx.project(value,backend.twists.s),true);
        if ForAny(value,x->x mod 2<>0) then Error("natural Bockstein boundary is not even"); fi;
        return List(value,x->QuoInt(x,2));
    fi;
    Error("unknown natural primary operation");
end);

InstallGlobalFunction(koAHSSNaturalSecondary,function(arg)
    local backend,n,input,options,inputType,data,ctx,s,omega,A,a,B,e,CB,u,v,w,
          primary,primaryR,family,bR,b,ss,omegaS,F,q,tau,K,t,x,y,G,L,psi,M,
          z,h,kappa,result,cochains,check,project,zero,correctedTau,correctedM,projected;
    if not Length(arg) in [3,4] then
        Error("koAHSSNaturalSecondary(backend,n,input[,options])");
    fi;
    backend:=arg[1]; n:=arg[2]; input:=arg[3]; options:=rec();
    if Length(arg)=4 then options:=arg[4]; fi;
    if not IsInt(n) or n<0 or not IsRecord(options) then Error("invalid natural secondary arguments"); fi;
    inputType:="integral";
    if IsBound(options.inputType) then inputType:=options.inputType; fi;
    if not inputType in ["integral","mod2"] then Error("invalid natural secondary input type"); fi;
    if n>7 then return rec(status:="unavailable",failedStage:="natural helper degree",
        reason:="the implemented calibrated helper family ends at input degree seven"); fi;
    KOAHSS_CC_CheckVector(input,backend.dimension(n),"natural secondary input");
    if inputType="integral" and ForAny(backend.coboundary(n,input,true),x->x<>0) then
        Error("natural integral secondary input is not a sign cocycle");
    fi;
    if ForAny(backend.coboundary(n,input,false),x->x mod 2<>0) then
        Error("natural secondary input is not a mod-two cocycle");
    fi;
    data:=KOAHSS_NaturalPrimaryData(backend,n,input,inputType);
    if data=fail then return rec(status:="unavailable",failedStage:="cochain model",
        reason:="natural secondary evaluation requires a simplicial transport"); fi;
    ctx:=data.ctx; s:=ctx.s; omega:=ctx.omega; A:=data.A; a:=data.a;
    B:=data.B; e:=data.e; u:=data.u; v:=data.v; w:=data.w; primary:=data.primary;
    primaryR:=List(ctx.project(primary,0),x->x mod 2);
    if IsBound(options.b) then
        KOAHSS_CC_CheckVector(options.b,backend.dimension(n+1),"natural nullhomotopy b");
        bR:=List(options.b,x->x mod 2);
        if List(backend.coboundary(n+1,bR,false),x->x mod 2)<>primaryR then
            Error("natural secondary b does not bound the transported primary cochain");
        fi;
    else
        family:=koAHSSSolveCochainEquation(backend,n+1,primaryR);
        if family=fail then return rec(status:="obstructed",failedStage:="b",
            reason:="the primary obstruction is nonzero"); fi;
        bR:=family.primitive;
    fi;
    # dH+Hd=1-gf supplies the comparison between the simplicial model and small-resolution
    # nullhomotopies. No space-dependent Adem or Cartan primitive is chosen.
    if ctx.isZero(primary) then
        b:=ctx.reduce(ctx.lift(n+1,bR,0));
    elif IsBound(ctx.transport.closedMod2HomotopyValue) then
        # primary is closed modulo two; reduction below retains exactly the
        # same binary defining cochain as the full normalized homotopy K.
        b:=ctx.reduce(ctx.add([ctx.lift(n+1,bR,0),ctx.make(n+1,
            sigma->ctx.transport.closedMod2HomotopyValue(n+2,primary.value,sigma))]));
    else
        b:=ctx.reduce(ctx.add([ctx.lift(n+1,bR,0),ctx.make(n+1,
            sigma->ctx.transport.homotopyValue(n+2,primary.value,sigma))]));
    fi;
    CB:=ctx.divide(ctx.add([B,e]),2);
    ss:=ctx.product(s,s); omegaS:=ctx.cup(1,omega,s,false);
    F:=ctx.reduce(ctx.add([ctx.q2(b),ctx.product(omega,b),ctx.zeta(2,omega,a),ctx.chi(a),
        ctx.cup(n+1,u,v,false),ctx.cup(n+1,u,w,false),ctx.cup(n+1,v,w,false),
        ctx.zeta(1,s,e),ctx.product(omegaS,e),ctx.product(s,u),ctx.product(ss,ctx.reduce(CB))]));
    q:=ctx.add([ctx.cup(0,omega,B,true),ctx.cup(n-1,B,B,true)]);
    tau:=ctx.divide(ctx.add([ctx.d(F,true),ctx.divide(ctx.d(q,true),2)]),2);
    # Calibrated eta=(1,0,1). Keep the integral sum of binary lifts so its
    # coboundary is exactly the sum of the two Bockstein representatives.
    z:=ctx.add([ctx.product(s,u),ctx.product(omega,e)]);
    correctedTau:=ctx.add([tau,ctx.divide(ctx.d(z,true),2)]);
    result:=rec(status:="candidate",inputDegree:=n,inputType:=inputType,
        referenceConvention:="chi-tail-degree7-epsilon100-eta101",
        isLocalCandidate:=false,isNormalizedOperation:=true,isStableFamily:=true,
        stableDegreeRange:=[0,7],normalizationCoefficients:=rec(Tau:=[1,0,0],Psi:=[1,0,1]),
        helpers:=rec(convention:="explicit tail-suspension words"),equations:=[],
        scope:="calibrated simplicial formulas transported to this resolution");
    cochains:=rec(b:=ShallowCopy(bR),primary:=primaryR);
    result.cochains:=cochains;
    project:=function(c,sign) return ctx.project(c,sign); end;
    zero:=degree->List([1..backend.dimension(degree)],j->0);
    check:=function(name,degree,lhs,rhs,modulus)
        local satisfied;
        if modulus=2 then lhs:=List(lhs,x->x mod 2); rhs:=List(rhs,x->x mod 2); fi;
        satisfied:=lhs=rhs;
        Add(result.equations,rec(name:=name,degree:=degree,coefficientModulus:=modulus,
            lhs:=lhs,rhs:=rhs,satisfied:=satisfied));
        if not satisfied then Error("natural secondary identity failed: ",name); fi;
    end;
    if inputType="mod2" then
        # tau' = d_s(2F+q+2z)/4. Project the numerator first, using the
        # integral chain map, and keep the divisibility check exact.
        cochains.tau:=backend.coboundary(n+3,
            project(ctx.add([ctx.scale(2,F),q,ctx.scale(2,z)]),backend.twists.s),true);
        if ForAny(cochains.tau,x->x mod 4<>0) then
            Error("natural secondary numerator is not divisible by four");
        fi;
        cochains.tau:=List(cochains.tau,x->QuoInt(x,4));
        check("d_s tau = 0",n+5,backend.coboundary(n+4,cochains.tau,true),zero(n+5),0);
        result.operation:="Psi"; result.noteOperation:="Theta";
        result.operationDegree:=n+4; result.cochain:=cochains.tau;
        return result;
    fi;
    K:=ctx.divide(ctx.add([A,ctx.scale(-1,a)]),2); t:=ctx.reduce(K);
    x:=ctx.reduce(ctx.d(t,false)); y:=ctx.product(s,a);
    G:=ctx.reduce(ctx.add([ctx.q2(t),ctx.product(omega,t),ctx.cup(n,x,y,false),
        ctx.zeta(1,s,a),ctx.product(omegaS,a),ctx.product(s,b)]));
    L:=ctx.divide(ctx.add([q,ctx.scale(-1,ctx.d(G,true))]),2);
    M:=ctx.add([F,L]);
    # For epsilon=(1,0,0), kappa=s^3 a. The exact mod-two
    # identity z+kappa=d(s b+s^2 t+omega t+(omega cup_1 s)a)
    # supplies a common integral
    # lift; no additional integer primitive has to be selected.
    h:=ctx.reduce(ctx.add([ctx.product(s,b),ctx.product(ss,t),ctx.product(omega,t),
        ctx.product(omegaS,a)]));
    kappa:=ctx.product(ctx.product(ss,s),a);
    psi:=ctx.reduce(ctx.add([M,kappa]));
    correctedM:=ctx.add([M,z,ctx.scale(-1,ctx.d(h,true))]);
    cochains.A:=ShallowCopy(input); cochains.a:=List(input,x->x mod 2);
    cochains.s:=ShallowCopy(backend.twists.s); cochains.omega:=ShallowCopy(backend.twists.omega);
    # Both expressions contain the same b faces and lower scalar operations.
    # Evaluate them on each section simplex together, retaining exact signed
    # values while the transport's bounded scalar memo scope is active.
    projected:=ctx.projectMany([psi,correctedM],[0,backend.twists.s]);
    cochains.psi:=List(projected[1],x->x mod 2);
    cochains.M:=projected[2];
    # The proved bar identity d_s M'=2 tau' and integral chain transfer
    # compute the matched value without evaluating the same formula again
    # one simplicial dimension higher.
    cochains.tau:=backend.coboundary(n+3,cochains.M,true);
    if ForAny(cochains.tau,x->x mod 2<>0) then Error("natural matched boundary is not even"); fi;
    cochains.tau:=List(cochains.tau,x->QuoInt(x,2));
    check("d_s tau = 0",n+5,backend.coboundary(n+4,cochains.tau,true),zero(n+5),0);
    check("rho M = psi",n+3,cochains.M,cochains.psi,2);
    check("d_s M = 2 tau",n+4,backend.coboundary(n+3,cochains.M,true),2*cochains.tau,0);
    check("d psi = 0",n+4,backend.coboundary(n+3,cochains.psi,false),zero(n+4),2);
    result.operation:="Tau"; result.noteOperation:="Psi";
    result.operationDegree:=n+3; result.cochain:=cochains.psi;
    # T must lift the same literal simplicial defining system. Coordinate vectors
    # alone lose the homotopy correcting gf on b and on the next cochain c.
    result.modelCochains:=rec(context:=ctx,A:=A,b:=b,psi:=psi,s:=s,omega:=omega);
    return result;
end);

InstallGlobalFunction(koAHSSNaturalOperations,function()
    local operations,callback;
    callback:=function(name)
        return function(ctx)
            local options,result;
            options:=rec(inputType:="mod2");
            if name="Tau" then options.inputType:="integral"; fi;
            result:=koAHSSNaturalSecondary(ctx.backend,ctx.degree,ctx.cochain,options);
            if result.status="candidate" then return result.cochain; fi;
            return rec(status:="unresolved",operations:=[name],reasons:=[result.reason]);
        end;
    end;
    operations:=rec(Tau:=callback("Tau"),Psi:=callback("Psi"),
        T:=koAHSSNaturalTCallback(),useNaturalPrimary:=true,
        usesNaturalT:=true,
        tertiaryReference:="Danus-degree0to3-chi7_tail-epsilon100-eta101-mu0");
    return operations;
end);
