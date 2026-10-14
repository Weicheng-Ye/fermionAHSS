# The higher page differentials d3 (Tau), d4 (Psi) and d5 (T) of the five-row
# ko AHSS, evaluated as the components of the generalized differential of the
# four-cochain stacking model (docs/differential.md): on a flat state
# (A,B,C,D) in package degree k,
#     d_k(A,B,C,D) = (delta_s A, delta B + P(A), delta C + tau_k(A,B), delta_s D + K_k(A,B,C)),
# and Tau_n(A) = tau_k(A,B) with delta B = P(A) (k = n+3), Psi_n(b) = K_k(0,b,c)
# with delta c = D(b) (k = n+2), T_n(A) = K_k(A,B,C) with delta C = tau_k(A,B)
# (k = n+3). The defining cochains are solved on the resolution R, the bar
# cochains Lambda r + H z are formed by the stacking worker with the raw
# comparison homotopy (python/light_rows.py, the page tasks), and the
# components are paired back to R along the comparison chains; K_k is
# delta_s^R of the paired potential. The same
# worker, defining data and potentials serve the extension problem, so a
# tertiary class is never evaluated twice. No separate secondary or tertiary
# formula exists in this package.
#
# Shortcuts: the page engine never asks for a differential with a trivial
# source or target (gap/pages.gi); an even integral input of Tau gives the
# primary rho Dtilde(rho(A/2)); Psi_0 vanishes on a class killed by D; and
# T_0 is the native Pontryagin-square Bockstein of KOAHSS_NativeT0Page.

# A pty ReadLine returns what has arrived so far, which is part of a line
# when the worker has not yet written the rest; read on to the newline.
BindGlobal("KOAHSS_ReadWorkerLine",function(stream)
    local line,more;
    line:=ReadLine(stream);
    while line<>fail and Last(line)<>'\n' do
        more:=ReadLine(stream);
        if more=fail then return fail; fi;
        Append(line,more);
    od;
    return line;
end);

# Page classes of T_0, with trivial sign. Odd positive inputs give zero
# when defined, and multiples of four give a primary Pontryagin Bockstein.
# In our R0 formula b=0 and O=E(c)/2-P(omega)/8. If dc=-d(omega)/2
# modulo two, z=omega+2c is closed modulo four and the same page class is
# -beta_8 P_4(z). The difference is a coboundary and beta_2 Sq^1(rho r),
# r=d(omega)/2, which vanishes because r is an integral cocycle. Changing
# the Z/4 lift adds Dtilde(h), killed in the target page. This uses only
# primary Pontryagin squares on R; it does not choose a universal primitive.
# See the native page-class reduction in docs/tertiary_operations.md.
BindGlobal("KOAHSS_NativeT0Page",function(backend,n,A)
    local W,dW,r,solution,c,z,dz,P,numerator,value,K,modulus,coefficient,
        defining,operation;
    if n<>0 or not ForAll(backend.twists.s,x->x=0) then return fail; fi;
    KOAHSS_CC_CheckVector(A,backend.dimension(0),"native T0 input");
    if IsEmpty(A) or not ForAll(A,x->x=A[1]) then return fail; fi;
    W:=backend.twists.omega; dW:=backend.coboundary(2,W,false);
    if ForAny(dW,x->x mod 2<>0) then Error("native T0 needs a binary omega cocycle"); fi;
    defining:=rec(A:=ShallowCopy(A));
    if IsOddInt(A[1]) then
        if A[1]<0 then return fail; fi;
        solution:=koAHSSSolveCochainEquation(backend,1,W);
        if solution=fail then return fail; fi;
        # On the bar, c=lambda_+ makes O a rational coboundary (R0o).
        # Only its zero page class, not that literal c, is returned here.
        defining.b:=solution.primitive; defining.bSolutions:=solution;
        numerator:=List([1..backend.dimension(4)],j->0);
        modulus:=1; operation:="degree-zero-odd";
    else
        K:=A[1]/2; defining.b:=List([1..backend.dimension(1)],j->0);
        if IsEvenInt(K) then
            # tau_0(4m,0)=0 and E((m mod 2)*omega)=0 on the bar.
            # Thus O=-m P_omega/4 with c=0, a primary operation.
            defining.c:=List([1..backend.dimension(2)],j->0);
            z:=W; dz:=dW; modulus:=4; coefficient:=-K/2;
            operation:="-beta4-P2";
        else
            r:=List(dW,x->x/2);
            solution:=koAHSSSolveCochainEquation(backend,2,List(r,x->(-x) mod 2));
            # Some surviving inputs need a nonzero b. Retain their calibrated
            # calculation rather than treating a failed mod-four lift as zero.
            if solution=fail then return fail; fi;
            c:=solution.primitive; z:=W+2*c; dz:=backend.coboundary(2,z,false);
            if ForAny(dz,x->x mod 4<>0) then Error("native T0 Z/4 lift is not closed"); fi;
            defining.c:=c; defining.cSolutions:=solution;
            defining.omegaLiftModFour:=z;
            modulus:=8; coefficient:=-K; operation:="-beta8-P4";
        fi;
        P:=backend.cupIntegral(0,2,z,2,z,false,false)
            +backend.cupIntegral(1,2,z,3,dz,false,false);
        numerator:=coefficient*P;
    fi;
    value:=backend.coboundary(4,numerator,false);
    if ForAny(value,x->x mod modulus<>0) then
        Error("native T0 phase has a nonintegral boundary");
    fi;
    value:=List(value,x->x/modulus);
    if ForAny(backend.coboundary(5,value,true),x->x<>0) then
        Error("native T0 output is not an integral cocycle");
    fi;
    return rec(status:="computed",cochain:=value,inputDegree:=0,
        phase:=numerator/modulus,phaseNumerator:=numerator,modulus:=modulus,
        definingSystem:=defining,
        referenceConvention:="Danus-degree0to3-chi7_tail-epsilon100-eta101-mu0",
        pageClassOnly:=true,isLocalChoice:=false,isNormalizedOperation:=true,
        transportModel:="native-R",inputGroupBarUsed:=false,modelSamples:=0,
        nativeOperation:=operation,oldCorrectionApplied:=false,muR:=0,
        oddPrimaryCoefficient:=2,
        kernelAudit:=rec(phaseBoundaryIntegral:=true,outputClosed:=true));
end);


# The coboundary matrix of R in degree n as rows (cached on the backend).
BindGlobal("KOAHSS_PageCoboundaryRows",function(backend,n,signed)
    local key,d;
    if not IsBound(backend.coboundaryRows) then backend.coboundaryRows:=rec(); fi;
    key:=Concatenation(String(n),"_",String(signed));
    if not IsBound(backend.coboundaryRows.(key)) then
        d:=backend.dimension(n);
        backend.coboundaryRows.(key):=List([1..d],function(j)
            local e;
            e:=List([1..d],i->0); e[j]:=1;
            return backend.coboundary(n,e,signed);
        end);
    fi;
    return backend.coboundaryRows.(key);
end);

# A binary primitive r on R of a binary coboundary z: delta r = z, or fail.
BindGlobal("KOAHSS_PageBinaryPrimitive",function(backend,n,z)
    local solution;
    z:=List(z,x->x mod 2);
    if backend.dimension(n)=0 then
        if ForAny(z,x->x<>0) then return fail; fi;
        return [];
    fi;
    solution:=koAHSSSolveMod2System(KOAHSS_PageCoboundaryRows(backend,n,false),z);
    if solution=fail then return fail; fi;
    return solution.particular;
end);

# The light task of the stacking worker on the backend's shared model.
BindGlobal("KOAHSS_PageTask",function(backend,k,task,data)
    local model,answer;
    model:=backend.pageModel();
    if model.status<>"computed" then
        return rec(status:="unresolved",reason:=model.reason);
    fi;
    answer:=model.light(k,task,data);
    if answer.status="error" then
        Error("koAHSS: the stacking worker failed the page task ",task,": ",answer.reason);
    fi;
    return answer;
end);

# The class exponents of the images of the generators of H^n(F2) under D,
# and the generator representatives, for the defining-system adjustment.
BindGlobal("KOAHSS_PageDImages",function(backend,n)
    local key,data,gens,reps,target;
    if not IsBound(backend.pageDImages) then backend.pageDImages:=rec(); fi;
    key:=Concatenation("n",String(n+10));
    if not IsBound(backend.pageDImages.(key)) then
        data:=backend.cohomologyData(n,-1);
        gens:=IndependentGeneratorsOfAbelianGroup(data.group);
        reps:=List(gens,g->List(data.represent(g),x->x mod 2));
        target:=backend.cohomologyData(n+2,-1);
        backend.pageDImages.(key):=rec(reps:=reps,target:=target,
            images:=List(reps,r->KOAHSS_PageClassF2(target,backend.nativePrimary("D",n,r))));
    fi;
    return backend.pageDImages.(key);
end);

BindGlobal("KOAHSS_PageClassF2",function(data,v)
    if IsEmpty(IndependentGeneratorsOfAbelianGroup(data.group)) then return []; fi;
    return List(IndependentGeneratorExponents(data.group,data.class(List(v,x->x mod 2))),x->x mod 2);
end);

# The defining system (B_R, C_R) of an integral sign cocycle A of degree n on
# R, in package degree k = n+3: delta B_R = D_R(rho A) and delta C_R = Pi tau_k,
# B_R changed by a cocycle h with [D h] = [Pi tau_k] when needed. The worker
# forms B = Lambda B_R + K_D and C = Lambda C_R + H tau_k. Returns the data
# record of the worker's A tasks, or an unresolved record.
BindGlobal("KOAHSS_PageDefiningSystem",function(backend,n,A,needC)
    local k,answer,BR,v,CR,images,class,solution,l;
    k:=n+3;
    answer:=KOAHSS_PageTask(backend,k,"page_qd",rec(A:=A));
    if answer.status<>"computed" then return answer; fi;
    BR:=KOAHSS_PageBinaryPrimitive(backend,k-2,answer.result.QD);
    if BR=fail then
        return rec(status:="obstructed",reason:="the primary obstruction of the input is nonzero");
    fi;
    answer:=KOAHSS_PageTask(backend,k,"page_tau",rec(A:=A,BR:=BR));
    if answer.status<>"computed" then return answer; fi;
    v:=answer.result.tau;
    if not needC then return rec(status:="computed",A:=A,BR:=BR,tau:=v); fi;
    CR:=KOAHSS_PageBinaryPrimitive(backend,k-1,v);
    if CR=fail then
        # Vanishing in the secondary quotient: change B_R within its class.
        images:=KOAHSS_PageDImages(backend,k-2);
        class:=KOAHSS_PageClassF2(images.target,v);
        if IsEmpty(images.reps) then solution:=fail;
        else solution:=koAHSSSolveMod2System(images.images,class); fi;
        if solution=fail then
            return rec(status:="obstructed",reason:="the input does not survive the secondary quotient");
        fi;
        for l in [1..Length(images.reps)] do
            if solution.particular[l] mod 2<>0 then BR:=BR+images.reps[l]; fi;
        od;
        BR:=List(BR,x->x mod 2);
        answer:=KOAHSS_PageTask(backend,k,"page_tau",rec(A:=A,BR:=BR));
        if answer.status<>"computed" then return answer; fi;
        v:=answer.result.tau;
        CR:=KOAHSS_PageBinaryPrimitive(backend,k-1,v);
        if CR=fail then
            Error("koAHSS: the secondary residual is not a coboundary after the class adjustment");
        fi;
    fi;
    return rec(status:="computed",A:=A,BR:=BR,CR:=CR,tau:=v);
end);

BindGlobal("KOAHSS_PageUnresolved",function(name,answer)
    local reason;
    reason:="unavailable";
    if IsBound(answer.reason) then reason:=answer.reason; fi;
    return rec(status:="unresolved",operations:=[name],
        reasons:=[Concatenation(name,": ",String(reason))]);
end);

# Tau_n(A), the C component of the differential on (A,B,0,0): tau_{n+3}(A,B).
BindGlobal("KOAHSS_PageTau",function(ctx)
    local backend,n,A,system;
    backend:=ctx.backend; n:=ctx.degree; A:=ctx.cochain;
    KOAHSS_CC_CheckVector(A,backend.dimension(n),"Tau input");
    if ForAny(backend.coboundary(n,A,true),x->x<>0) then
        Error("koAHSS: a Tau input must be a sign-integral cocycle");
    fi;
    # An even integral input: the page class is rho Dtilde(rho(A/2)).
    if ForAll(A,x->x mod 2=0) then
        return List(backend.nativePrimary("Dtilde",n,List(A,x->(x/2) mod 2)),x->x mod 2);
    fi;
    system:=KOAHSS_PageDefiningSystem(backend,n,A,false);
    if system.status<>"computed" then return KOAHSS_PageUnresolved("Tau",system); fi;
    return List(system.tau,x->x mod 2);
end);

# Psi_n(b), the D component of the differential on (0,b,c,0): K_{n+2}(0,b,c).
BindGlobal("KOAHSS_PagePsi",function(ctx)
    local backend,n,b,k,answer,c;
    backend:=ctx.backend; n:=ctx.degree; b:=List(ctx.cochain,x->x mod 2);
    KOAHSS_CC_CheckVector(b,backend.dimension(n),"Psi input");
    if ForAny(backend.coboundary(n,b,false),x->x mod 2<>0) then
        Error("koAHSS: a Psi input must be a mod-two cocycle");
    fi;
    # In input degree zero a class killed by D has Psi_0 = 0 on BG.
    if n=0 and IsOne(backend.data(2,-2).class(backend.nativePrimary("D",0,b))) then
        return List([1..backend.dimension(4)],i->0);
    fi;
    k:=n+2;
    answer:=KOAHSS_PageTask(backend,k,"page_qd",rec(b:=b));
    if answer.status<>"computed" then return KOAHSS_PageUnresolved("Psi",answer); fi;
    c:=KOAHSS_PageBinaryPrimitive(backend,k-1,answer.result.QD);
    if c=fail then
        return rec(status:="unresolved",operations:=["Psi"],
            reasons:=["Psi: the primary obstruction of the input is nonzero"]);
    fi;
    # The D component on (0,b,c) in its closed A=0 form (task psi); the
    # general potential (task atom_curvature) gives the same page class at
    # the cost of a theta prism on every comparison simplex.
    answer:=KOAHSS_PageTask(backend,k,"psi",rec(a:=b,c:=c));
    if answer.status<>"computed" then return KOAHSS_PageUnresolved("Psi",answer); fi;
    return answer.result.psi;
end);

# T_n(A), the D component of the differential on the full defining system
# (A,B,C,0): K_{n+3}(A,B,C) = delta_s^R Pi Omega_{n+3}(A,B,C).
BindGlobal("KOAHSS_PageT",function(ctx)
    local backend,n,A,native,system,answer;
    backend:=ctx.backend; n:=ctx.degree; A:=ctx.cochain;
    KOAHSS_CC_CheckVector(A,backend.dimension(n),"T input");
    if ForAny(backend.coboundary(n,A,true),x->x<>0) then
        Error("koAHSS: a T input must be a sign-integral cocycle");
    fi;
    if n<0 or n>3 then
        return rec(status:="unresolved",operations:=["T"],
            reasons:=["T: the tertiary family covers input degrees 0 to 3"]);
    fi;
    native:=KOAHSS_NativeT0Page(backend,n,A);
    if native<>fail then return native.cochain; fi;
    system:=KOAHSS_PageDefiningSystem(backend,n,A,true);
    if system.status<>"computed" then return KOAHSS_PageUnresolved("T",system); fi;
    answer:=KOAHSS_PageTask(backend,n+3,"page_pot",rec(A:=A,BR:=system.BR,CR:=system.CR));
    if answer.status<>"computed" then return KOAHSS_PageUnresolved("T",answer); fi;
    return answer.result.pot;
end);

BindGlobal("KOAHSS_PageOperations",rec(Tau:=KOAHSS_PageTau,Psi:=KOAHSS_PagePsi,T:=KOAHSS_PageT));
