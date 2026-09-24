# Finite diagrams of tertiary defining systems. All equations are explicit
# integral equations, including witnesses for the specified quotient classes.
# Solving a finite diagram is not certification of a universal stable class.
BindGlobal("KOAHSS_TE_Matrix",function(backend,degree)
    local d,e,rows,i,v,row;
    d:=backend.dimension(degree); e:=backend.dimension(degree+1); rows:=[];
    for i in [1..d] do
        v:=List([1..d],j->0); v[i]:=1;
        row:=backend.coboundary(degree,v,true);
        KOAHSS_CC_CheckVector(row,e,"tertiary coboundary matrix row");
        Add(rows,row);
    od;
    return rows;
end);

BindGlobal("KOAHSS_TE_Secondary",function(backend,n,a,options,b)
    local opts,evaluate;
    opts:=rec(inputType:="mod2");
    if IsBound(options.secondaryOptions) then opts:=ShallowCopy(options.secondaryOptions); fi;
    opts.inputType:="mod2";
    if b<>fail then opts.b:=b; fi;
    evaluate:=koAHSSSecondaryReference;
    if IsBound(backend.usesNaturalPrimary) and backend.usesNaturalPrimary then
        evaluate:=koAHSSNaturalSecondary;
    fi;
    if IsBound(options.psiEvaluator) then evaluate:=options.psiEvaluator; fi;
    if not IsFunction(evaluate) then Error("koAHSS: psiEvaluator must be a function"); fi;
    return evaluate(backend,n,a,opts);
end);

# Rows generating integral boundaries, J images, and (when requested) Psi
# images. Witness coefficients are unrestricted integers, not binary phases.
BindGlobal("KOAHSS_TE_Denominator",function(backend,n,quotient,options)
    local rows,H,K,generators,g,matrix,kernel,weights,h,i,secondary;
    rows:=KOAHSS_TE_Matrix(backend,n+4);
    if quotient="cohomology" then return rec(status:="ok",rows:=rows); fi;
    if not quotient in ["J","full"] then
        Error("koAHSS: tertiary quotient must be cohomology, J, or full");
    fi;
    if not IsBound(backend.data) or not IsBound(backend.primary) then
        return rec(status:="unavailable",reason:="quotient equations need cohomology representatives and primary operations");
    fi;
    H:=backend.data(n+2,-1);
    for g in GeneratorsOfGroup(H.group) do
        Add(rows,backend.primary("Dtilde",n+2,H.represent(g)));
    od;
    if quotient="J" then return rec(status:="ok",rows:=rows); fi;
    H:=backend.data(n+1,-1); K:=backend.data(n+3,-1);
    generators:=GeneratorsOfGroup(H.group);
    matrix:=List(generators,g->Exponents(K.class(backend.primary("D",n+1,H.represent(g)))));
    kernel:=koAHSSSolveMod2System(matrix,List(GeneratorsOfGroup(K.group),g->0));
    for weights in kernel.homogeneousGenerators do
        h:=List([1..backend.dimension(n+1)],j->0);
        for i in [1..Length(weights)] do
            if weights[i]=1 then h:=h+H.represent(generators[i]); fi;
        od;
        h:=List(h,x->x mod 2);
        secondary:=KOAHSS_TE_Secondary(backend,n+1,h,options,fail);
        if not IsRecord(secondary) or not IsBound(secondary.status)
           or secondary.status<>"candidate" then
            return rec(status:="unavailable",reason:="Psi image needed for the full tertiary quotient",
                       secondary:=secondary);
        fi;
        Add(rows,secondary.cochain);
    od;
    return rec(status:="ok",rows:=rows);
end);

InstallGlobalFunction(koAHSSTertiaryEquations,function(arg)
    local nodes,relations,options,modulus,builder,prepared,i,j,node,previous,
          backend,n,A,b,c,psi,localFamily,delta,phase,output,shared,preparedNode,
          relation,type,from,to,first,second,third,h,z,value,secondary,quotient,
          denominatorCache,denominator,addClassEquation,failure,terms,zero,
          identity,matrix,lowerMatrix,phaseWitness,liftWitness,equations,
          index,checkNode,cochainCheck,closedCheck,counts,autoBoundary;
    if not Length(arg) in [1,2,3] then
        Error("koAHSSTertiaryEquations(nodes[,relations[,options]])");
    fi;
    nodes:=arg[1]; relations:=[]; options:=rec();
    if Length(arg)>=2 then relations:=arg[2]; fi;
    if Length(arg)=3 then options:=arg[3]; fi;
    if not IsList(nodes) or IsEmpty(nodes) or not IsList(relations) or not IsRecord(options) then
        Error("koAHSS: tertiary equations need nonempty nodes, a relation list, and options");
    fi;
    modulus:=4;
    if IsBound(options.modulus) then modulus:=options.modulus; fi;
    if not IsInt(modulus) or modulus<2 or modulus mod 2<>0 then
        Error("koAHSS: tertiary equation modulus must be positive and even");
    fi;
    autoBoundary:=true;
    if IsBound(options.autoBoundary) then
        if not IsBool(options.autoBoundary) then Error("koAHSS: autoBoundary must be boolean"); fi;
        autoBoundary:=options.autoBoundary;
    fi;
    builder:=KOAHSS_IntegralEquationBuilder(); prepared:=[];
    counts:=rec(residual:=0,sharedR:=0,boundary:=0,first:=0,last:=0,
                naturality:=0,stability:=0,additivity:=0,representative:=0,calibration:=0,phase:=0);
    cochainCheck:=function(model,degree,v,label)
        KOAHSS_CC_CheckVector(v,model.dimension(degree),label);
    end;
    closedCheck:=function(model,degree,v)
        if ForAny(model.coboundary(degree,v,true),x->x<>0) then
            Error("koAHSS: a prescribed tertiary comparison value is not an integral cocycle");
        fi;
    end;
    for i in [1..Length(nodes)] do
        node:=nodes[i];
        if not IsRecord(node) or not ForAll(["backend","n","A","b","c","psi"],
                                              name->IsBound(node.(name))) then
            Error("koAHSS: each tertiary node needs backend,n,A,b,c,psi");
        fi;
        backend:=node.backend; n:=node.n;
        KOAHSS_EquationBackend(backend,n);
        A:=ShallowCopy(node.A); b:=ShallowCopy(node.b);
        c:=ShallowCopy(node.c); psi:=ShallowCopy(node.psi);
        cochainCheck(backend,n,A,"tertiary A");
        cochainCheck(backend,n+1,b,"tertiary b");
        cochainCheck(backend,n+2,c,"tertiary c");
        cochainCheck(backend,n+3,psi,"tertiary psi");
        if not ForAll(Concatenation(b,c,psi),x->x in [0,1]) then
            Error("koAHSS: tertiary b,c,psi must be binary");
        fi;
        if not IsBound(backend.primary) or not IsFunction(backend.primary) then
            Error("koAHSS: tertiary nodes require the primary operation backend");
        fi;
        if List(backend.coboundary(n+1,b,false),x->x mod 2)<>backend.primary("Dbar",n,A) then
            Error("koAHSS: tertiary b must satisfy db=Dbar(A)");
        fi;
        # Reuse the verified upper-phase identity and the explicit k=1 term.
        localFamily:=koAHSSTertiaryCandidates(backend,n,psi,c,modulus,A);
        if localFamily=fail then
            return rec(status:="unavailable",node:=i,reason:="signed cups needed for the coefficient-one correction");
        fi;
        delta:=KOAHSS_TE_Matrix(backend,n+4); shared:=fail;
        for j in [1..Length(prepared)] do
            previous:=prepared[j];
            if IsIdenticalObj(previous.backend,backend) and previous.n=n
               and previous.A=A and previous.b=b then
                if previous.psi<>psi then
                    Error("koAHSS: equal (A,b) must use the same secondary cochain psi");
                fi;
                shared:=j; break;
            fi;
        od;
        if shared=fail then phase:=builder.variable(Concatenation("R",String(i)),backend.dimension(n+4));
        else phase:=prepared[shared].phaseVariable; counts.sharedR:=counts.sharedR+1; fi;
        output:=builder.variable(Concatenation("Tref",String(i)),backend.dimension(n+5));
        builder.equation(Concatenation("residual/integrality node ",String(i)),
            [[phase,delta],[output,-modulus*IdentityMat(output.size)]],
            -backend.coboundary(n+4,localFamily.upperPhase,true));
        counts.residual:=counts.residual+1;
        preparedNode:=rec(backend:=backend,n:=n,A:=A,b:=b,c:=c,psi:=psi,
            phaseVariable:=phase,outputVariable:=output,upperPhase:=localFamily.upperPhase,
            correction:=localFamily.correction,sharedWith:=shared);
        Add(prepared,preparedNode);
    od;
    checkNode:=function(number)
        if not IsInt(number) or not number in [1..Length(prepared)] then
            Error("koAHSS: tertiary relation has an invalid node index");
        fi;
        return prepared[number];
    end;
    denominatorCache:=[]; failure:=fail;
    denominator:=function(target,mode)
        local cached,result;
        for cached in denominatorCache do
            if IsIdenticalObj(cached.backend,target.backend) and cached.n=target.n and cached.mode=mode then
                return cached.result;
            fi;
        od;
        result:=KOAHSS_TE_Denominator(target.backend,target.n,mode,options);
        Add(denominatorCache,rec(backend:=target.backend,n:=target.n,mode:=mode,result:=result));
        return result;
    end;
    # Entries [node index, integer coefficient, cochain matrix] describe a
    # combination of final T values. The fixed correction enters the RHS.
    addClassEquation:=function(label,entries,target,value,mode)
        local rows,witness,parts,rhs,entry,source,termMatrix,known;
        cochainCheck(target.backend,target.n+5,value,"tertiary comparison value");
        closedCheck(target.backend,target.n+5,value);
        rows:=denominator(target,mode);
        if rows.status<>"ok" then failure:=rows; return; fi;
        parts:=[]; rhs:=ShallowCopy(value);
        for entry in entries do
            source:=checkNode(entry[1]); termMatrix:=entry[2]*entry[3];
            if not IsList(entry[3]) or Length(entry[3])<>source.outputVariable.size
               or not ForAll(entry[3],r->IsList(r) and Length(r)=target.outputVariable.size and ForAll(r,IsInt)) then
                Error("koAHSS: tertiary comparison map has incorrect dimensions");
            fi;
            Add(parts,[source.outputVariable,termMatrix]);
            known:=KOAHSS_CC_Multiply(source.correction,termMatrix,target.outputVariable.size);
            rhs:=rhs-known;
        od;
        witness:=builder.variable(Concatenation("witness ",label),Length(rows.rows));
        Add(parts,[witness,List(rows.rows,r->-r)]);
        builder.equation(label,parts,rhs);
    end;
    # The boundary is an exact cohomology class, before either tertiary
    # quotient. A specified boundary callback may supply its normalized value.
    if autoBoundary then
        for i in [1..Length(prepared)] do
            node:=prepared[i];
            if ForAll(node.A,x->x=0) then
                if IsBound(nodes[i].boundaryValue) then value:=nodes[i].boundaryValue;
                else
                    secondary:=KOAHSS_TE_Secondary(node.backend,node.n+1,node.b,options,node.c);
                    if not IsRecord(secondary) or not IsBound(secondary.status) or secondary.status<>"candidate" then
                        return rec(status:="unavailable",reason:="the exact A=0 secondary boundary could not be evaluated",
                                   node:=i,secondary:=secondary);
                    fi;
                    value:=secondary.cochain;
                fi;
                addClassEquation(Concatenation("A=0 boundary node ",String(i)),
                    [[i,1,IdentityMat(node.outputVariable.size)]],node,value,"cohomology");
                counts.boundary:=counts.boundary+1;
            fi;
        od;
    fi;
    # Shared R automatically enforces independence of c at the phase level;
    # include the resulting last-nullhomotopy law as an independently checked
    # cohomology equation as well.
    relations:=ShallowCopy(relations);
    for i in [1..Length(prepared)] do
        if prepared[i].sharedWith<>fail then
            Add(relations,rec(type:="last",from:=prepared[i].sharedWith,to:=i));
        fi;
    od;
    for index in [1..Length(relations)] do
        relation:=relations[index];
        if not IsRecord(relation) or not IsBound(relation.type) then
            Error("koAHSS: tertiary relations need a type");
        fi;
        type:=relation.type;
        if not type in RecNames(counts) or type in ["residual","sharedR"] then
            Error("koAHSS: unknown tertiary relation type");
        fi;
        counts.(type):=counts.(type)+1;
        if type in ["boundary","calibration","phase"] then
            if not IsBound(relation.node) or not IsBound(relation.value) then
                Error("koAHSS: boundary/calibration/phase relations need node and value");
            fi;
            node:=checkNode(relation.node);
            if type="phase" then
                cochainCheck(node.backend,node.n+4,relation.value,"prescribed R numerator");
                liftWitness:=builder.variable(Concatenation("phase lift ",String(index)),node.phaseVariable.size);
                builder.equation(Concatenation("phase value ",String(index)),
                    [[node.phaseVariable,IdentityMat(node.phaseVariable.size)],
                     [liftWitness,-modulus*IdentityMat(node.phaseVariable.size)]],relation.value);
            else
                quotient:="cohomology";
                if type="calibration" then quotient:="full"; fi;
                if IsBound(relation.quotient) then quotient:=relation.quotient; fi;
                if type="boundary" and quotient<>"cohomology" then
                    Error("koAHSS: a secondary boundary is an exact cohomology class");
                fi;
                if type="boundary" and ForAny(node.A,x->x<>0) then
                    Error("koAHSS: secondary boundary constraints require A=0");
                fi;
                addClassEquation(Concatenation(type," ",String(index)),
                    [[relation.node,1,IdentityMat(node.outputVariable.size)]],node,relation.value,quotient);
            fi;
        elif type="additivity" then
            if not ForAll(["left","right","sum"],name->IsBound(relation.(name))) then
                Error("koAHSS: additivity needs left,right,sum node indices");
            fi;
            first:=checkNode(relation.left); second:=checkNode(relation.right); node:=checkNode(relation.sum);
            if not IsIdenticalObj(first.backend,node.backend) or not IsIdenticalObj(second.backend,node.backend)
               or first.n<>node.n or second.n<>node.n or first.A+second.A<>node.A then
                Error("koAHSS: additivity requires the same model and A_sum=A_left+A_right");
            fi;
            identity:=IdentityMat(node.outputVariable.size);
            addClassEquation(Concatenation("additivity ",String(index)),
                [[relation.sum,1,identity],[relation.left,-1,identity],[relation.right,-1,identity]],
                node,List([1..node.outputVariable.size],i->0),"full");
        else
            if not IsBound(relation.from) or not IsBound(relation.to) then
                Error("koAHSS: tertiary variation/map relations need from and to");
            fi;
            from:=relation.from; to:=relation.to; first:=checkNode(from); node:=checkNode(to);
            identity:=IdentityMat(node.outputVariable.size);
            value:=List([1..node.outputVariable.size],i->0); quotient:="cohomology";
            if type in ["first","last","representative"] then
                if not IsIdenticalObj(first.backend,node.backend) or first.n<>node.n then
                    Error("koAHSS: defining-choice relations must use the same model and input degree");
                fi;
                matrix:=identity;
                if type in ["first","last"] and first.A<>node.A then
                    Error("koAHSS: nullhomotopy changes require the same integral input A");
                fi;
                if type="last" then
                    if first.b<>node.b then Error("koAHSS: last-nullhomotopy changes require the same b"); fi;
                    z:=List(node.c-first.c,x->x mod 2);
                    value:=node.backend.primary("Dtilde",node.n+2,z);
                elif type="first" then
                    h:=List(node.b-first.b,x->x mod 2); quotient:="J";
                    if koAHSSSolveCochainEquation(node.backend,node.n+2,
                           node.backend.primary("D",node.n+1,h))=fail then
                        Error("koAHSS: first-nullhomotopy comparison requires D[h]=0");
                    fi;
                    if IsBound(relation.value) then value:=relation.value;
                    else
                        secondary:=KOAHSS_TE_Secondary(node.backend,node.n+1,h,options,fail);
                        if not IsRecord(secondary) or not IsBound(secondary.status) or secondary.status<>"candidate" then
                            return rec(status:="unavailable",reason:="first-nullhomotopy Psi value could not be evaluated",
                                       relation:=index,secondary:=secondary);
                        fi;
                        value:=secondary.cochain;
                    fi;
                else
                    quotient:="full";
                    if node.n=0 then
                        if node.A<>first.A then Error("koAHSS: distinct degree-zero cocycles are not cohomologous"); fi;
                    elif koAHSSSolveIntegerSystem(KOAHSS_TE_Matrix(node.backend,node.n-1),node.A-first.A)=fail then
                        Error("koAHSS: representative relation inputs are not cohomologous");
                    fi;
                fi;
            else
                # A finite naturality/suspension diagram supplies its actual
                # cochain maps. Check commutation with d in the relevant degrees.
                if not IsBound(relation.matrix) or not IsBound(relation.lowerMatrix) then
                    Error("koAHSS: naturality/stability need output and lower cochain matrices");
                fi;
                if type="naturality" and first.n<>node.n then
                    Error("koAHSS: naturality preserves the input degree");
                fi;
                matrix:=relation.matrix; lowerMatrix:=relation.lowerMatrix;
                if not IsList(matrix) or Length(matrix)<>first.outputVariable.size
                   or not ForAll(matrix,r->IsList(r) and Length(r)=node.outputVariable.size and ForAll(r,IsInt))
                   or not IsList(lowerMatrix) or Length(lowerMatrix)<>first.phaseVariable.size
                   or not ForAll(lowerMatrix,r->IsList(r) and Length(r)=node.phaseVariable.size and ForAll(r,IsInt)) then
                    Error("koAHSS: cochain map has incorrect dimensions");
                fi;
                for i in [1..first.phaseVariable.size] do
                    if KOAHSS_CC_Multiply(KOAHSS_TE_Matrix(first.backend,first.n+4)[i],matrix,node.outputVariable.size)
                       <>KOAHSS_CC_Multiply(lowerMatrix[i],KOAHSS_TE_Matrix(node.backend,node.n+4),node.outputVariable.size) then
                        Error("koAHSS: supplied cochain map does not commute with the coboundary");
                    fi;
                od;
                quotient:="full";
                # Optional coherent pullback of R, allowing a phase coboundary.
                if IsBound(relation.phaseCompatible) and not IsBool(relation.phaseCompatible) then
                    Error("koAHSS: phaseCompatible must be boolean");
                fi;
                if IsBound(relation.phaseCompatible) and relation.phaseCompatible=true then
                    phaseWitness:=builder.variable(Concatenation("phase homotopy ",String(index)),node.backend.dimension(node.n+3));
                    liftWitness:=builder.variable(Concatenation("phase integer lift ",String(index)),node.phaseVariable.size);
                    builder.equation(Concatenation("R map coherence ",String(index)),
                        [[node.phaseVariable,IdentityMat(node.phaseVariable.size)],
                         [first.phaseVariable,-lowerMatrix],
                         [phaseWitness,List(KOAHSS_TE_Matrix(node.backend,node.n+3),r->-r)],
                         [liftWitness,-modulus*IdentityMat(node.phaseVariable.size)]],
                        List([1..node.phaseVariable.size],j->0));
                fi;
            fi;
            addClassEquation(Concatenation(type," ",String(index)),
                [[to,1,identity],[from,-1,matrix]],node,value,quotient);
        fi;
        if failure<>fail then return failure; fi;
    od;
    if failure<>fail then return failure; fi;
    equations:=builder.assemble(); equations.status:="equations";
    equations.nodes:=prepared; equations.modulus:=modulus; equations.constraints:=counts;
    equations.autoBoundary:=autoBoundary; equations.tertiaryCoefficient:=1;
    equations.isFiniteDiagram:=true; equations.isNormalizedOperation:=false;
    equations.assumptions:=[
        "supplied psi values are the selected common lower cochains",
        "the supplied diagram maps describe the intended inputs, twists, and suspensions",
        "Psi values on a basis generate its image only when the selected Psi is additive modulo J"];
    return equations;
end);

InstallGlobalFunction(koAHSSSolveTertiaryEquations,function(problem)
    local solution,result,decode,particular,homogeneous,verify,item,family;
    if not IsRecord(problem) or not IsBound(problem.status) or problem.status<>"equations"
       or not IsBound(problem.nodes) or not IsBound(problem.modulus) then
        Error("koAHSS: solve a problem returned by koAHSSTertiaryEquations");
    fi;
    solution:=koAHSSSolveIntegerSystem(problem.matrix,problem.rhs);
    if solution=fail then
        return rec(status:="inconsistent",problem:=problem,
                   reason:="the supplied finite compatibility and normalization equations have no solution at this modulus");
    fi;
    decode:=function(vector,constant)
        local outputs,node,r,t,part;
        outputs:=[];
        for node in problem.nodes do
            r:=vector{[node.phaseVariable.offset+1..node.phaseVariable.offset+node.phaseVariable.size]};
            t:=vector{[node.outputVariable.offset+1..node.outputVariable.offset+node.outputVariable.size]};
            if constant then t:=t+node.correction; fi;
            part:=rec(phaseNumerator:=r,phase:=List(r,x->x mod problem.modulus),T:=t);
            if ForAny(node.backend.coboundary(node.n+5,t,true),x->x<>0) then
                Error("koAHSS: solved tertiary value is not a closed integral cochain");
            fi;
            Add(outputs,part);
        od;
        return outputs;
    end;
    particular:=decode(solution.particular,true);
    homogeneous:=List(solution.homogeneousGenerators,v->decode(v,false));
    result:=rec(status:="solved",problem:=problem,solutionLattice:=solution,
        particular:=particular,homogeneousGenerators:=homogeneous,
        modulus:=problem.modulus,tertiaryCoefficient:=1,
        isFiniteDiagram:=true,isNormalizedOperation:=false);
    return result;
end);
