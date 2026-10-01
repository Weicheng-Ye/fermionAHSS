# Light extension rows (doc/extensions.md, "Light rows").
#
# The light frame of one prime part of one degree computes the rows of the
# relations that the target-layer shortcut and the primary operations do not
# settle, from fixed residue cochains of transported defining data: every
# non-closed defining cochain is P(z;r) = Lambda r + H z with delta r = Pi z
# solved here on R, defining data are selected by linear algebra on the
# cohomology of R (the actual classes, never a zero test of a cochain), and
# the residues are read in the E6 cells. The bar cochains are evaluated by the
# worker of the transferred model (python/extension_light.py); every request
# carries its data, so the worker keeps no frame state.
#
# Markings. A generator of a recorded row denotes one group element: a C
# generator its marked flat state (C_j, D_j); a B generator a marking, an
# integer combination of A=0 atoms (b, c, D) and other B generators; a row
# states what it is exact for (witness.truncatedBelow, witness.shiftedWithin).
# A complete request (KOAHSS_ExtensionPrimeRows) asks for the row of the
# current marking. The first A request builds the A plan: the A-over-D
# relations are adapted (a_i -> a_i + 2^(e_j-e_i) a_j) until their leading B
# vectors are independent, each pivot refers to one atom (or the B basis is
# changed so that it does), the atom is re-marked by the pure-C state of its
# relative C class, and the changed B generators are listed in
# witness.light.remarked, which PrimeRows records again.
#
# env: rec(backend, k, layers, context, part, coboundaryRows(n,signed),
# call(task,data)), where call returns rec(status:="computed",result) or an
# unresolved refusal.
BindGlobal("KOAHSS_LightFrame",function(env)
    local backend,k,layers,context,part,frame,step,call,rowsZ,rowsF2,zeros,
        cochainF2,cochainZ,f2Cache,f2,classF2,zCache,zdata,classZ,imageCache,images,
        solveF2,solveClassZ,combineF2,coords,liftModular,leadC,nB,nC,nD,
        atoms,atomCache,makeAtom,freshAtom,cMarks,cMark,cExact,gaugeData,bGauge,
        pageRow,exactRow,atomRow,EX,DS,ETA,combine,precisionName,rank,satisfies,
        plan,mode,roles,ensurePlan,marks,rows,needs,issued,remarkedSet,singleAtom,
        defaultMarking,ensureMarking,markingRow,generatorRow,switchExact,flatten,
        response,witnessOf,bRequest,cRequest,leading,kerCache,kerColumns,aLower,lowerCache,
        aEntries,aBuilt,buildAPlan,refData,refB,refD,aGauge,adRow,
        adConverted,acRow,aRequest,p3Request,shortcut;
    backend:=env.backend; k:=env.k; layers:=env.layers; context:=env.context; part:=env.part;
    frame:=rec(lastStep:="start",refusal:=fail,shortcuts:=[],tasks:=0);
    step:=function(text) frame.lastStep:=text; end;
    shortcut:=function(name) AddSet(frame.shortcuts,name); end;
    call:=function(task,data)
        local answer;
        step(Concatenation("light worker task ",task));
        frame.tasks:=frame.tasks+1;
        answer:=env.call(task,data);
        if answer.status="error" then
            step(Concatenation("light worker task ",task," failed: ",answer.reason));
            if IsBound(GAPInfo.SystemEnvironment.FERMIONAHSS_LIGHT_DEBUG) then
                Print("#I light worker task ",task," failed:\n",answer.traceback,"\n");
            fi;
            Error("koFull light: worker task ",task," failed: ",answer.reason);
        fi;
        if answer.status<>"computed" then
            frame.refusal:=answer;
            Error("koFull light: the worker refused task ",task,": ",answer.reason);
        fi;
        return answer.result;
    end;
    rowsZ:=n->env.coboundaryRows(n,true);
    rowsF2:=n->env.coboundaryRows(n,false);
    zeros:=n->List([1..n],i->0);
    nB:=Length(layers.B.orders); nC:=Length(layers.C.orders); nD:=Length(layers.D.orders);

    # ------------------------------------------------------------------
    # Linear algebra on R.
    cochainF2:=function(n,rhs)
        local solution;
        rhs:=List(rhs,x->x mod 2);
        if backend.dimension(n)=0 then
            if ForAny(rhs,x->x<>0) then Error("koFull light: a binary source is not a coboundary on R"); fi;
            return [];
        fi;
        solution:=koAHSSSolveMod2System(rowsF2(n),rhs);
        if solution=fail then Error("koFull light: a binary source is not a coboundary on R in degree ",n+1); fi;
        return solution.particular;
    end;
    cochainZ:=function(n,rhs)
        local solution;
        if backend.dimension(n)=0 then
            if ForAny(rhs,x->x<>0) then Error("koFull light: an integral source is not a coboundary on R"); fi;
            return [];
        fi;
        solution:=koAHSSSolveIntegerSystem(rowsZ(n),rhs);
        if solution=fail then Error("koFull light: an integral source is not a coboundary on R in degree ",n+1); fi;
        return solution.particular;
    end;
    f2Cache:=rec();
    f2:=function(n)
        local key,data,gens;
        key:=Concatenation("d",String(n+10));
        if not IsBound(f2Cache.(key)) then
            data:=backend.cohomologyData(n,-1);
            gens:=IndependentGeneratorsOfAbelianGroup(data.group);
            f2Cache.(key):=rec(data:=data,gens:=gens,
                reps:=List(gens,g->List(data.represent(g),x->x mod 2)));
        fi;
        return f2Cache.(key);
    end;
    classF2:=function(n,v)
        local data;
        data:=f2(n);
        if IsEmpty(data.gens) then return []; fi;
        return List(IndependentGeneratorExponents(data.data.group,
            data.data.class(List(v,x->x mod 2))),x->x mod 2);
    end;
    zCache:=rec();
    zdata:=function(n)
        local key,data,gens;
        key:=Concatenation("d",String(n+10));
        if not IsBound(zCache.(key)) then
            data:=backend.cohomologyData(n,0);
            gens:=IndependentGeneratorsOfAbelianGroup(data.group);
            zCache.(key):=rec(data:=data,gens:=gens,reps:=List(gens,data.represent),
                orders:=List(gens,function(g) if Order(g)=infinity then return 0; fi; return Order(g); end));
        fi;
        return zCache.(key);
    end;
    classZ:=function(n,v)
        local data;
        data:=zdata(n);
        if IsEmpty(data.gens) then return []; fi;
        return IndependentGeneratorExponents(data.data.group,data.data.class(v));
    end;
    # Images of the generators of H^n under a primary operation, as class
    # exponents: D and Dbar in H^(n+2)(F2), Dtilde in H^(n+3)(Z_s).
    imageCache:=rec();
    images:=function(name,n)
        local key;
        key:=Concatenation(name,String(n+10));
        if not IsBound(imageCache.(key)) then
            if name="Dbar" then
                imageCache.(key):=List(zdata(n).reps,r->classF2(n+2,backend.nativePrimary("Dbar",n,r)));
            elif name="D" then
                imageCache.(key):=List(f2(n).reps,r->classF2(n+2,backend.nativePrimary("D",n,r)));
            else
                imageCache.(key):=List(f2(n).reps,r->classZ(n+3,backend.nativePrimary("Dtilde",n,r)));
            fi;
        fi;
        return imageCache.(key);
    end;
    solveF2:=function(rows,target)
        local solution;
        if ForAll(target,x->x mod 2=0) then return List(rows,r->0); fi;
        if IsEmpty(rows) then return fail; fi;
        solution:=koAHSSSolveMod2System(rows,List(target,x->x mod 2));
        if solution=fail then return fail; fi;
        return solution.particular;
    end;
    # Integer coefficients x with sum x_i rows_i = target in H^n(Z_s).
    solveClassZ:=function(n,rows,target)
        local orders,relations,j,unit,solution;
        if ForAll(target,x->x=0) then return List(rows,r->0); fi;
        orders:=zdata(n).orders; relations:=[];
        for j in [1..Length(orders)] do
            if orders[j]<>0 then unit:=zeros(Length(orders)); unit[j]:=orders[j]; Add(relations,unit); fi;
        od;
        if IsEmpty(rows) and IsEmpty(relations) then return fail; fi;
        solution:=koAHSSSolveIntegerSystem(Concatenation(rows,relations),target);
        if solution=fail then return fail; fi;
        return solution.particular{[1..Length(rows)]};
    end;
    combineF2:=function(n,coefficients)
        local value,i;
        value:=zeros(backend.dimension(n));
        for i in [1..Length(coefficients)] do
            if coefficients[i] mod 2<>0 then value:=value+f2(n).reps[i]; fi;
        od;
        return List(value,x->x mod 2);
    end;
    # E6 cell coordinates of a cocycle of row q in the layer name.
    coords:=function(name,n,q,v)
        local layer,element;
        layer:=layers.(name);
        if IsEmpty(layer.orders) then return []; fi;
        if q in [-1,-2] then v:=List(v,x->x mod 2); fi;
        element:=layer.cell.project(backend.cohomologyData(n,q).class(v));
        return KOAHSS_ExtensionGroupCoordinates(layer.group,layer.generators,element);
    end;
    # An integral cocycle t with rho_m t in the class of v (coefficient
    # reduction, light-completion-transport Lemma 4); determined modulo m.
    liftModular:=function(n,v,m)
        local vt,y,t;
        vt:=List(v,x->x mod m);
        y:=backend.coboundary(n,vt,true);
        if ForAny(y,x->x mod m<>0) then Error("koFull light: a residue is not closed modulo ",m); fi;
        t:=vt-m*cochainZ(n,List(y,x->x/m));
        if ForAny(backend.coboundary(n,t,true),x->x<>0) then Error("koFull light: the integral lift is not closed"); fi;
        return t;
    end;
    # The C coordinates of (Sq^1+s)b, the leading part of a B cocycle.
    leadC:=function(b)
        local twice;
        if nC=0 then return []; fi;
        twice:=backend.coboundary(k-2,List(b,x->x mod 2),true);
        if ForAny(twice,IsOddInt) then Error("koFull light: a B cocycle has an odd Bockstein numerator"); fi;
        return coords("C",k-1,-2,List(twice,x->(x/2) mod 2));
    end;

    # ------------------------------------------------------------------
    # Flat A=0 atoms (0,b,c,D) and marked C states (C_j,D_j).
    atoms:=[]; atomCache:=rec();
    # c=fail: solve delta c = Pi Q_D(b) and correct c by the Dtilde image so
    # that the curvature class vanishes (flat-admissible); otherwise c is a
    # re-marked admissible choice. D solves delta_s D = -Pi J_k(0,b,c).
    makeAtom:=function(b,c)
        local J,class,solution,l;
        b:=List(b,x->x mod 2);
        step("B atom");
        if c=fail then
            c:=cochainF2(k-1,call("qd",rec(b:=b)).QD);
            J:=call("atom_curvature",rec(b:=b,c:=c)).J;
            class:=classZ(k+2,J);
            if ForAny(class,x->x<>0) then
                solution:=solveClassZ(k+2,images("Dtilde",k-1),-class);
                if solution=fail then Error("koFull light: a B cocycle has no flat-admissible C"); fi;
                for l in [1..Length(solution)] do
                    if solution[l] mod 2<>0 then c:=c+f2(k-1).reps[l]; fi;
                od;
                c:=List(c,x->x mod 2);
                J:=call("atom_curvature",rec(b:=b,c:=c)).J;
            fi;
        else
            c:=List(c,x->x mod 2);
            J:=call("atom_curvature",rec(b:=b,c:=c)).J;
        fi;
        Add(atoms,rec(b:=b,c:=c,D:=cochainZ(k+1,-J),L:=leadC(b),rows:=rec()));
        return Length(atoms);
    end;
    freshAtom:=function(b)
        local key;
        b:=List(b,x->x mod 2); key:=Concatenation("b",String(b));
        if not IsBound(atomCache.(key)) then atomCache.(key):=makeAtom(b,fail); fi;
        return atomCache.(key);
    end;
    cMarks:=[];
    cMark:=function(j)
        local c,r;
        if not IsBound(cMarks[j]) then
            step("C marking");
            c:=List(layers.C.cochains[j],x->x mod 2);
            r:=call("c_mark",rec(c:=c));
            cMarks[j]:=rec(c:=c,D:=cochainZ(k+1,-r.J),gamma:=r.gamma);
        fi;
        return cMarks[j];
    end;
    # The exact row of the marked C state: its double is the D state 2D_c+gamma(c,c).
    cExact:=function(j)
        local mark,t;
        mark:=cMark(j); t:=2*mark.D+mark.gamma;
        if ForAny(backend.coboundary(k+1,t,true),x->x<>0) then Error("koFull light: an exact C row is not closed"); fi;
        return coords("D",k+1,-4,t);
    end;

    # ------------------------------------------------------------------
    # B over D: gauge selection and residues.
    gaugeData:=function(gauge)
        local data;
        data:=ShallowCopy(gauge); Unbind(data.kind); return data;
    end;
    # The gauge (u,y,pi) with delta pi = x + C + P (light-transport-proof (8)):
    # pure, a D gauge y0, or a Tau gauge u whose E3 class the page Tau map
    # selects; the remaining class is the actual one, killed by D.
    bGauge:=function(atom,eps)
        local data,v,class,cell3,v3,solution,hom,u0,uR,yR,w,gauge;
        step("B gauge");
        data:=rec(b:=atom.b,Cref:=List(eps,j->cMark(j).c),gauge:=rec(),want:=["xC"]);
        v:=call("gauge",data).xC;
        class:=backend.cohomologyData(k-1,-2).class(v);
        gauge:=rec(kind:="pure");
        if not IsOne(class) then
            v3:=fail;
            if k>=4 then
                cell3:=context.getCell(3,k-1,-2);
                if KOAHSS_IsUnresolved(cell3) then Error("koFull light: the E3 C cell is unresolved"); fi;
                v3:=cell3.project(class);
            fi;
            if v3=fail or IsOne(v3) then
                solution:=solveF2(images("D",k-3),classF2(k-1,v));
                if solution=fail then Error("koFull light: a B gauge class is not in the image of D"); fi;
                gauge:=rec(kind:="D",y0:=combineF2(k-3,solution));
            else
                hom:=context.getMap(3,k-4,0);
                if IsRecord(hom) then Error("koFull light: the page Tau map is unavailable"); fi;
                u0:=PreImagesRepresentative(hom,v3);
                if u0=fail then Error("koFull light: a B gauge class is not in the image of Tau"); fi;
                uR:=backend.cohomologyData(k-4,0).represent(context.getCell(3,k-4,0).lift(u0));
                data.gauge:=rec(u:=uR); data.want:=["QDu"];
                yR:=cochainF2(k-3,call("gauge",data).QDu);
                data.gauge:=rec(u:=uR,yR:=yR); data.want:=["xCp"];
                w:=call("gauge",data).xCp;
                solution:=solveF2(images("D",k-3),classF2(k-1,w));
                if solution=fail then Error("koFull light: a Tau gauge leaves a class outside the image of D"); fi;
                gauge:=rec(kind:="tau",u:=uR,yR:=yR,y1:=combineF2(k-3,solution));
                shortcut("PF-B-tau");
            fi;
            data.gauge:=gaugeData(gauge); data.want:=["target"];
            v:=call("gauge",data).target;
        fi;
        gauge.pi:=cochainF2(k-2,v);
        return gauge;
    end;
    EX:=rec(eta:=false,level:=2); DS:=rec(eta:=false,level:=1);
    ETA:=rec(eta:=true,level:=2);
    combine:=function(x,y) return rec(eta:=x.eta or y.eta,level:=Minimum(x.level,y.level)); end;
    precisionName:=function(x)
        if x.level=0 then return "TR"; fi;
        if x.eta then if x.level=2 then return "ETA"; fi; return "ETA-D"; fi;
        if x.level=2 then return "EX"; fi; return "DS";
    end;
    rank:=need->Position(["basis","lower","element"],need)-1;
    satisfies:=function(precision,need) return precision.level>=rank(need); end;
    # Page form rho[T_b] = [Phi_P(b,P;c,pi)] of a kernel atom (light-completion-tau
    # (9), theta_k = 0): exact for (b,c) with some D completion.
    pageRow:=function(atom)
        local gauge,phi,t;
        gauge:=bGauge(atom,[]);
        step("B over D page form");
        phi:=call("bd_page",rec(b:=atom.b,c:=atom.c,gauge:=gaugeData(gauge))).Phi;
        t:=liftModular(k+1,phi,2);
        shortcut("PF-B");
        return rec(entries:=rec(C:=zeros(nC),D:=coords("D",k+1,-4,t)),precision:=DS,
            method:=Concatenation("B/D page form Phi_P, ",gauge.kind," gauge [Tau 9]"));
    end;
    # The exact residue 2D_b + Z_B (light-transport-proof (11)-(12)) relative
    # to the pure-C reference of the leading C part, with the marked C states.
    exactRow:=function(atom)
        local eps,gauge,Z,T,j;
        eps:=Filtered([1..Length(atom.L)],j->atom.L[j] mod 2<>0);
        gauge:=bGauge(atom,eps);
        step("B over D exact residue");
        Z:=call("bd_exact",rec(b:=atom.b,c:=atom.c,Cref:=List(eps,j->cMark(j).c),
            gauge:=gaugeData(gauge))).Z;
        T:=2*atom.D+Z;
        for j in eps do T:=T-cMark(j).D; od;
        if ForAny(backend.coboundary(k+1,T,true),x->x<>0) then Error("koFull light: an exact B residue is not closed"); fi;
        return rec(entries:=rec(C:=List(atom.L,x->x mod 2),D:=coords("D",k+1,-4,T)),precision:=EX,
            method:="B/D exact residue 2D_b + Z_B [P 11-12]");
    end;
    atomRow:=function(index,need)
        local atom;
        atom:=atoms[index];
        if k>=4 and ForAll(atom.L,x->x mod 2=0) and rank(need)<=1 then
            if not IsBound(atom.rows.page) then atom.rows.page:=pageRow(atom); fi;
            return atom.rows.page;
        fi;
        if not IsBound(atom.rows.exact) then atom.rows.exact:=exactRow(atom); fi;
        return atom.rows.exact;
    end;

    # ------------------------------------------------------------------
    # The B plan and the markings of the B generators.
    plan:=fail; mode:="exact"; roles:=[];
    marks:=[]; rows:=[]; needs:=[]; issued:=[]; remarkedSet:=[];
    ensurePlan:=function()
        local L,i,pivots,pivotRows,solution;
        if plan<>fail then return plan; fi;
        step("B plan");
        L:=[];
        for i in [1..nB] do L[i]:=List(leadC(layers.B.cochains[i]),x->x mod 2); od;
        # Absorption (light-transport-proof (27)-(28)) in its strong regime:
        # every two-primary or free D generator has order two.
        if KOAHSS_LightAbsorptionEnabled() and ForAll(layers.D.orders,
               o->o=2 or (o<>0 and KOAHSS_RelationPrime(o)<>2)) then
            mode:="strong";
        fi;
        pivots:=[]; pivotRows:=[];
        for i in [1..nB] do
            if ForAll(L[i],x->x=0) then roles[i]:=rec(role:="kernel");
            else
                solution:=solveF2(pivotRows,L[i]);
                if solution=fail then
                    Add(pivots,i); Add(pivotRows,L[i]); roles[i]:=rec(role:="pivot");
                else
                    roles[i]:=rec(role:="nonpivot",P:=pivots{Filtered([1..Length(solution)],
                        j->solution[j] mod 2<>0)});
                fi;
            fi;
        od;
        plan:=rec(L:=L,pivots:=pivots);
        return plan;
    end;
    # The default marking of a B generator: its own atom, and in the strong
    # regime the formal sum (kernel atom of b_i + sum b_p) + sum p of a
    # non-pivot. The atom is made only when a row or a reference needs it.
    defaultMarking:=function(i)
        local beta,q;
        if mode="strong" and roles[i].role="nonpivot" then
            beta:=List(layers.B.cochains[i],x->x mod 2);
            for q in roles[i].P do beta:=beta+layers.B.cochains[q]; od;
            return Concatenation([[1,"atom",freshAtom(beta)]],List(roles[i].P,q->[1,"gen",q]));
        fi;
        return [[1,"atom",freshAtom(layers.B.cochains[i])]];
    end;
    ensureMarking:=function(i)
        ensurePlan();
        if not IsBound(marks[i]) then marks[i]:=rec(terms:=fail,version:=1); fi;
        if marks[i].terms=fail then marks[i].terms:=defaultMarking(i); fi;
        return marks[i];
    end;
    singleAtom:=function(i)
        ensurePlan();
        if not IsBound(marks[i]) or marks[i].terms=fail then
            return not (mode="strong" and roles[i].role="nonpivot");
        fi;
        return Length(marks[i].terms)=1 and marks[i].terms[1][1]=1 and marks[i].terms[1][2]="atom";
    end;
    # The row of a marking: an absorbed pivot is (L,0) in the eta frame; a
    # formal marking sums the rows of its terms.
    markingRow:=function(i,need)
        local mark,entries,precision,term,r,method;
        if not IsBound(marks[i]) then marks[i]:=rec(terms:=fail,version:=1); fi;
        if mode="strong" and roles[i].role="pivot" and singleAtom(i) then
            shortcut("absorption-strong");
            return rec(entries:=rec(C:=ShallowCopy(plan.L[i]),D:=zeros(nD)),precision:=ETA,
                method:="B/D absorbed pivot (L,0) in the eta frame [P 27-28]");
        fi;
        mark:=ensureMarking(i);
        entries:=rec(C:=zeros(nC),D:=zeros(nD)); precision:=EX; method:=fail;
        for term in mark.terms do
            if term[2]="atom" then r:=atomRow(term[3],need);
            else r:=generatorRow(term[3],"basis"); fi;
            entries.C:=entries.C+term[1]*r.entries.C;
            entries.D:=entries.D+term[1]*r.entries.D;
            precision:=combine(precision,r.precision);
            if method=fail then method:=r.method; fi;
        od;
        if Length(mark.terms)>1 or mark.terms[1][1]<>1 then
            method:="B/D formal combination of atom and generator rows [P 27-29]";
        fi;
        return rec(entries:=entries,precision:=precision,method:=method);
    end;
    generatorRow:=function(i,need)
        local mark,r;
        ensurePlan();
        if not IsBound(marks[i]) then marks[i]:=rec(terms:=fail,version:=1); fi;
        mark:=marks[i];
        if IsBound(rows[i]) and rows[i].version=mark.version and satisfies(rows[i].precision,need) then
            return rows[i];
        fi;
        r:=markingRow(i,need); r.version:=mark.version; rows[i]:=r;
        if IsBound(issued[i]) then AddSet(remarkedSet,i); fi;
        return r;
    end;
    switchExact:=function()
        local i;
        if mode="exact" then return; fi;
        mode:="exact";
        for i in [1..nB] do
            if IsBound(marks[i]) then
                marks[i]:=rec(terms:=fail,version:=marks[i].version+1);
                if IsBound(issued[i]) then AddSet(remarkedSet,i); fi;
            fi;
        od;
    end;
    # The atoms of a marking with their integer coefficients (generator terms expanded).
    flatten:=function(i,coefficient)
        local result,term,inner,found,entry;
        result:=[];
        for term in ensureMarking(i).terms do
            if term[2]="atom" then inner:=[[coefficient*term[1],term[3]]];
            else inner:=flatten(term[3],coefficient*term[1]); fi;
            for entry in inner do
                found:=First(result,x->x[2]=entry[2]);
                if found=fail then Add(result,ShallowCopy(entry));
                else found[1]:=found[1]+entry[1]; fi;
            od;
        od;
        return Filtered(result,x->x[1]<>0);
    end;

    # ------------------------------------------------------------------
    # Responses.
    response:=function(lower,entries,witness)
        local row,name,stored,j;
        row:=zeros(lower.generatorCount);
        for name in ["D","C","B"] do
            if IsBound(entries.(name)) then
                stored:=First(lower.layers,l->l.name=name);
                for j in [1..Length(entries.(name))] do
                    row[stored.startColumn+j-1]:=entries.(name)[j];
                od;
            fi;
        od;
        return rec(status:="computed",lowerPresentationId:=lower.presentationId,
            lowerCoordinates:=row,witness:=witness);
    end;
    witnessOf:=function(order,prime,method,precision,target,light)
        local witness;
        witness:=rec(operation:="xtimes",power:=order,prime:=prime,model:="light-R",method:=method,
            sufficiency:=target,light:=light);
        light.precision:=precisionName(precision);
        light.shortcuts:=ShallowCopy(frame.shortcuts);
        light.tasks:=frame.tasks;
        if precision.eta then witness.truncatedBelow:="D"; witness.measuredLayers:=["C"];
        else witness.truncatedBelow:=fail; witness.measuredLayers:=["C","D"]; fi;
        if precision.level=2 then witness.shiftedWithin:=fail;
        elif precision.level=1 then witness.shiftedWithin:="D"; fi;
        return witness;
    end;
    bRequest:=function(index,lower,target,complete,through)
        local need,r;
        ensurePlan();
        need:="lower";
        if complete then
            if IsBound(needs[index]) then need:=needs[index];
            elif through>=3 then need:="element"; fi;
        fi;
        r:=generatorRow(index,need);
        issued[index]:=true; RemoveSet(remarkedSet,index);
        return response(lower,r.entries,witnessOf(2,2,r.method,r.precision,target,
            rec(rowType:=roles[index].role,mode:=mode,version:=r.version)));
    end;
    cRequest:=function(index,lower,target)
        local witness;
        witness:=witnessOf(layers.C.orders[index],2,"C/D exact row 2D_c + gamma_C(c,c) [R 11.7]",EX,target,
            rec(rowType:="C-exact"));
        witness.measuredLayers:=["D"];
        return response(lower,rec(D:=cExact(index)),witness);
    end;

    # ------------------------------------------------------------------
    # A relations: lower systems, gauges, plans and residues.
    # The B coordinates of D(rho u), delta_s u = m A: the leading part.
    leading:=function(A,m)
        local u;
        u:=cochainZ(k-4,m*A);
        return coords("B",k-2,-1,backend.nativePrimary("D",k-4,List(u,x->x mod 2)));
    end;
    # The A=0 systems (b_i,c_i) of a basis of ker D on H^(k-2)(F2), with
    # their curvature classes J_i in H^(k+2)(Z_s) (light-transport-proof (3)).
    kerCache:=fail;
    kerColumns:=function()
        local imgs,null,w,b,c;
        if kerCache<>fail then return kerCache; fi;
        imgs:=images("D",k-2); kerCache:=[];
        if IsEmpty(imgs) then return kerCache; fi;
        if IsEmpty(imgs[1]) then null:=IdentityMat(Length(imgs));
        else null:=List(NullspaceMat(List(imgs,r->List(r,x->x*Z(2)))),IntVecFFE); fi;
        for w in null do
            b:=combineF2(k-2,w);
            c:=cochainF2(k-1,call("qd",rec(b:=b)).QD);
            Add(kerCache,rec(b:=b,c:=c,J:=classZ(k+2,call("atom_curvature",rec(b:=b,c:=c)).J)));
        od;
        return kerCache;
    end;
    # A flat-admissible lower system (A,B_A',C_A') (light-transport-proof
    # Section 2.1): B by Q_D(rho A) and the D image, C by f#(A,B), and the
    # curvature class killed by the A=0 systems and the Dtilde image.
    # Without needC only B_A' is needed (it is B_A when ker D is zero).
    lowerCache:=rec();
    aLower:=function(A,needC)
        local system,v,solution,jA,cols,n,l,terms,zs,key;
        key:=Concatenation("A",String(A),String(needC));
        if IsBound(lowerCache.(key)) then return lowerCache.(key); fi;
        step("A lower system");
        system:=rec(A:=A);
        system.BR:=cochainF2(k-2,call("a_step",rec(A:=A,want:=["QDa"])).QDa);
        v:=call("a_step",rec(A:=A,BR:=system.BR,want:=["fsharp"])).fsharp;
        solution:=solveF2(images("D",k-2),classF2(k,v));
        if solution=fail then Error("koFull light: f#(A,B) is not killed by the D image"); fi;
        if ForAny(solution,x->x mod 2<>0) then
            system.BR:=List(system.BR+combineF2(k-2,solution),x->x mod 2);
            v:=call("a_step",rec(A:=A,BR:=system.BR,want:=["fsharp"])).fsharp;
        fi;
        cols:=kerColumns();
        lowerCache.(key):=system;
        if not needC and IsEmpty(cols) then return system; fi;
        system.CR:=cochainF2(k-1,v);
        jA:=classZ(k+2,call("a_step",rec(A:=A,BR:=system.BR,CR:=system.CR,want:=["pot"])).pot);
        solution:=solveClassZ(k+2,Concatenation(List(cols,c->c.J),images("Dtilde",k-1)),-jA);
        if solution=fail then Error("koFull light: no flat-admissible lower system of an A relation"); fi;
        n:=solution{[1..Length(cols)]}; l:=solution{[Length(cols)+1..Length(solution)]};
        terms:=List(Filtered([1..Length(cols)],i->n[i]<>0),i->[n[i],cols[i].b,cols[i].c]);
        zs:=List(Filtered([1..Length(l)],j->l[j] mod 2<>0),j->f2(k-1).reps[j]);
        if not IsEmpty(terms) or not IsEmpty(zs) then system.star:=rec(terms:=terms,z:=zs); fi;
        return system;
    end;
    # The worker data, the B cocycle and the D marking of a reference.
    refData:=function(ref)
        if ref.kind="zero" then return rec(kind:="zero"); fi;
        if ref.kind="atom" then return rec(kind:="atom",b:=atoms[ref.atom].b,c:=atoms[ref.atom].c); fi;
        if ref.kind="pureC" then return rec(kind:="pureC",C:=List(ref.eps,j->cMark(j).c)); fi;
        return rec(kind:="lower",terms:=List(ref.terms,t->[t[1],atoms[t[2]].b,atoms[t[2]].c]));
    end;
    refB:=function(ref)
        local value,t;
        value:=zeros(backend.dimension(k-2));
        if ref.kind="atom" then value:=ShallowCopy(atoms[ref.atom].b);
        elif ref.kind="lower" then
            for t in ref.terms do if t[1] mod 2<>0 then value:=value+atoms[t[2]].b; fi; od;
        fi;
        return List(value,x->x mod 2);
    end;
    refD:=function(ref)
        local value,j;
        value:=zeros(backend.dimension(k+1));
        if ref.kind="atom" then value:=ShallowCopy(atoms[ref.atom].D);
        elif ref.kind="pureC" then for j in ref.eps do value:=value+cMark(j).D; od; fi;
        return value;
    end;
    # The relation gauge U, Y of (relative note (4.2), (11.8)): delta_s U = mA,
    # U corrected by the Dbar image so that Q_D(rho U)+X_B+B0 is exact, and Y.
    aGauge:=function(x,system,ref)
        local data,u,ell,solution,g;
        step("A relation gauge");
        data:=ShallowCopy(system);
        data.m:=x.m; data.e:=x.e; data.ref:=refData(ref);
        u:=cochainZ(k-4,x.m*x.A);
        ell:=classF2(k-2,backend.nativePrimary("D",k-4,List(u,v->v mod 2)))+classF2(k-2,refB(ref));
        solution:=solveF2(images("Dbar",k-4),ell);
        if solution=fail then Error("koFull light: the leading class of an A relation is not in the Dbar image"); fi;
        for g in [1..Length(solution)] do
            if solution[g] mod 2<>0 then u:=u+zdata(k-4).reps[g]; fi;
        od;
        data.U:=rec(u:=u);
        return data;
    end;
    # The relation row of an adapted A-over-D relation a' (light-transport-proof
    # (14) with (25)), in B, C and D coordinates.
    adRow:=function(x)
        local system,data,v,class,cell3,v3,hom,u0,vR,solution,t,entries;
        if IsBound(x.row) then return x.row; fi;
        system:=aLower(x.A,true);
        data:=aGauge(x,system,x.ref);
        data.want:=["Ysrc"];
        data.U.YR:=cochainF2(k-3,call("a_step",data).Ysrc);
        data.want:=["RC"];
        v:=call("a_step",data).RC;
        class:=backend.cohomologyData(k-1,-2).class(v);
        if not IsOne(class) then
            step("A over D final gauge");
            cell3:=context.getCell(3,k-1,-2);
            if KOAHSS_IsUnresolved(cell3) then Error("koFull light: the E3 C cell is unresolved"); fi;
            v3:=cell3.project(class);
            if not IsOne(v3) then
                hom:=context.getMap(3,k-4,0);
                if IsRecord(hom) then Error("koFull light: the page Tau map is unavailable"); fi;
                u0:=PreImagesRepresentative(hom,v3);
                if u0=fail then Error("koFull light: the relative C class of an A relation is not killed"); fi;
                vR:=backend.cohomologyData(k-4,0).represent(context.getCell(3,k-4,0).lift(u0));
                data.U.tau:=rec(v:=vR,yR:=[]);
                data.want:=["QDv"];
                data.U.tau.yR:=cochainF2(k-3,call("a_step",data).QDv);
                data.want:=["RC"];
                v:=call("a_step",data).RC;
            fi;
            solution:=solveF2(images("D",k-3),classF2(k-1,v));
            if solution=fail then Error("koFull light: the relative C class of an A relation is not in the D image"); fi;
            if ForAny(solution,y->y mod 2<>0) then
                data.U.yD:=combineF2(k-3,solution);
                v:=call("a_step",data).RC;
            fi;
        fi;
        data.U.WR:=cochainF2(k-2,v);
        step("A over D residue");
        data.want:=["residue"];
        v:=call("a_step",data).residue;
        t:=liftModular(k+1,List(v-refD(x.ref),y->y mod x.m),x.m);
        entries:=rec(B:=ShallowCopy(x.Lp),C:=zeros(nC),D:=coords("D",k+1,-4,t));
        if x.ref.kind="pureC" then for v in x.ref.eps do entries.C[v]:=1; od; fi;
        if k=6 then shortcut("unary-q"); fi;
        x.row:=entries;
        return entries;
    end;
    # The row of an original A generator: m_x a_x = m_x a'_x - sum_l m_l a'_l.
    adConverted:=function(x)
        local entries,l,r,name;
        entries:=ShallowCopy(adRow(x));
        for name in ["B","C","D"] do entries.(name):=ShallowCopy(entries.(name)); od;
        for l in x.adapt do
            r:=adRow(aEntries[l]);
            for name in ["B","C","D"] do entries.(name):=entries.(name)-r.(name); od;
        od;
        return entries;
    end;
    aEntries:=[]; aBuilt:=false;
    # The A plan (light-transport-proof Section 5, (29); relative note Section 8).
    buildAPlan:=function(lower)
        local i,x,AD,pivots,pv,L,A,coef,simple,tier2,J,M,Minv,q,eps,Q,Labs,
            data,system,v,j,terms,ok;
        step("A plan");
        ensurePlan();
        for i in part.generators.A do
            x:=rec(index:=i,m:=layers.A.orders[i],A:=layers.A.cochains[i]);
            x.e:=Log2Int(x.m);
            x.target:=CallFuncList(ValueGlobal("KOAHSS_ExtensionTargetLayer"),[lower,x.m]);
            x.T:=x.target.layer;
            if x.T in ["C","D"] then x.ell:=List(leading(x.A,x.m),y->y mod 2); fi;
            aEntries[i]:=x;
        od;
        AD:=Filtered(List(part.generators.A,i->aEntries[i]),x->x.T="D");
        AD:=Concatenation(List(Reversed(Set(List(AD,x->x.e))),e->Filtered(AD,x->x.e=e)));
        pivots:=[];
        for x in AD do
            L:=ShallowCopy(x.ell); A:=ShallowCopy(x.A); coef:=[];
            for pv in pivots do
                if L[pv.col]<>0 then
                    L:=List(L+pv.Lp,y->y mod 2); A:=A+2^(pv.e-x.e)*pv.A; Add(coef,pv.index);
                fi;
            od;
            x.Lp:=L; x.A:=A; x.adapt:=coef; x.col:=fail;
            if ForAny(L,y->y<>0) then x.col:=PositionProperty(L,y->y<>0); Add(pivots,x); fi;
        od;
        # A pivot refers to the atom of its one B generator; otherwise (tier
        # two: several B generators, or a formal marking) to a fresh atom after a
        # unitriangular change of the B basis, which is done without absorption.
        simple:=x->Number(x.Lp,y->y<>0)=1 and singleAtom(x.col);
        repeat
            ok:=true;
            tier2:=ForAny(pivots,x->not simple(x));
            if tier2 and mode="strong" then
                switchExact(); tier2:=ForAny(pivots,x->not simple(x));
            fi;
            # References and preliminary relative C classes.
            for x in AD do
                if x.col=fail then x.ref:=rec(kind:="zero");
                elif not tier2 then x.ref:=rec(kind:="atom",atom:=flatten(x.col,1)[1][2]);
                else
                    v:=zeros(backend.dimension(k-2));
                    for j in [1..nB] do
                        if x.Lp[j]<>0 then v:=v+layers.B.cochains[j]; fi;
                    od;
                    x.ref:=rec(kind:="atom",atom:=freshAtom(v));
                fi;
                system:=aLower(x.A,true);
                data:=aGauge(x,system,x.ref);
                data.want:=["Ysrc"]; data.U.YR:=cochainF2(k-3,call("a_step",data).Ysrc);
                data.want:=["RC"];
                x.q:=List(coords("C",k-1,-2,call("a_step",data).RC),y->y mod 2);
            od;
            # Pure-C references must keep the eta frame of the absorbed pivots
            # (light-transport-proof (28)): otherwise the B rows are not absorbed.
            Q:=List(Filtered(AD,x->x.col=fail and ForAny(x.q,y->y<>0)),x->x.q);
            if mode="strong" and not IsEmpty(Q) then
                Labs:=List(Filtered([1..nB],i->roles[i].role="pivot"),i->plan.L[i]);
                if not IsEmpty(Labs) and RankMat(List(Concatenation(Q,Labs),r->r*Z(2)))
                     <>RankMat(List(Q,r->r*Z(2)))+RankMat(List(Labs,r->r*Z(2))) then
                    switchExact(); ok:=false;
                fi;
            fi;
        until ok;
        for x in AD do
            if x.col=fail then
                eps:=Filtered([1..nC],j->x.q[j]<>0);
                if IsEmpty(eps) then x.ref:=rec(kind:="zero"); else x.ref:=rec(kind:="pureC",eps:=eps); fi;
            elif ForAny(x.q,y->y<>0) then
                # Re-mark the reference atom by the pure-C state of its relative C class.
                v:=ShallowCopy(atoms[x.ref.atom].c);
                for j in [1..nC] do if x.q[j]<>0 then v:=v+cMark(j).c; fi; od;
                x.ref:=rec(kind:="atom",atom:=makeAtom(atoms[x.ref.atom].b,v));
                shortcut("re-marking");
            fi;
        od;
        # New markings of the B generators that the pivots refer to.
        if not tier2 then
            for x in pivots do
                if marks[x.col].terms[1][3]<>x.ref.atom then
                    marks[x.col]:=rec(terms:=[[1,"atom",x.ref.atom]],version:=marks[x.col].version+1);
                fi;
                needs[x.col]:="element";
                generatorRow(x.col,"element");
                AddSet(remarkedSet,x.col);
            od;
        else
            # The basis: the pivot atoms, then the other columns J. Its matrix is
            # unitriangular, so the B generators of the pivot columns are integer
            # combinations of the pivot atoms and of the generators in J, and each
            # pivot atom is exactly the combination sum_i Lp_i b_i.
            J:=Filtered([1..nB],j->ForAll(pivots,x->x.col<>j));
            M:=Concatenation(List(pivots,x->ShallowCopy(x.Lp)),List(J,function(j)
                local r; r:=zeros(nB); r[j]:=1; return r; end));
            Minv:=Inverse(M);
            if not ForAll(Minv,r->ForAll(r,IsInt)) then Error("koFull light: the B basis change is not unimodular"); fi;
            for x in pivots do
                terms:=[];
                for q in [1..Length(pivots)] do
                    if Minv[x.col][q]<>0 then Add(terms,[Minv[x.col][q],"atom",pivots[q].ref.atom]); fi;
                od;
                for q in [1..Length(J)] do
                    if Minv[x.col][Length(pivots)+q]<>0 then
                        Add(terms,[Minv[x.col][Length(pivots)+q],"gen",J[q]]);
                        if not IsBound(needs[J[q]]) then needs[J[q]]:="basis"; fi;
                        generatorRow(J[q],"lower");
                        AddSet(remarkedSet,J[q]);
                    fi;
                od;
                marks[x.col]:=rec(terms:=terms,version:=marks[x.col].version+1);
                needs[x.col]:="basis";
            od;
            for x in pivots do
                generatorRow(x.col,"element");
                AddSet(remarkedSet,x.col);
            od;
            shortcut("tier-2");
        fi;
        aBuilt:=true;
    end;
    # A over C (target C): t_C relative to the lower part of the reference
    # sum_i ell_i b_i (light-transport-proof (10); page form PF-AC2 and
    # Theorem C where they apply).
    acRow:=function(x,lower)
        local terms,i,t,found,ref,system,data,v,form;
        terms:=[];
        for i in [1..nB] do
            if x.ell[i]<>0 then
                if not IsBound(needs[i]) or rank(needs[i])<1 then needs[i]:="lower"; fi;
                for t in flatten(i,1) do
                    found:=First(terms,y->y[2]=t[2]);
                    if found=fail then Add(terms,ShallowCopy(t)); else found[1]:=found[1]+t[1]; fi;
                od;
            fi;
        od;
        terms:=Filtered(terms,t->t[1]<>0);
        if IsEmpty(terms) then ref:=rec(kind:="zero"); else ref:=rec(kind:="lower",terms:=terms); fi;
        system:=rec(A:=x.A);
        if x.m=2 then system:=aLower(x.A,false); fi;
        data:=aGauge(x,system,ref);
        form:="exact";
        if x.m=2 and ref.kind="zero" and k in [5,6] then form:="page";
        elif x.m=2 and ref.kind="zero" and k=4 and ForAll(backend.twists.omega,w->w mod 2=0) then form:="theoremC";
        fi;
        step(Concatenation("A over C (",form,")"));
        if form="theoremC" then
            data.want:=["thmC"]; v:=call("a_step",data).thmC; shortcut("theorem-C");
        elif form="page" then
            data.want:=["ypp"]; data.U.ypp:=cochainF2(k-3,call("a_step",data).ypp);
            data.want:=["page"]; v:=call("a_step",data).page; shortcut("PF-AC2");
        else
            data.want:=["Ysrc"]; data.U.YR:=cochainF2(k-3,call("a_step",data).Ysrc);
            data.want:=["RC"]; v:=call("a_step",data).RC;
        fi;
        return rec(entries:=rec(B:=ShallowCopy(x.ell),C:=coords("C",k-1,-2,v)),
            method:=rec(exact:="A/C exact R_C [P 10]",page:="A/C page form Phi' + omega a [AS AC-5]",
                theoremC:="A/C Theorem C: hD(B_A,B_A) [AS AC-6]").(form));
    end;
    aRequest:=function(index,order,lower,target)
        local x,entries,witness,light,method,precision,list;
        if not aBuilt then buildAPlan(lower); fi;
        x:=aEntries[index];
        if x.T<>target.layer then Error("koFull light: an A target layer changed"); fi;
        if target.layer="C" then
            entries:=acRow(x,lower); method:=entries.method; entries:=entries.entries;
            precision:=rec(eta:=true,level:=2);
        else
            entries:=adConverted(x);
            method:="A/D residue Psi_rel [P 14] with the bounded tail [P 25]";
            if k=6 then method:=Concatenation(method,", unary q [P 20]"); fi;
            precision:=EX;
        fi;
        list:=List(remarkedSet,i->["B",i]); remarkedSet:=[];
        light:=rec(rowType:=Concatenation("A over ",target.layer),mode:=mode);
        if not IsEmpty(list) then light.remarked:=list; fi;
        witness:=witnessOf(order,2,method,precision,target,light);
        Unbind(witness.shiftedWithin);
        if target.layer="C" then
            witness.measuredLayers:=["B","C"]; light.precision:="through C";
        else witness.measuredLayers:=["B","C","D"]; fi;
        return response(lower,entries,witness);
    end;
    # The prime three in degrees five and six: t_D = 2*3^(e-1) Y with
    # rho_3 Y = P^1_s rho_3 A (light-transport-proof (31)).
    p3Request:=function(index,order,lower,target)
        local A,e,Y,witness;
        A:=layers.A.cochains[index]; e:=Length(Factors(order));
        step("p=3 row");
        if k=5 then Y:=call("p3_power",rec(A:=A)).Y;
        else Y:=liftModular(k+1,call("p3_power",rec(A:=A)).P1,3); fi;
        witness:=rec(operation:="xtimes",power:=order,prime:=3,model:="light-R",
            method:="p=3: 2*3^(e-1) Y, rho_3 Y = P^1 rho_3 A [P 31]",
            measuredLayers:=["D"],truncatedBelow:=fail,sufficiency:=target,
            light:=rec(rowType:="p3",shortcuts:=[],tasks:=frame.tasks));
        return response(lower,rec(D:=coords("D",k+1,-4,2*3^(e-1)*Y)),witness);
    end;
    frame.row:=function(name,index,order,lower,target,options)
        local through;
        through:=3; if IsBound(options.through) then through:=options.through; fi;
        if part.prime=3 then return p3Request(index,order,lower,target); fi;
        if name="C" then return cRequest(index,lower,target); fi;
        if name="B" then
            return bRequest(index,lower,target,IsBound(options.complete) and options.complete=true,through);
        fi;
        return aRequest(index,order,lower,target);
    end;
    return frame;
end);
