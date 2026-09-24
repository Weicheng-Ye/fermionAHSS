# Sparse native coherent chain operators on an integral HAP resolution.
# This module supplies chain data, not calibrated secondary/tertiary operations.
# A tensor is a list [coefficient, [[degree,basis,R!.elts index],...]].
# Coendomorphism composition uses root-first order and its Koszul sign.
# HAP may leave its final contraction level incomplete: reserve one extra
# resolution degree, so every evaluated total tensor degree is < R length.

InstallGlobalFunction(koAHSSNativeCoherence, function(arg)
    local R, options, length, dimension, index, identityIndex, capacity,
          validate, counters, cache, nextId, owner, reduce, scale, add, act,
          boundary, contract, permute, checkBasis, checkOperation, makeOperation,
          operationBoundary, combination, compose, fillCell, engine,
          identityOp, augmentationOp, diagonals, diagonal, cyclics, cyclic,
          associatorOp, associator, pentagonOp, pentagon, leftUnitOp,
          rightUnitOp, unit, character, evaluateCochains, property;
    if Length(arg)<1 or Length(arg)>2 then
        Error("koAHSSNativeCoherence(resolution[,options])");
    fi;
    R:=arg[1]; options:=rec();
    if Length(arg)=2 then options:=arg[2]; fi;
    if not IsRecord(options) or
       not ForAll(RecNames(options),x->x in ["cacheEntries","validate"]) then
        Error("native coherence options are cacheEntries and validate");
    fi;
    capacity:=256; validate:=true;
    if IsBound(options.cacheEntries) then capacity:=options.cacheEntries; fi;
    if IsBound(options.validate) then validate:=options.validate; fi;
    if not IsInt(capacity) or capacity<0 or not IsBool(validate) then
        Error("native coherence requires nonnegative cacheEntries and boolean validate");
    fi;
    if not IsBoundGlobal("IsHapResolution") or
       not CallFuncList(ValueGlobal("IsHapResolution"),[R]) or
       not IsBound(R!.homotopy) then
        Error("native coherence requires a HAP resolution with contraction");
    fi;
    property:=ValueGlobal("EvaluateProperty");
    length:=property(R,"length");
    if property(R,"characteristic")<>0 or not IsInt(length) or length<1 then
        Error("native coherence requires an integral resolution of known positive length");
    fi;
    dimension:=function(n)
        if not IsInt(n) or n<0 or n>length then
            Error("native coherence resolution degree outside 0..",length);
        fi;
        return R!.dimension(n);
    end;
    if dimension(0)<1 then Error("native coherence requires an augmented vertex"); fi;
    index:=function(g)
        local k;
        k:=Position(R!.elts,g);
        if k=fail then
            if IsBound(R!.appendToElts) then R!.appendToElts(g);
            else Add(R!.elts,g); fi;
            k:=Position(R!.elts,g);
            if k=fail then Error("native coherence could not index group element"); fi;
        fi;
        return k;
    end;
    identityIndex:=index(One(R!.group));
    checkBasis:=function(n,j)
        if not IsInt(j) or j<1 or j>dimension(n) then
            Error("native coherence invalid resolution basis element");
        fi;
    end;
    reduce:=function(terms)
        local result,t,f,k;
        if not IsList(terms) then Error("native tensor must be a list"); fi;
        result:=[];
        for t in terms do
            if not IsList(t) or Length(t)<>2 or not IsInt(t[1]) or
               not IsList(t[2]) then Error("invalid native tensor term"); fi;
            for f in t[2] do
                if not IsList(f) or Length(f)<>3 then Error("invalid native tensor factor"); fi;
                checkBasis(f[1],f[2]);
                if not IsInt(f[3]) or f[3]<1 or f[3]>Length(R!.elts) then
                    Error("invalid native tensor group index");
                fi;
            od;
            if t[1]<>0 then Add(result,[t[1],List(t[2],ShallowCopy)]); fi;
        od;
        Sort(result,function(a,b) return a[2]<b[2]; end);
        terms:=result; result:=[];
        for t in terms do
            k:=Length(result);
            if k>0 and result[k][2]=t[2] then
                result[k][1]:=result[k][1]+t[1];
            else Add(result,t); fi;
        od;
        result:=Filtered(result,t->t[1]<>0);
        MakeImmutable(result);
        return result;
    end;
    scale:=function(c,terms)
        if not IsInt(c) then Error("native tensor coefficient must be integral"); fi;
        return reduce(List(terms,t->[c*t[1],t[2]]));
    end;
    add:=function(a,b) return reduce(Concatenation(a,b)); end;
    act:=function(terms,g)
        return reduce(List(terms,t->[t[1],List(t[2],f->
            [f[1],f[2],index(g*R!.elts[f[3]])])]));
    end;
    boundary:=function(terms)
        local result,t,k,f,u,fs,parity;
        result:=[];
        for t in reduce(terms) do
            parity:=0;
            for k in [1..Length(t[2])] do
                f:=t[2][k];
                if f[1]>0 then
                    for u in R!.boundary(f[1],f[2]) do
                        fs:=List(t[2],ShallowCopy);
                        fs[k]:=[f[1]-1,AbsInt(u[1]),
                            index(R!.elts[f[3]]*R!.elts[u[2]])];
                        Add(result,[t[1]*(-1)^parity*SignInt(u[1]),fs]);
                    od;
                fi;
                parity:=parity+f[1];
            od;
        od;
        return reduce(result);
    end;
    contract:=function(terms)
        local result,t,k,f,u,fs,j;
        result:=[];
        for t in reduce(terms) do
            for k in [1..Length(t[2])] do
                f:=t[2][k];
                if f[1]>=length-1 then
                    Error("native tensor contraction requires resolution degree ",f[1]+2,
                          " including a spare HAP contraction degree; available through ",length);
                fi;
                for u in R!.homotopy(f[1],[f[2],f[3]]) do
                    fs:=List(t[2],ShallowCopy);
                    for j in [1..k-1] do fs[j]:=[0,1,identityIndex]; od;
                    fs[k]:=[f[1]+1,AbsInt(u[1]),u[2]];
                    Add(result,[t[1]*SignInt(u[1]),fs]);
                od;
                if f[1]<>0 then break; fi;
            od;
        od;
        return reduce(result);
    end;
    permute:=function(terms,p)
        local result,t,k,j,parity;
        if not IsList(p) or not ForAll(p,IsInt) or Set(p)<>[1..Length(p)] then
            Error("native tensor permutation must list 1..arity once");
        fi;
        result:=[];
        for t in reduce(terms) do
            if Length(t[2])<>Length(p) then Error("native permutation arity mismatch"); fi;
            parity:=0;
            for k in [1..Length(p)] do
                for j in [k+1..Length(p)] do
                    if p[k]>p[j] then parity:=parity+t[2][p[k]][1]*t[2][p[j]][1]; fi;
                od;
            od;
            Add(result,[t[1]*(-1)^parity,t[2]{p}]);
        od;
        return reduce(result);
    end;
    counters:=rec(requests:=0,hits:=0,misses:=0,createdOperations:=0,
        retainedTerms:=0,peakRetainedTerms:=0,maxEntryTerms:=0);
    cache:=[]; nextId:=0; owner:=[];
    checkOperation:=function(op)
        if not IsRecord(op) or not IsBound(op.arity) or not IsInt(op.arity) or op.arity<0
           or not IsBound(op.degree) or not IsInt(op.degree) or op.degree<0
           or not IsBound(op.evaluate) or not IsFunction(op.evaluate) then
            Error("invalid native operation record");
        fi;
        if IsBound(op.nativeCoherenceOwner) and not IsIdenticalObj(op.nativeCoherenceOwner,owner) then
            Error("native operation belongs to a different coherence engine");
        fi;
    end;
    makeOperation:=function(name,arity,degree,fn)
        local op,id;
        if not IsInt(arity) or arity<0 or not IsInt(degree) or degree<0 or degree>=length then
            Error("native operation arity/degree invalid or resolution too short");
        fi;
        nextId:=nextId+1; id:=nextId;
        counters.createdOperations:=counters.createdOperations+1;
        op:=rec(name:=name,arity:=arity,degree:=degree,nativeCoherenceOwner:=owner);
        op.evaluate:=function(n,j)
            local key,pos,item,value,t;
            checkBasis(n,j);
            if n+degree>=length then
                Error("native operation ",name," needs total tensor degree ",n+degree,
                      " plus one spare HAP contraction degree; resolution available through ",length);
            fi;
            counters.requests:=counters.requests+1;
            key:=[id,n,j]; pos:=PositionProperty(cache,x->x[1]=key);
            if pos<>fail then
                counters.hits:=counters.hits+1;
                item:=Remove(cache,pos); Add(cache,item); return item[2];
            fi;
            counters.misses:=counters.misses+1;
            value:=reduce(fn(n,j));
            if not ForAll(value,t->Length(t[2])=arity and Sum(t[2],f->f[1])=n+degree) then
                Error("native operation produced an incorrect arity or degree");
            fi;
            counters.maxEntryTerms:=Maximum(counters.maxEntryTerms,Length(value));
            if capacity>0 then
                if Length(cache)>=capacity then
                    item:=Remove(cache,1);
                    counters.retainedTerms:=counters.retainedTerms-Length(item[2]);
                fi;
                Add(cache,[key,value]);
                counters.retainedTerms:=counters.retainedTerms+Length(value);
                counters.peakRetainedTerms:=Maximum(counters.peakRetainedTerms,counters.retainedTerms);
            fi;
            return value;
        end;
        return op;
    end;
    operationBoundary:=function(op,n,j)
        local result,u;
        checkOperation(op); checkBasis(n,j);
        result:=boundary(op.evaluate(n,j));
        if n>0 then
            for u in R!.boundary(n,j) do
                result:=add(result,scale(-(-1)^op.degree*SignInt(u[1]),
                    act(op.evaluate(n-1,AbsInt(u[1])),R!.elts[u[2]])));
            od;
        fi;
        return result;
    end;
    combination:=function(terms)
        local arity,degree,t;
        if not IsList(terms) or Length(terms)=0 then Error("native combination must be nonempty"); fi;
        for t in terms do
            if not IsList(t) or Length(t)<>2 or not IsInt(t[1]) then Error("invalid operation summand"); fi;
            checkOperation(t[2]);
        od;
        arity:=terms[1][2].arity; degree:=terms[1][2].degree;
        if not ForAll(terms,t->t[2].arity=arity and t[2].degree=degree) then
            Error("native operation summands must have equal arity and degree");
        fi;
        return makeOperation("linear combination",arity,degree,function(n,j)
            local result,t;
            result:=[];
            for t in terms do result:=add(result,scale(t[1],t[2].evaluate(n,j))); od;
            return result;
        end);
    end;
    compose:=function(f,slot,g)
        checkOperation(f); checkOperation(g);
        if not IsInt(slot) or slot<1 or slot>f.arity then Error("native composition invalid slot"); fi;
        return makeOperation("composition",f.arity+g.arity-1,f.degree+g.degree,function(n,j)
            local result,t,u,v,prefix,suffix,parity;
            result:=[];
            for t in f.evaluate(n,j) do
                v:=t[2][slot];
                prefix:=t[2]{[1..slot-1]}; suffix:=t[2]{[slot+1..Length(t[2])]};
                parity:=f.degree*g.degree+g.degree*Sum(prefix,x->x[1]);
                for u in act(g.evaluate(v[1],v[2]),R!.elts[v[3]]) do
                    Add(result,[(-1)^parity*t[1]*u[1],Concatenation(prefix,u[2],suffix)]);
                od;
            od;
            return result;
        end);
    end;
    fillCell:=function(name,arity,degree,B)
        local op;
        checkOperation(B);
        if not IsString(name) or not IsInt(arity) or arity<0 or
           not IsInt(degree) or degree<1 or B.arity<>arity or B.degree<>degree-1 then
            Error("native cell must have boundary of matching arity and degree minus one");
        fi;
        op:=makeOperation(name,arity,degree,function(n,j)
            local rhs,u,value;
            if validate and operationBoundary(B,n,j)<>[] then
                Error("native cell boundary is not closed: ",name);
            fi;
            rhs:=B.evaluate(n,j);
            if n>0 then
                for u in R!.boundary(n,j) do
                    rhs:=add(rhs,scale((-1)^degree*SignInt(u[1]),
                        act(op.evaluate(n-1,AbsInt(u[1])),R!.elts[u[2]])));
                od;
            fi;
            if boundary(rhs)<>[] then Error("native filler input is not a tensor cycle: ",name); fi;
            if n+degree=1 and Sum(rhs,t->t[1])<>0 then
                Error("native filler degree-zero cycle has nonzero augmentation: ",name);
            fi;
            value:=contract(rhs);
            if boundary(value)<>rhs then Error("native contraction did not fill cell: ",name); fi;
            return value;
        end);
        op.boundaryOperation:=B;
        return op;
    end;
    identityOp:=makeOperation("identity",1,0,function(n,j)
        return [[1,[[n,j,identityIndex]]]];
    end);
    augmentationOp:=makeOperation("augmentation",0,0,function(n,j)
        if n=0 then return [[1,[]]]; fi;
        return [];
    end);
    diagonals:=[];
    diagonal:=function(i)
        local previous,B,op;
        if not IsInt(i) or i<0 or i>=length then
            Error("native cup index must be nonnegative and below resolution length");
        fi;
        if IsBound(diagonals[i+1]) then return diagonals[i+1]; fi;
        if i=0 then
            op:=makeOperation("diagonal 0",2,0,function(n,j)
                local rhs,u,value;
                if n=0 then return [[1,[[0,j,identityIndex],[0,j,identityIndex]]]]; fi;
                rhs:=[];
                for u in R!.boundary(n,j) do
                    rhs:=add(rhs,scale(SignInt(u[1]),
                        act(op.evaluate(n-1,AbsInt(u[1])),R!.elts[u[2]])));
                od;
                value:=contract(rhs);
                if validate and boundary(value)<>rhs then Error("native diagonal contraction failed"); fi;
                return value;
            end);
        else
            previous:=diagonal(i-1);
            B:=makeOperation("binary boundary",2,i-1,function(n,j)
                return add(permute(previous.evaluate(n,j),[2,1]),
                           scale((-1)^i,previous.evaluate(n,j)));
            end);
            op:=fillCell(Concatenation("diagonal ",String(i)),2,i,B);
        fi;
        diagonals[i+1]:=op; return op;
    end;
    associatorOp:=fail;
    associator:=function()
        local m,B;
        if associatorOp=fail then
            m:=diagonal(0);
            B:=combination([[1,compose(m,1,m)],[-1,compose(m,2,m)]]);
            associatorOp:=fillCell("associator",3,1,B);
        fi;
        return associatorOp;
    end;
    pentagonOp:=fail;
    pentagon:=function()
        local m,a,B;
        if pentagonOp=fail then
            m:=diagonal(0); a:=associator();
            B:=combination([[1,compose(m,1,a)],[1,compose(a,2,m)],
                [1,compose(m,2,a)],[-1,compose(a,3,m)],[-1,compose(a,1,m)]]);
            pentagonOp:=fillCell("associativity pentagon",4,2,B);
        fi;
        return pentagonOp;
    end;
    leftUnitOp:=fail; rightUnitOp:=fail;
    unit:=function(slot)
        return fillCell("unit homotopy",1,1,combination([
            [1,compose(diagonal(0),slot,augmentationOp)],[-1,identityOp]]));
    end;
    cyclics:=[];
    cyclic:=function(i)
        local previous,B;
        if not IsInt(i) or i<0 or i>2 then Error("native cyclic index must lie in 0..2"); fi;
        if IsBound(cyclics[i+1]) then return cyclics[i+1]; fi;
        if i=0 then cyclics[1]:=compose(diagonal(0),1,diagonal(0));
        else
            previous:=cyclic(i-1);
            B:=makeOperation("cyclic boundary",3,i-1,function(n,j)
                local value;
                value:=previous.evaluate(n,j);
                if IsOddInt(i) then return add(permute(value,[3,1,2]),scale(-1,value)); fi;
                return add(add(value,permute(value,[3,1,2])),permute(value,[2,3,1]));
            end);
            cyclics[i+1]:=fillCell(Concatenation("cyclic diagonal ",String(i)),3,i,B);
        fi;
        return cyclics[i+1];
    end;
    character:=function(s)
        local j,u;
        if IsInt(s) and s=0 then return g->1; fi;
        if not IsList(s) or Length(s)<>dimension(1) or not ForAll(s,x->x in [0,1]) then
            Error("native sign twist must be a binary degree-one cochain");
        fi;
        if length<2 then Error("native sign cocycle validation requires resolution degree two"); fi;
        for j in [1..dimension(2)] do
            if Sum(R!.boundary(2,j),u->SignInt(u[1])*s[AbsInt(u[1])]) mod 2<>0 then
                Error("native sign twist is not a cocycle");
            fi;
        od;
        s:=ShallowCopy(s);
        return function(g)
            return (-1)^(Sum(R!.homotopy(0,[1,index(g)]),
                u->SignInt(u[1])*s[AbsInt(u[1])]) mod 2);
        end;
    end;
    evaluateCochains:=function(op,degrees,vectors,characters)
        local k,n,result,j,t,value,f,signs,c;
        checkOperation(op);
        if not IsList(degrees) or not IsList(vectors) or not IsList(characters) or
           Length(degrees)<>op.arity or Length(vectors)<>op.arity or Length(characters)<>op.arity then
            Error("native cochain evaluation arity mismatch");
        fi;
        signs:=[];
        for k in [1..op.arity] do
            if not IsList(vectors[k]) or Length(vectors[k])<>dimension(degrees[k]) or
               not ForAll(vectors[k],IsInt) then Error("invalid native integral cochain vector"); fi;
            if IsInt(characters[k]) and characters[k]=0 then Add(signs,g->1);
            elif IsFunction(characters[k]) then Add(signs,characters[k]);
            else Error("native coefficient characters must be functions or zero"); fi;
        od;
        n:=Sum(degrees)-op.degree;
        if n<0 then return []; fi;
        result:=List([1..dimension(n)],j->0);
        for j in [1..Length(result)] do
            for t in op.evaluate(n,j) do
                if List(t[2],f->f[1])=degrees then
                    value:=t[1];
                    for k in [1..op.arity] do
                        f:=t[2][k]; c:=signs[k](R!.elts[f[3]]);
                        if not c in [-1,1] then Error("native sign character must take values +1/-1"); fi;
                        value:=value*c*vectors[k][f[2]];
                    od;
                    result[j]:=result[j]+value;
                fi;
            od;
        od;
        return result;
    end;
    engine:=rec(noBar:=true,calibratedHigherOperations:=false,
        resolution:=R,length:=length,maxTotalDegree:=length-1,dimension:=dimension,
        reduce:=reduce,scale:=scale,add:=add,act:=act,boundary:=boundary,
        contract:=contract,permute:=permute,compose:=compose,
        linearCombination:=combination,fillCell:=fillCell,
        operationBoundary:=operationBoundary,identity:=identityOp,
        augmentation:=augmentationOp,diagonal:=diagonal,associator:=associator,
        pentagon:=pentagon,cyclic:=cyclic,character:=character,
        evaluateCochains:=evaluateCochains);
    engine.leftUnit:=function()
        if leftUnitOp=fail then leftUnitOp:=unit(1); fi;
        return leftUnitOp;
    end;
    engine.rightUnit:=function()
        if rightUnitOp=fail then rightUnitOp:=unit(2); fi;
        return rightUnitOp;
    end;
    # HAP uses D_i|P_n = (-1)^(i*n+i*(i+1)/2) Delta_i|P_n.
    # The seed diagonal agrees, so this converts exact integral pairings.
    engine.hapCupSign:=function(i,n)
        if not IsInt(i) or i<0 or not IsInt(n) or n<0 then Error("invalid HAP cup sign degrees"); fi;
        return (-1)^(i*n+QuoInt(i*(i+1),2));
    end;
    engine.clearCache:=function()
        cache:=[]; counters.retainedTerms:=0;
    end;
    engine.stats:=function()
        local result;
        result:=ShallowCopy(counters);
        result.cachedEntries:=Length(cache); result.cacheEntries:=capacity;
        result.noBar:=true; result.calibratedHigherOperations:=false;
        return result;
    end;
    return engine;
end);
