# Cochain comparison identity, including nonclosed inputs and the cell bar.
CallFuncList(function()
    local patterns,i,n,a,b,word,evaluator,groupOrder,R,backend,native,tr,compare,
        lift,K,p,q,da,db,source,cup,vertices,sigma,lhs,rhs,cells,engine,
        labels,integral,binary,reduceEncoded,chain,limited;
    patterns:=[];
    for i in [0..1] do
        patterns[i+1]:=[];
        word:=List([0..i+1],j->1+j mod 2);
        for n in [0..3] do
            patterns[i+1][n+1]:=[];
            for a in Filtered(Combinations([1..n+1]),x->not IsEmpty(x)) do
                for b in Filtered(Combinations([1..n+1]),x->not IsEmpty(x)) do
                    if Length(a)+Length(b)<>n+i+2 then continue; fi;
                    evaluator:=KOAHSS_NaturalWordEvaluator(word,[Length(a)-1,Length(b)-1],
                        [function(x) if x=a then return 1; else return 0; fi; end,
                         function(x) if x=b then return 1; else return 0; fi; end]);
                    if evaluator([1..n+1])=1 then Add(patterns[i+1][n+1],[a,b]); fi;
                od;
            od;
        od;
    od;
    for groupOrder in [3,4] do
        R:=ResolutionFiniteGroup(CyclicGroup(groupOrder),7);
        backend:=koAHSSHAPSpace(R,koAHSSNaturalOperations()).koAHSS(0,0,4);
        native:=backend.nativeCoherence();
        for cells in [false,true] do
            if cells then tr:=KOAHSS_ExtensionCellTransport(R,5);
            else tr:=backend.naturalBar(); fi;
            compare:=KOAHSS_ExtensionPrimaryTransport(tr,native,2000000);
            lift:=function(v,sigma) return Sum(tr.f(sigma),t->t[3]*v[t[1]]) mod 2; end;
            K:=function(i,p,a,q,b,sigma)
                local answer;
                if i<0 then return 0; fi;
                answer:=compare(i,sigma,patterns); Assert(0,answer.status="computed");
                return Sum(Filtered(answer.terms,t->t[2]=p and t[4]=q),
                    t->a[t[3]+1]*b[t[5]+1]) mod 2;
            end;
            vertices:=AsList(R!.group);
            if cells then vertices:=List(vertices,g->[g,1]); fi;
            for i in [0..1] do
                p:=i+1; q:=p; n:=p+q-i; a:=[1]; b:=[1];
                da:=List(backend.coboundary(p,a,false),x->x mod 2);
                db:=List(backend.coboundary(q,b,false),x->x mod 2);
                source:=backend.cupMod2(i,p,a,q,b);
                word:=List([0..i+1],j->1+j mod 2);
                cup:=KOAHSS_NaturalWordEvaluator(word,[p,q],[x->lift(a,x),x->lift(b,x)]);
                for sigma in Tuples(vertices,n+1) do
                    lhs:=Sum([1..n+1],j->K(i,p,a,q,b,sigma{Filtered([1..n+1],h->h<>j)}))
                        +K(i,p+1,da,q,b,sigma)+K(i,p,a,q+1,db,sigma);
                    rhs:=lift(source,sigma)+cup(sigma)
                        +K(i-1,p,a,q,b,sigma)+K(i-1,q,b,p,a,sigma);
                    Assert(0,(lhs-rhs) mod 2=0);
                od;
            od;
        od;
    od;
    limited:=KOAHSS_ExtensionPrimaryTransport(tr,native,1);
    limited:=limited(0,vertices{[1,2]},patterns);
    Assert(0,limited.status="unresolved" and limited.code="primary-comparison-term-budget");
    # h2 is the same fixed normalized chain modulo two, including the
    # normalization of its output vertices and the signed input fibers.
    reduceEncoded:=function(terms)
        local keys;
        keys:=Set(List(terms,t->t[3]));
        return Filtered(List(keys,v->[Sum(Filtered(terms,t->t[3]=v),t->t[1]) mod 2,v]),t->t[1]<>0);
    end;
    backend:=koAHSSHAPSpace(R,koAHSSNaturalOperations()).koAHSS([1],0,4);
    engine:=KOAHSS_ExtensionNormalizedTransport(backend,3);
    chain:=engine.chain(2,0); Assert(0,chain.status="computed");
    for labels in List(chain.terms,t->t[3]) do
        integral:=engine.terms("h",labels); binary:=engine.terms("h2",labels);
        Assert(0,integral.status="computed" and binary.status="computed");
        Assert(0,reduceEncoded(integral.terms)=reduceEncoded(binary.terms));
    od;
end,[]);
