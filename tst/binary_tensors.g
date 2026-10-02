# Compare the binary engine to the independent integral coherent diagonals,
# before evaluating cochains. This checks group actions as well as pairings.
CallFuncList(function()
    local groups,G,R,native,binary,i,n,j,full,expected,cap,actual,t,p,q,a,b,
        backend,source,memo,word;
    groups:=[CyclicGroup(2),CyclicGroup(3),CyclicGroup(4),AbelianGroup([2,2])];
    for G in groups do
        R:=ResolutionFiniteGroup(G,9);
        native:=koAHSSNativeCoherence(R); binary:=native.binaryTensor();
        backend:=koAHSSHAPSpace(R).koAHSS(0,0,6);
        for i in [0..3] do
            for n in [0..5] do
                if n+i>6 and Size(G)<>2 then continue; fi;
                if Size(G)=4 and not IsCyclic(G) and n+i>4 then continue; fi;
                for j in [1..R!.dimension(n)] do
                    full:=native.diagonal(i).evaluate(n,j);
                    expected:=List(Filtered(full,t->IsOddInt(t[1])),t->
                        Concatenation(t[2][1],t[2][2]));
                    Sort(expected);
                    for cap in [0..n+i] do
                        actual:=binary.diagonal(i,n,j,cap);
                        Assert(0,actual=Filtered(expected,t->t[1]<=cap and t[4]<=cap));
                    od;
                od;
            od;
        od;
        # Nonclosed and nonbinary integer vectors ensure agreement is literal
        # modulo two, rather than only on cohomology classes.
        for p in [0..3] do
            for q in [0..3] do
                if not IsCyclic(G) and p+q>4 then continue; fi;
                a:=List([1..R!.dimension(p)],j->(-1)^j*(j+2));
                b:=List([1..R!.dimension(q)],j->2-j);
                for i in [0..Minimum(p,q)] do
                    source:=native.evaluateCochains(native.diagonal(i),[p,q],[a,b],[0,0]);
                    Assert(0,backend.cupMod2(i,p,a,q,b)=List(source,x->x mod 2));
                od;
            od;
        od;
        Assert(0,binary.stats.contractionHits>0 and binary.stats.diagonalHits>0);
        memo:=KOAHSS_ResolutionMemo(R);
        for n in [0..2] do
            for j in [1..R!.dimension(n)] do
                for t in [1..Length(R!.elts)] do
                    word:=StructuralCopy(R!.homotopy(n,[j,t]));
                    Assert(0,memo.contract(n,j,t)=word);
                    Assert(0,memo.contract(n,j,t)=word);
                od;
            od;
        od;
    od;
    # Explicit lists retain their first occurrence after an external append.
    R!.elts:=ShallowCopy(R!.elts); t:=R!.elts[1]; Add(R!.elts,t);
    Assert(0,memo.index(t)=Position(R!.elts,t));
    # Infinite HAP lazy lists supply an inverse index; never enumerate them.
    R:=ResolutionAbelianGroup([0],5); binary:=KOAHSS_BinaryTensorEngine(R);
    G:=R!.group; t:=GeneratorsOfGroup(G)[1]^17;
    Assert(0,binary.index(t)=Position(R!.elts,t));
    Assert(0,binary.index(t^-1)=Position(R!.elts,t^-1));
end,[]);
