# Compare maps and marked page representatives, not only invariant factors.
# Direct calibrated cochain APIs remain the independent bar reference.
CallFuncList(function()
    local compare,configs,config,R,w,ops,native,B,H,K,n,q,name,g,a,
        gens,N,O,z,dz,P,h,shifted;
    compare:=function(R,s,w,k)
        local ops,native,bar,key,cell,other,map,otherMap,g;
        ops:=koAHSSNaturalOperations(); ops.useNativePages:=true;
        native:=KOAHSS_ComputePageContext([koAHSSHAPSpace(R,ops),s,w,k,5]);
        ops:=koAHSSNaturalOperations(); ops.useNativePages:=false;
        bar:=KOAHSS_ComputePageContext([koAHSSHAPSpace(R,ops),s,w,k,5]);
        Assert(0,RecNames(native.cells)=RecNames(bar.cells));
        Assert(0,RecNames(native.maps)=RecNames(bar.maps));
        for key in RecNames(native.cells) do
            cell:=native.cells.(key); other:=bar.cells.(key);
            Assert(0,not KOAHSS_IsUnresolved(cell) and not KOAHSS_IsUnresolved(other));
            Assert(0,AbelianInvariants(cell.group)=AbelianInvariants(other.group));
            Assert(0,List(GeneratorsOfGroup(cell.group),g->Exponents(cell.lift(g)))=
                List(GeneratorsOfGroup(other.group),g->Exponents(other.lift(g))));
            # The actual integral/binary representatives in R also agree.
            Assert(0,List(GeneratorsOfGroup(cell.group),g->
                native.backend.data(cell.degree,cell.q).represent(cell.lift(g)))=
                List(GeneratorsOfGroup(other.group),g->
                bar.backend.data(other.degree,other.q).represent(other.lift(g))));
        od;
        for key in RecNames(native.maps) do
            map:=native.maps.(key); otherMap:=bar.maps.(key);
            if IsRecord(map) then Assert(0,map=otherMap);
            else
                Assert(0,IsGeneralMapping(otherMap));
                Assert(0,List(GeneratorsOfGroup(Source(map)),g->Exponents(Image(map,g)))=
                    List(GeneratorsOfGroup(Source(otherMap)),g->Exponents(Image(otherMap,g))));
            fi;
        od;
        return native;
    end;
    configs:=[[[2],0,0,4],[[2],[1],0,4],[[2],0,[1],4],[[2],[1],[1],4],
        [[3],0,0,4],[[4],0,0,4],[[4],[1],0,6],[[4],0,[1],4],
        [[2,2],0,0,3]];
    for config in configs do
        R:=ResolutionAbelianGroup(config[1],config[4]+3);
        native:=compare(R,config[2],config[3],config[4]);
        B:=native.backend;
        # Primary operations agree in cohomology, including nonzero arrows,
        # without requiring their different cochain representatives to agree.
        for n in [0..Minimum(3,config[4]-1)] do
            for name in ["Dbar","D","Dtilde"] do
                q:=-1; if name="Dbar" then q:=0; fi;
                H:=B.data(n,q); K:=B.data(n+2,-2);
                if name="Dtilde" then K:=B.data(n+3,-4); fi;
                for g in GeneratorsOfGroup(H.group) do
                    a:=H.represent(g);
                    Assert(0,K.class(B.nativePrimary(name,n,a))=K.class(B.primary(name,n,a)));
                od;
            od;
        od;
    od;
    # A nonzero T0 class: integral H^5(C4 x C4) has order-four coordinates.
    R:=ResolutionAbelianGroup([4,4],7);
    ops:=koAHSSNaturalOperations(); ops.useNativePages:=true;
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,0,6);
    H:=B.data(1,-1); gens:=GeneratorsOfGroup(H.group);
    w:=B.cupMod2(0,1,H.represent(gens[1]),1,H.represent(gens[2]));
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,w,6);
    N:=KOAHSS_NativeT0Page(B,0,[2]); O:=koAHSSNaturalTertiary(B,0,[2]);
    Assert(0,N.status="computed" and O.status="computed");
    K:=B.data(5,-4);
    Assert(0,not IsOne(K.class(N.cochain)));
    Assert(0,K.class(N.cochain)=K.class(O.cochain));
    Assert(0,N.inputGroupBarUsed=false and N.pageClassOnly=true and O.inputGroupBarUsed=true);
    # All other Z/4 lifts change the result by the advertised Dtilde class.
    H:=B.data(2,-1);
    for g in GeneratorsOfGroup(H.group) do
        h:=H.represent(g); z:=N.definingSystem.omegaLiftModFour+2*h;
        dz:=B.coboundary(2,z,false);
        P:=B.cupIntegral(0,2,z,2,z,false,false)+B.cupIntegral(1,2,z,3,dz,false,false);
        shifted:=-B.coboundary(4,P,false)/8;
        Assert(0,ForAll(shifted,IsInt));
        Assert(0,K.class(shifted-N.cochain)=K.class(B.nativePrimary("Dtilde",2,h)));
    od;
    for n in [4,6,8] do
        O:=koAHSSNaturalTertiary(B,0,[n]); shifted:=KOAHSS_NativeT0Page(B,0,[n]);
        Assert(0,K.class(shifted.cochain)=K.class(O.cochain));
    od;
    # The native shortcut must not even construct a bar comparison.
    B.naturalTransport:=function() Error("unexpected bar transport"); end;
    B.naturalBar:=B.naturalTransport;
    B.primary:=function(arg) Error("page primary used the defining-system representative"); end;
    KOAHSS_ComputePageContext([rec(koAHSS:=function(arg) return B; end),0,w,4,2]);
    Assert(0,KOAHSS_NativeT0Page(B,0,[2]).cochain=N.cochain);
    ops:=koAHSSNaturalOperations(); ops.useNativePages:=true;
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,w,6);
    B.naturalTransport:=function() Error("unexpected bar transport"); end;
    B.naturalBar:=B.naturalTransport;
    Assert(0,ops.T(rec(backend:=B,degree:=0,cochain:=[2]))=N.cochain);
    a:=2*B.data(2,0).represent(GeneratorsOfGroup(B.data(2,0).group)[1]);
    Assert(0,ops.Tau(rec(backend:=B,degree:=2,cochain:=a))=
        List(B.nativePrimary("Dtilde",2,List(a,x->(x/2) mod 2)),x->x mod 2));
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,0,6);
    B.naturalTransport:=function() Error("unexpected bar transport"); end;
    B.naturalBar:=B.naturalTransport;
    Assert(0,ops.Psi(rec(backend:=B,degree:=0,cochain:=[1]))=
        List([1..B.dimension(4)],j->0));
    # Refuse the native T0 branch when omega has no closed mod-four lift.
    R:=ResolutionAbelianGroup([2,2],7);
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,0,6);
    H:=B.data(1,-1); gens:=GeneratorsOfGroup(H.group);
    w:=B.cupMod2(0,1,H.represent(gens[1]),1,H.represent(gens[2]));
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,w,6);
    Assert(0,KOAHSS_NativeT0Page(B,0,[2])=fail);
    Assert(0,KOAHSS_IsUnresolved(ops.Psi(rec(backend:=B,degree:=0,cochain:=[1]))));
    Assert(0,KOAHSS_NativeT0Page(B,1,List([1..B.dimension(1)],j->2))=fail);
    N:=KOAHSS_NativeT0Page(B,0,[4]); O:=koAHSSNaturalTertiary(B,0,[4]);
    K:=B.data(5,-4); Assert(0,K.class(N.cochain)=K.class(O.cochain));
    # A nonzero b can rescue the input 2: the full callback must fall back.
    ops.T(rec(backend:=B,degree:=0,cochain:=[2]));
    Assert(0,Last(B.tertiaryEvaluations).inputGroupBarUsed=true);
    native:=compare(R,0,w,3);
    # Odd inputs with an exact, nonzero omega use the R0o cancellation.
    R:=ResolutionAbelianGroup([3],7);
    B:=koAHSSHAPSpace(R,ops).koAHSS(0,[1],6);
    for n in [1,3,5] do
        N:=KOAHSS_NativeT0Page(B,0,[n]); O:=koAHSSNaturalTertiary(B,0,[n]);
        Assert(0,N.nativeOperation="degree-zero-odd" and ForAll(N.cochain,x->x=0));
        Assert(0,IsOne(B.data(5,-4).class(O.cochain)));
    od;
    native:=compare(R,0,[1],3);
end,[]);
