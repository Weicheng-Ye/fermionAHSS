# Synthetic linear algebra fixture: primary B rows followed by a tier-two
# A reference, both with uninitialized markings and with issued dependent rows.
# Formula tasks are stubbed; this checks frame bookkeeping, not identities.
LightFrameFixture:=function(prior,kind)
    local dim,zero,mkData,dataB,dataC,dataD,trivial,backend,layers,name,layer,
        oracle,lower,target,frame,rows,lowerB,recorded,before,i,answer;
    dim:=function(n) if n in [3,4] then return 2; fi; return 1; end;
    zero:=n->List([1..dim(n)],i->0);
    mkData:=function(orders,n)
        local group,gens;
        group:=AbelianGroup(IsPcpGroup,orders);
        gens:=IndependentGeneratorsOfAbelianGroup(group);
        return rec(group:=group,
            class:=function(v)
                local x,i; x:=One(group);
                for i in [1..Length(gens)] do x:=x*gens[i]^v[i]; od;
                return x;
            end,
            represent:=function(x)
                local v,e,i; v:=zero(n); e:=IndependentGeneratorExponents(group,x);
                for i in [1..Length(e)] do v[i]:=e[i]; od;
                return v;
            end);
    end;
    dataB:=mkData([2,2],3); dataC:=mkData([2,2],4);
    dataD:=mkData([2],6); trivial:=mkData([],1);
    backend:=rec(dimension:=dim,twists:=rec(omega:=[0]));
    backend.cohomologyData:=function(n,q)
        if n=3 then return dataB; fi;
        if n=4 then return dataC; fi;
        if n=6 then return dataD; fi;
        return trivial;
    end;
    backend.coboundary:=function(n,v,signed)
        if signed and n=1 then return [4*v[1]]; fi;
        if signed and n=3 then return [0,2*Sum(v)]; fi;
        return zero(n+1);
    end;
    backend.nativePrimary:=function(op,n,v)
        if op="D" and n=1 then
            if kind="pureC" then return [0,0]; fi;
            if kind="remark" then return [v[1] mod 2,0]; fi;
            return [v[1] mod 2,v[1] mod 2];
        fi;
        if op="Dtilde" then return zero(n+3); fi;
        return zero(n+2);
    end;
    layers:=rec(A:=rec(orders:=[4],cochains:=[[1]]),
        B:=rec(orders:=[2,2],cochains:=[[1,0],[0,1]],group:=dataB.group),
        C:=rec(orders:=[2,2],cochains:=[[1,0],[0,1]],group:=dataC.group),
        D:=rec(orders:=[2],cochains:=[[1]],group:=dataD.group));
    for name in ["A","B","C","D"] do
        layer:=layers.(name); layer.name:=name;
        if IsBound(layer.group) then
            layer.generators:=IndependentGeneratorsOfAbelianGroup(layer.group);
            layer.cell:=rec(project:=x->x);
        fi;
    od;
    oracle:=function(l,i,m,h)
        local v;
        if l.name="C" then v:=[0]; if i=1 then v:=[1]; fi;
        else v:=[0,0,1]; fi;
        return rec(status:="computed",lowerPresentationId:=h.presentationId,
            lowerCoordinates:=v,witness:=rec());
    end;
    Assert(0,KOAHSS_ExtensionTargetLayer(KOAHSS_ExtensionLowerPresentation(layers,["D","C"],oracle),2).layer="C");
    lower:=KOAHSS_ExtensionLowerPresentation(layers,["D","C","B"],oracle);
    target:=KOAHSS_ExtensionTargetLayer(lower,4);
    Assert(0,target.layer="D");
    rows:=function(n,signed)
        return List(IdentityMat(dim(n)),v->backend.coboundary(n,v,signed));
    end;
    frame:=KOAHSS_LightFrame(rec(backend:=backend,k:=5,layers:=layers,context:=rec(),
        part:=rec(prime:=2,generators:=rec(A:=[1],B:=[1,2],C:=[],D:=[])),coboundaryRows:=rows,
        call:=function(task,d)
            local out,w;
            out:=rec();
            if task="qd" then out.QD:=zero(5);
            elif task="atom_curvature" then out.J:=zero(7);
            elif task="c_mark" then out.J:=zero(7); out.gamma:=zero(6);
            elif task="gauge" then out.xC:=zero(4);
            elif task="bd_exact" then out.Z:=zero(6);
            elif task="bd_page" then out.Phi:=zero(6);
            elif task="a_step" then
                for w in d.want do
                    if w="residue" then out.(w):=zero(6);
                    elif w="QDa" then out.(w):=zero(4);
                    elif w="fsharp" then out.(w):=zero(5);
                    elif w="pot" then out.(w):=zero(7);
                    elif w="Ysrc" then out.(w):=zero(3);
                    elif w="RC" then
                        out.(w):=zero(4);
                        if kind="pureC" and d.ref.kind="zero" then out.(w):=[0,1]; fi;
                        if kind="remark" and d.ref.kind="atom" and ForAll(d.ref.c,IsZero) then
                            out.(w):=[1,0];
                        fi;
                    else Error("unexpected A task: ",w); fi;
                od;
            else Error("unexpected task: ",task); fi;
            return rec(status:="computed",result:=out);
        end));
    lowerB:=KOAHSS_ExtensionLowerPresentation(layers,["D","C"],oracle);
    before:=[];
    if prior then
        for i in [1,2] do
            Add(before,frame.row("B",i,2,lowerB,fail,rec()).witness.light.version);
        od;
    fi;
    answer:=frame.row("A",1,4,lower,target,rec());
    if kind="tier2" then
        Assert(0,answer.lowerCoordinates=[0,0,0,1,1]);
        Assert(0,"tier-2" in answer.witness.light.shortcuts);
    elif kind="pureC" then
        Assert(0,answer.lowerCoordinates=[0,0,1,0,0]);
        Assert(0,answer.witness.light.mode="exact");
    else
        Assert(0,answer.lowerCoordinates=[0,0,0,1,0]);
        Assert(0,"re-marking" in answer.witness.light.shortcuts);
        Assert(0,answer.witness.light.mode="strong");
    fi;
    Assert(0,answer.witness.light.remarked=[["B",1],["B",2]]);
    recorded:=rec(A:=[answer],B:=[],C:=[],D:=[]);
    for i in [1,2] do
        recorded.B[i]:=frame.row("B",i,2,lowerB,fail,rec(complete:=true,through:=3));
        if prior then Assert(0,recorded.B[i].witness.light.version>before[i]); fi;
    od;
    Assert(0,frame.audit(recorded,n->lower)=fail);
    recorded.B[1].witness.light.references[1][3]:=0;
    Assert(0,frame.audit(recorded,n->lower)<>fail);
end;
LightFrameFixture(false,"tier2");
LightFrameFixture(true,"tier2");
LightFrameFixture(true,"pureC");
LightFrameFixture(true,"remark");
Unbind(LightFrameFixture);
