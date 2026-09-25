# Reproduce the C2h higher-gauge relation without rerunning the point-group batch:
#   gap -q --quitonbreak examples/c2h_gauge.g
# Paths are relative to this file, so ReadPackage and absolute Read also work.
# The saved marked C/D representatives and lower coordinate come from the dated
# fixture. The differential, product, and boundary equality are recomputed.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    CallFuncList(function()
        local file, slash, prefix;
        file:=INPUT_FILENAME(); slash:=Positions(file,'/'); prefix:="";
        if not IsEmpty(slash) then prefix:=file{[1..Last(slash)]}; fi;
        Read(Concatenation(prefix,"../load.g"));
    end,[]);
fi;
CallFuncList(function()
    local file,slash,prefix,fixture,saved,G,R,backend,model,kernelCache,
        C,D,target,layers,lower,result,proof,check,oldBreak,caught;
    file:=INPUT_FILENAME(); slash:=Positions(file,'/'); prefix:="";
    if not IsEmpty(slash) then prefix:=file{[1..Last(slash)]}; fi;
    fixture:=JsonStringToGap(StringFile(Concatenation(
        prefix,"../data/extension-gauge-results-20260926.json")));
    saved:=fixture.savedCheck;
    Reset(GlobalMersenneTwister,1); Reset(GlobalRandomSource,1);
    G:=AbelianGroup([2,2]); R:=ResolutionFiniteGroup(G,7);
    backend:=koAHSSHAPSpace(R,koAHSSNaturalOperations()).koAHSS([1,0],[1,1,0],7);
    model:=KOAHSS_ExtensionBarModel(backend,4);
    Assert(0,model.status="computed");
    kernelCache:=rec();
    model.gaugeKernelRepresentatives:=function(n,signed,family)
        local key,q,data;
        key:=Concatenation(String(n),"_",String(signed));
        if not IsBound(kernelCache.(key)) then
            q:=-1; if signed then q:=0; fi;
            data:=backend.cohomologyData(n,q);
            kernelCache.(key):=List(IndependentGeneratorsOfAbelianGroup(data.group),
                g->model.lift(n,data.represent(g),signed));
        fi;
        return kernelCache.(key);
    end;
    check:=function()
        # Check the reconstructed complete-bar twist before reusing the
        # fixture cochains; the following equations also check their states.
        Assert(0,model.s=saved.s and model.omega=saved.omega);
        C:=saved.C; D:=saved.D;
        Assert(0,KOAHSS_ExtensionStateIsZero(model.d(4,C)));
        Assert(0,KOAHSS_ExtensionStateIsZero(model.d(4,D)));
        target:=model.xtimes(4,C,C);
        Assert(0,target=saved.target);
        Assert(0,koAHSSSolveIntegerSystem(
            Concatenation([D.D],model.matrix(4,true)),target.D)=fail);
        layers:=rec(D:=rec(orders:=[2],fullLifts:=[rec(state:=D)]));
        lower:=rec(generatorCount:=1,generatorIds:=["D:1"],
            presentationId:="c2h-regression:D",
            layers:=[rec(name:="D",orders:=[2],startColumn:=1)]);
        # Recheck the saved lower coordinate directly. Supplying it bounds
        # the search but does not certify the relation: only the full ordered
        # cochain equality below can establish that the coordinate is correct.
        Assert(0,saved.result.lowerCoordinates=[1]);
        result:=KOAHSS_ExtensionGaugeReduce(model,4,target,layers,lower,
            saved.result.lowerCoordinates);
        Assert(0,result.status="computed" and result.lowerCoordinates=[1]);
        proof:=result.canonicalComparison;
        Assert(0,proof.winningStage="CD" and proof.equalityVerified
            and proof.boundaryFlatnessVerified);
        Assert(0,result.canonicalLowerProduct=D);
        Assert(0,ForAll(proof.gauge.A,x->x=0) and ForAll(proof.gauge.B,x->x=0));
        Assert(0,model.d(3,proof.gauge)=proof.boundary);
        Assert(0,KOAHSS_ExtensionStateIsZero(model.d(4,proof.boundary)));
        Assert(0,model.xtimes(4,proof.boundary,D)=target);
    end;
    oldBreak:=BreakOnError; BreakOnError:=false;
    caught:=CALL_WITH_CATCH(check,[]);
    BreakOnError:=oldBreak; model.close();
    if not caught[1] then Error("C2h higher-gauge regression failed; see the preceding error"); fi;
    Print("C2h: 2C = D via a C/D gauge; exact cochain equality verified.\n");
end,[]);
