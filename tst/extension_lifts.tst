# A triangular model forces both lower choices to change: the first B
# primitive does not admit C, and the first C primitive does not admit D.
# The oracle is synthetic so failures isolate the defining-system solver.
gap> flatModel := rec(modelId:="affine-defining-system-test");;
gap> flatModel.dimension := function(n) if n<0 then return 0; fi; return 1; end;;
gap> flatModel.zero := k->rec(A:=List([1..flatModel.dimension(k-3)],i->0),B:=List([1..flatModel.dimension(k-2)],i->0),C:=List([1..flatModel.dimension(k-1)],i->0),D:=List([1..flatModel.dimension(k+1)],i->0));;
gap> flatModel.coboundary := function(n,v,signed) if n=4 and signed then return 2*v; fi; return List([1..flatModel.dimension(n+1)],i->0); end;;
gap> flatModel.d := function(k,x) Assert(0,k=3); return rec(A:=[0],B:=[0],C:=[x.A[1]*(1-x.B[1]) mod 2],D:=[x.C[1]-x.A[1]+2*x.D[1]]); end;;
gap> flatLift := KOAHSS_ExtensionFlatLift(flatModel,3,"A",[1]);;
gap> Assert(0,flatLift.status="computed" and flatLift.state=rec(A:=[1],B:=[1],C:=[1],D:=[0]));
gap> Assert(0,flatLift.witness.flatnessVerified and flatLift.witness.curvature=rec(A:=[0],B:=[0],C:=[0],D:=[0]));
gap> Assert(0,List(flatLift.witness.definingEquations,e->e.chosenPrimitive)=[[1],[1],[0]]);
gap> Assert(0,List(flatLift.witness.definingEquations,e->e.adjustment)=[[1],[1],[0]]);
gap> Assert(0,not IsMutable(flatLift.state));
gap> # Exhausting an explicit bound is not an obstruction or a split result.
gap> flatLimited := KOAHSS_ExtensionFlatLift(flatModel,3,"A",[1],rec(maxChoices:=2));;
gap> Assert(0,flatLimited.status="unresolved" and not IsBound(flatLimited.state));
gap> # A verified unchanged tuple consumes only its first differential audit.
gap> flatZero := KOAHSS_ExtensionFlatLift(flatModel,3,"A",[0],rec(maxChoices:=1));;
gap> Assert(0,flatZero.status="computed" and flatZero.witness.differentialEvaluations=1);
gap> # Integral D is solved over Z, retaining its nonzero signed primitive.
gap> flatIntegralModel := ShallowCopy(flatModel);;
gap> flatIntegralModel.d := function(k,x) return rec(A:=[0],B:=[0],C:=[0],D:=2*x.C+2*x.D); end;;
gap> flatIntegral := KOAHSS_ExtensionFlatLift(flatIntegralModel,3,"C",[1]);;
gap> Assert(0,flatIntegral.state.D=[-1] and flatIntegral.witness.definingEquations[1].modulus=0);
gap> # An actual obstruction remains distinct from a resource-limited search.
gap> flatObstructedModel := ShallowCopy(flatModel);;
gap> flatObstructedModel.d := function(k,x) return rec(A:=[0],B:=[0],C:=[0],D:=x.C+2*x.D); end;;
gap> flatObstructed := KOAHSS_ExtensionFlatLift(flatObstructedModel,3,"C",[1]);;
gap> Assert(0,flatObstructed.status="obstructed" and flatObstructed.lastFailure.layer="D");
gap> # The chosen full tuple is shared, rather than reconstructed with new B,C.
gap> flatModel.lift := function(n,v,signed) return ShallowCopy(v); end;;
gap> flatLayers := [rec(name:="A",p:=0,cochains:=[[1]])];;
gap> flatGet := KOAHSS_ExtensionLiftSolver(rec(),3,flatLayers,flatModel);;
gap> flatFirst := flatGet(flatLayers[1],1);;
gap> Assert(0,IsIdenticalObj(flatFirst,flatGet(flatLayers[1],1)));
gap> Assert(0,IsIdenticalObj(flatFirst,flatLayers[1].fullLifts[1]));
gap> # Production check: lift every marked generator, including the p+ip A
gap> # generator, in the same complete bar basis used by the actual xtimes.
gap> flatActualChecks := function(k,s,omega)
> local ahss,layers,model,get,name,layer,i,lift,curvature,count;
> ahss:=koAHSS(CyclicGroup(2),s,omega,k,rec(details:=true));
> layers:=KOAHSS_ExtensionLayers(ahss._context,k);
> model:=KOAHSS_ExtensionBarModel(ahss._context.backend,k);
> get:=KOAHSS_ExtensionLiftSolver(ahss._context.backend,k,layers,model);
> count:=0;
> for name in ["D","C","B","A"] do
>   layer:=layers.(name);
>   for i in [1..Length(layer.generators)] do
>     lift:=get(layer,i);
>     Assert(0,lift.status="computed" and lift.witness.flatnessVerified);
>     curvature:=model.d(k,lift.state);
>     Assert(0,ForAll(["A","B","C","D"],f->ForAll(curvature.(f),x->x=0)));
>     Assert(0,IsIdenticalObj(lift,get(layer,i)));
>     Assert(0,IsIdenticalObj(lift,layer.fullLifts[i]));
>     Assert(0,lift.state.(name)=model.lift(layer.p,layer.cochains[i],name in ["A","D"]));
>     count:=count+1;
>   od;
> od;
> model.close();
> return count;
> end;;
gap> Assert(0,flatActualChecks(3,0,0)=4);
gap> Assert(0,flatActualChecks(4,[1],[1])=4);
