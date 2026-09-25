# Production xtimes restrictions, exact relations and retained context reuse.
gap> stackUnitary := koFull(CyclicGroup(2),0,0,2);;
gap> Assert(0,stackUnitary.status="computed" and stackUnitary.invariants=[[0],[],[2,2],[2,2]]);
gap> stackSignedAHSS := koAHSS(CyclicGroup(2),[1],0,2,5,rec(details:=true));;
gap> stackSigned := koFull(stackSignedAHSS);;
gap> Assert(0,IsIdenticalObj(stackSigned.ahss,stackSignedAHSS));
gap> Assert(0,IsIdenticalObj(stackSigned.ahss._context,stackSignedAHSS._context));
gap> Assert(0,stackSigned.invariants=[[],[2],[2],[8]]);
gap> Assert(0,IsIdenticalObj(stackSigned.pages,stackSignedAHSS.pages));
gap> Assert(0,koAHSSFormat(stackSigned)=koAHSSFormat(stackSignedAHSS.pages,6));
gap> stackEight := stackSigned.degreeResults[4];;
gap> Assert(0,stackEight.relationMatrix=[[2,0,0],[-1,2,0],[0,-1,2]]);
gap> Assert(0,stackEight.smith.U*stackEight.relationMatrix*stackEight.smith.V=stackEight.smith.S);
gap> stackRelations := Filtered(stackEight.extensionVectors,v->v.layer in ["B","C"]);;
gap> Assert(0,List(stackRelations,v->v.result.lowerCoordinates)=[[1],[0,1]]);
gap> Assert(0,ForAll(stackRelations,v->v.result.witness.xtimesCalls=1));
gap> Assert(0,ForAll(stackRelations,v->v.result.witness.calibration=KOAHSS_STACKING_LOW_PHASE.calibration));
gap> # The omega twist changes the actual C-doubling cochain, producing Z/4.
gap> stackTwisted := koFull(CyclicGroup(2),0,[1],1);;
gap> Assert(0,stackTwisted.invariants=[[0],[],[4]]);
gap> stackTwistedWitness := First(stackTwisted.degreeResults[3].extensionVectors,v->v.layer="C").result.witness;;
gap> Assert(0,stackTwistedWitness.lowerD=[1] and stackTwistedWitness.upperPrimitive=[0]);
gap> Assert(0,stackTwistedWitness.upperObstruction=[0] and stackTwistedWitness.xtimesCalls=1);
gap> # Two independent C generators retain separate lower D coordinates.
gap> # B doubles into their sum; reduce by stacking both chosen C lifts.
gap> stackTwo := koFull(AbelianGroup([2,2]),[1,1],0,2);;
gap> Assert(0,stackTwo.invariants=[[],[2],[2,2],[4,8]]);
gap> stackTwoC := Filtered(stackTwo.degreeResults[4].extensionVectors,v->v.layer="C");;
gap> Assert(0,List(stackTwoC,v->v.result.lowerCoordinates)=[[1,0],[0,1]]);
gap> stackTwoB := First(stackTwo.degreeResults[4].extensionVectors,v->v.layer="B").result;;
gap> Assert(0,stackTwoB.lowerCoordinates=[0,0,1,1]);
gap> Assert(0,stackTwoB.witness.lowerCCoordinates=[1,1] and stackTwoB.witness.lowerReductionProducts=2);
gap> # Unsupported omega transport remains an explicit missing relation.
gap> stackUnsupportedAHSS := koAHSS(CyclicGroup(2),0,[1],2,rec(details:=true));;
gap> stackUnsupportedLayers := KOAHSS_ExtensionLayers(stackUnsupportedAHSS._context,2);;
gap> stackUnsupportedOracle := KOAHSS_StackingExtensionOracle(stackUnsupportedAHSS._context.backend,2,stackUnsupportedLayers);;
gap> stackUnsupported := stackUnsupportedOracle(stackUnsupportedLayers.C,1,2,rec());;
gap> Assert(0,stackUnsupported.status="unresolved" and not IsBound(stackUnsupported.lowerCoordinates));
gap> Assert(0,Length(KOAHSS_STACKING_LOW_PHASE.values)=256 and not IsMutable(KOAHSS_STACKING_LOW_PHASE));
