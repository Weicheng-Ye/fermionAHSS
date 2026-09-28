# Degrees one and two use the native transferred model, as degrees 3-5 do.
gap> stackUnitary := koFull_batch(CyclicGroup(2),0,0,2);;
gap> Assert(0,stackUnitary.status="computed" and stackUnitary.invariants=[[0],[],[2,2],[2,2]]);
gap> Assert(0,ForAll(stackUnitary.degreeResults{[3,4]},r->r.modelId="transferred-normalized-bar"
>     and r.certificateLevel="transfer-R" and r.gaugeCompletenessAssumed));
gap> Assert(0,ForAll(stackUnitary.degreeResults{[1,2]},r->not IsBound(r.modelId)));
gap> stackSignedAHSS := koAHSS_batch(CyclicGroup(2),[1],0,2,5,rec(details:=true));;
gap> stackSigned := koFull_batch(stackSignedAHSS);;
gap> Assert(0,IsIdenticalObj(stackSigned.ahss,stackSignedAHSS));
gap> Assert(0,IsIdenticalObj(stackSigned.ahss._context,stackSignedAHSS._context));
gap> Assert(0,stackSigned.invariants=[[],[2],[2],[8]]);
gap> Assert(0,IsIdenticalObj(stackSigned.pages,stackSignedAHSS.pages));
gap> Assert(0,koAHSSFormat(stackSigned)=koAHSSFormat(stackSignedAHSS.pages,6));
gap> stackEight := stackSigned.degreeResults[4];;
gap> Assert(0,HermiteNormalFormIntegerMat(stackEight.relationMatrix)
>     =HermiteNormalFormIntegerMat([[2,0,0],[-1,2,0],[0,-1,2]]));
gap> Assert(0,stackEight.smith.U*stackEight.relationMatrix*stackEight.smith.V=stackEight.smith.S);
gap> stackRelations := Filtered(stackEight.extensionVectors,v->v.layer in ["B","C"]);;
gap> Assert(0,ForAll(stackRelations,v->v.result.witness.power=2
>     and v.result.witness.modelId="transferred-normalized-bar"));
gap> # The omega twist changes the actual C-doubling cochain, producing Z/4.
gap> stackTwisted := koFull_batch(CyclicGroup(2),0,[1],1);;
gap> Assert(0,stackTwisted.invariants=[[0],[],[4]]);
gap> stackTwistedC := First(stackTwisted.degreeResults[3].extensionVectors,v->v.layer="C").result;;
gap> Assert(0,stackTwistedC.lowerCoordinates=[1] and stackTwistedC.witness.power=2);
gap> # Two independent C generators retain separate lower D coordinates;
gap> # B doubles into their sum.
gap> stackTwo := koFull_batch(AbelianGroup([2,2]),[1,1],0,2);;
gap> Assert(0,stackTwo.invariants=[[],[2],[2,2],[4,8]]);
gap> stackTwoC := Filtered(stackTwo.degreeResults[4].extensionVectors,v->v.layer="C");;
gap> Assert(0,List(stackTwoC,v->v.result.lowerCoordinates)=[[1,0],[0,1]]);
gap> stackTwoB := First(stackTwo.degreeResults[4].extensionVectors,v->v.layer="B").result;;
gap> Assert(0,stackTwoB.lowerCoordinates{[3,4]}=[1,1]
>     and ForAll(stackTwoB.lowerCoordinates{[1,2]},IsEvenInt));
gap> # Nonzero omega in degree two, and the cell comparison in degrees one and two.
gap> stackOmega := koFull_batch(CyclicGroup(2),[1],[1],2);;
gap> Assert(0,stackOmega.invariants=[[],[2],[],[2]]);
gap> KOAHSS_EXTENSION_TRANSPORT_OVERRIDE.cells:=true;;
gap> stackCellModel := KOAHSS_ExtensionTransferredModel(stackSignedAHSS._context.backend,1);;
gap> Assert(0,stackCellModel.status="computed" and stackCellModel.comparison="cells" and stackCellModel.supports(1));
gap> stackCells := koFull_batch(stackSignedAHSS);;
gap> KOAHSS_EXTENSION_TRANSPORT_OVERRIDE.cells:=false;;
gap> Assert(0,stackCells.invariants=stackSigned.invariants);
gap> Assert(0,List(stackCells.degreeResults,r->r.relationMatrix)=List(stackSigned.degreeResults,r->r.relationMatrix));
