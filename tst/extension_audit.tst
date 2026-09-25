# Refuse a Smith prediction when the actual product has no boundary witness.
gap> auditModel := rec();;
gap> auditModel.dimension := function(n) if n<0 then return 0; fi; return 1; end;;
gap> auditModel.zero := k->rec(A:=List([1..auditModel.dimension(k-3)],i->0),B:=List([1..auditModel.dimension(k-2)],i->0),C:=List([1..auditModel.dimension(k-1)],i->0),D:=[0]);;
gap> auditModel.d := function(k,state) return auditModel.zero(k+1); end;;
gap> auditModel.coboundary := function(n,v,signed) return List([1..auditModel.dimension(n+1)],i->0); end;;
gap> auditModel.xtimes := function(k,x,y) local z; z:=auditModel.zero(k); z.D:=x.D+y.D; return z; end;;
gap> auditLift := rec(status:="computed",state:=rec(A:=[0],B:=[0],C:=[0],D:=[1]));;
gap> auditLayers := rec(A:=rec(orders:=[],fullLifts:=[]),B:=rec(orders:=[],fullLifts:=[]),C:=rec(orders:=[],fullLifts:=[]),D:=rec(orders:=[2],fullLifts:=[auditLift]));;
gap> auditCandidate := koAHSSExtensionFromLayers(rec(D:=[2]),fail);;
gap> auditFailed := KOAHSS_ExtensionFiniteAudit(auditModel,3,auditLayers,auditCandidate);;
gap> Assert(0,auditFailed.status="unresolved" and auditFailed.left=2 and auditFailed.right=2);
gap> # A literal cyclic product does pass every actual-pair comparison.
gap> auditModel.xtimes := function(k,x,y) local z; z:=auditModel.zero(k); z.D:=List(x.D+y.D,v->v mod 2); return z; end;;
gap> auditGood := KOAHSS_ExtensionFiniteAudit(auditModel,3,auditLayers,auditCandidate);;
gap> Assert(0,auditGood.status="computed" and auditGood.table=[[1,2],[2,1]] and auditGood.associativityVerified);
gap> auditLayers.D:=rec(orders:=List([1..6],i->2),fullLifts:=List([1..6],i->auditLift));;
gap> Assert(0,KOAHSS_ExtensionFiniteAudit(auditModel,3,auditLayers,auditCandidate).normalFormCount=64);
gap> auditLayers.D:=rec(orders:=[],fullLifts:=[]);;
gap> auditEmpty := KOAHSS_ExtensionFiniteAudit(auditModel,3,auditLayers,koAHSSExtensionFromLayers(rec(),fail));;
gap> Assert(0,auditEmpty.status="computed" and auditEmpty.table=[[1]]);
gap> # Real higher-degree lifts and all 64 torsion normal-form products.
gap> auditReal := koFull(CyclicGroup(2),0,0,3);;
gap> Assert(0,auditReal.invariants[5]=[0,8]);
gap> Assert(0,Length(auditReal.degreeResults[5].algebraAudit.table)=8);
gap> Assert(0,ForAll(auditReal.degreeResults[5].extensionVectors,v->v.layer="D" or v.result.witness.reduction.canonicalComparison.status="computed"));
