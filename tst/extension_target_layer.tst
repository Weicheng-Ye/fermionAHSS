# The target layer of a relation: the lowest layer of the recorded lower
# presentation with a generator outside n*G_low. Generators below it cannot
# change the extension class, and the certificate exhibits their multiples.
gap> # Pin-: lower group <b,c,d | 2b=c, 2c=d, 2d=0> = Z/8 (layers D, C, B added in this order).
gap> targetLowerPin := rec(generatorCount:=3,generatorIds:=["D:1","C:1","B:1"],relationMatrix:=[[2,0,0],[-1,2,0],[0,-1,2]],layers:=[rec(name:="D",startColumn:=1,orders:=[2]),rec(name:="C",startColumn:=2,orders:=[2]),rec(name:="B",startColumn:=3,orders:=[2])]);;
gap> targetPin := KOAHSS_ExtensionTargetLayer(targetLowerPin,2);;
gap> Assert(0,targetPin.layer="B" and targetPin.index=1 and targetPin.freeGenerators=["B:1"]);
gap> Assert(0,List(targetPin.certificate,c->c.generatorId)=["D:1","C:1"]);
gap> Assert(0,ForAll(targetPin.certificate,c->2*c.combination+c.relationCoefficients*targetLowerPin.relationMatrix=List([1..3],function(i) if i=c.column then return 1; fi; return 0; end)));
gap> # Order four on the same lower group: 4*G_low = <d>, so c is free as well and the target is C.
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerPin,4).layer="C");
gap> # Z/2 + Z/4 with 2b=0, 2c=d: b and c are free, d=2c is not; the target is C, not B.
gap> targetLowerMixed := rec(generatorCount:=3,generatorIds:=["D:1","C:1","B:1"],relationMatrix:=[[2,0,0],[-1,2,0],[0,0,2]],layers:=targetLowerPin.layers);;
gap> targetMixed := KOAHSS_ExtensionTargetLayer(targetLowerMixed,2);;
gap> Assert(0,targetMixed.layer="C" and targetMixed.index=2 and targetMixed.freeGenerators=["C:1","B:1"]);
gap> Assert(0,List(targetMixed.certificate,c->c.generatorId)=["D:1"]);
gap> # A split D layer (Z/2)^3: every generator is free and the complete measurement is needed.
gap> targetLowerSplit := rec(generatorCount:=3,generatorIds:=["D:1","D:2","D:3"],relationMatrix:=[[2,0,0],[0,2,0],[0,0,2]],layers:=[rec(name:="D",startColumn:=1,orders:=[2,2,2])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerSplit,2).layer="D");
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerSplit,2).index=3);
gap> # An odd lower group and an even order: nothing is free, the class vanishes.
gap> targetLowerOdd := rec(generatorCount:=1,generatorIds:=["D:1"],relationMatrix:=[[3]],layers:=[rec(name:="D",startColumn:=1,orders:=[3])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerOdd,2).layer=fail);
gap> Assert(0,Length(KOAHSS_ExtensionTargetLayer(targetLowerOdd,2).certificate)=1);
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerOdd,3).layer="D");
gap> # A free lower generator of infinite order is free for every order.
gap> targetLowerFree := rec(generatorCount:=1,generatorIds:=["B:1"],relationMatrix:=[],layers:=[rec(name:="B",startColumn:=1,orders:=[0])]);;
gap> Assert(0,KOAHSS_ExtensionTargetLayer(targetLowerFree,2).layer="B");
gap> # No lower generators: no target.
gap> Assert(0,KOAHSS_ExtensionTargetLayer(rec(generatorCount:=0,generatorIds:=[],relationMatrix:=[],layers:=[]),2).layer=fail);
gap> Assert(0,KOAHSS_ExtensionStateIsZeroThrough(rec(A:=[0],B:=[0],C:=[1],D:=[3]),1));
gap> Assert(0,not KOAHSS_ExtensionStateIsZeroThrough(rec(A:=[0],B:=[0],C:=[1],D:=[3]),2));
gap> Assert(0,KOAHSS_ExtensionStatesAgreeThrough(rec(A:=[1],B:=[1],C:=[0],D:=[5]),rec(A:=[1],B:=[1],C:=[1],D:=[0]),1));
gap> Assert(0,IsBool(KOAHSS_LayeredRelationsEnabled()));
