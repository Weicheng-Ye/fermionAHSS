# A relation whose target layer lies right below its generator's layer is read
# from a primary operation on the resolution (doc/extension_cup_i_formulas.md):
# D(rho u) with delta_s u = m A for A over B, the reduced twisted Bockstein for
# B over C, and an integral lift of D(c) for C over D. In its target layer the
# row agrees with the measurement in the transferred model, exactly in the
# binary layers and modulo two in D; every other row is the same, and so is
# the group. The other rows are measured in the model in both runs (light
# rows off; tst/extension_light_rows.tst compares those).
gap> PrimaryCompare := function(group,s,omega,k)
>     local ahss,on,off,rows,primary,v,w,j,id,ok;
>     ahss:=koAHSS_batch(group,s,omega,k,rec(details:=true));
>     KOAHSS_EXTENSION_RELATION_OVERRIDE.light:=false;
>     KOAHSS_EXTENSION_RELATION_OVERRIDE.native:=true;
>     on:=koFull(ahss);
>     KOAHSS_EXTENSION_RELATION_OVERRIDE.native:=false;
>     off:=koFull(ahss);
>     Unbind(KOAHSS_EXTENSION_RELATION_OVERRIDE.native);
>     Unbind(KOAHSS_EXTENSION_RELATION_OVERRIDE.light);
>     ok:=on.status="computed" and off.status="computed" and on.invariants=off.invariants
>         and not IsBound(on.degreeResult.primaryFallbacks);
>     for v in on.degreeResult.extensionVectors do
>         w:=First(off.degreeResult.extensionVectors,x->x.generatorId=v.generatorId);
>         if w=fail or w.lowerGeneratorIds<>v.lowerGeneratorIds then ok:=false; continue; fi;
>         primary:=IsBound(v.result.witness.model) and v.result.witness.model="primary-R";
>         for j in [1..Length(v.lowerGeneratorIds)] do
>             id:=v.lowerGeneratorIds[j];
>             if not primary then
>                 if v.result.lowerCoordinates[j]<>w.result.lowerCoordinates[j] then ok:=false; fi;
>             elif not [id[1]] in v.result.witness.measuredLayers then
>                 if v.result.lowerCoordinates[j]<>0 then ok:=false; fi;
>             elif id[1]='D' then
>                 if (v.result.lowerCoordinates[j]-w.result.lowerCoordinates[j]) mod 2<>0 then ok:=false; fi;
>             elif v.result.lowerCoordinates[j]<>w.result.lowerCoordinates[j] then ok:=false;
>             fi;
>         od;
>     od;
>     rows:=Filtered(on.degreeResult.extensionVectors,v->IsBound(v.result.witness.model)
>         and v.result.witness.model="primary-R");
>     return rec(ok:=ok,on:=on,off:=off,primary:=List(rows,v->v.layer));
> end;;
gap> # C2 with every twist in degrees one to three: C over D and B over C.
gap> primaryLow := List(Cartesian([0,[1]],[0,[1]],[1,2,3]),c->PrimaryCompare(CyclicGroup(2),c[1],c[2],c[3]));;
gap> Assert(0,ForAll(primaryLow,r->r.ok));
gap> Assert(0,ForAny(primaryLow,r->"C" in r.primary) and ForAny(primaryLow,r->"B" in r.primary));
gap> # Table VII, (0,[1]) in degree one: 2C=D from D(1)=omega, so Z/4.
gap> primaryTwisted := PrimaryCompare(CyclicGroup(2),0,[1],1);;
gap> Assert(0,primaryTwisted.ok and primaryTwisted.on.invariants=[4] and primaryTwisted.primary=["C"]);
gap> # Pin-: C2 with s at degree six. The lower group of the A generator is
gap> # <b,c,d | 2b=c, 2c=d, 2d=0>, every relation has an adjacent target, and
gap> # Z/16 is determined without the stacking model.
gap> primaryPin := PrimaryCompare(CyclicGroup(2),[1],0,6);;
gap> Assert(0,primaryPin.ok and primaryPin.on.invariants=[16]);
gap> Assert(0,SortedList(primaryPin.primary)=["A","B","C"]);
gap> Assert(0,primaryPin.on.degreeResult.certificateLevel="primary-R"
>     and not IsBound(primaryPin.on.degreeResult.modelId));
gap> Assert(0,primaryPin.off.degreeResult.certificateLevel="transfer-R");
gap> # Z4^Tf and Cs spinless: C2 with s=omega=[1] in degree four, 2A=B.
gap> primaryZ4T := PrimaryCompare(CyclicGroup(2),[1],[1],4);;
gap> Assert(0,primaryZ4T.ok and "A" in primaryZ4T.primary);
gap> # Two C generators with separate D coordinates, and C4: 2C=D from
gap> # D(y)=y^2, the reduction of a generator of H^4(Z/4;Z), in degree three.
gap> Assert(0,PrimaryCompare(AbelianGroup([2,2]),[1,1],0,2).ok);
gap> primaryC4 := PrimaryCompare(CyclicGroup(4),0,0,3);;
gap> Assert(0,primaryC4.ok and "C" in primaryC4.primary);
gap> Assert(0,PrimaryCompare(CyclicGroup(4),[1],0,2).ok);
