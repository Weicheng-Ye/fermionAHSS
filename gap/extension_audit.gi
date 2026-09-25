# Finite quotient audit using complete flat normal forms and actual boundaries.
BindGlobal("KOAHSS_ExtensionFiniteAudit",function(model,k,layers,result)
    local generators,ranges,free,layer,name,j,choices,forms,keys,coordinates,
        digits,state,power,smithCoordinate,key,table,proofs,i,target,product,proof,
        order,identity,inverses,a,b,c;
    generators:=[]; ranges:=[]; free:=[];
    for name in ["D","C","B","A"] do
        layer:=layers.(name);
        for j in [1..Length(layer.orders)] do
            Add(generators,layer.fullLifts[j].state);
            if layer.orders[j]=0 then Add(ranges,[0]); Add(free,Length(generators));
            else Add(ranges,[0..layer.orders[j]-1]); fi;
        od;
    od;
    order:=Product(List(ranges,Length));
    if order>32 then return rec(status:="unresolved",
        reason:="the explicit finite stacking quotient audit exceeds 32 normal forms",normalFormCount:=order); fi;
    smithCoordinate:=function(row)
        local value,i;
        if IsEmpty(row) then return []; fi;
        value:=(row*result.smith.V){result.smith.activeIndices};
        for i in [1..Length(value)] do
            if result.basis.orders[i]<>0 then value[i]:=value[i] mod result.basis.orders[i]; fi;
        od;
        return value;
    end;
    choices:=Cartesian(ranges); forms:=[]; keys:=[];
    for digits in choices do
        state:=model.zero(k);
        for j in [1..Length(generators)] do
            power:=KOAHSS_ExtensionPower(model,k,generators[j],digits[j]);
            if not KOAHSS_ExtensionStateIsZero(power) then
                if KOAHSS_ExtensionStateIsZero(state) then state:=power;
                else state:=model.xtimes(k,state,power); fi;
            fi;
        od;
        Add(forms,state); Add(keys,smithCoordinate(digits));
    od;
    if Length(Set(keys))<>order then
        Error("koFull: measured Smith relations contradict the marked filtration normal forms");
    fi;
    table:=[]; proofs:=[];
    for i in [1..order] do
        table[i]:=[]; proofs[i]:=[];
        for j in [1..order] do
            key:=smithCoordinate(choices[i]+choices[j]); target:=Position(keys,key);
            if target=fail then return rec(status:="unresolved",
                reason:="the finite digit representatives do not form a closed torsion subgroup"); fi;
            product:=model.xtimes(k,forms[i],forms[j]);
            if product=forms[target] then
                proof:=rec(status:="computed",kind:="literal-equality",equalityVerified:=true);
            else
                proof:=KOAHSS_ExtensionGaugeCompare(model,k,product,forms[target]);
            fi;
            if proof.status<>"computed" then return rec(status:="unresolved",
                reason:="a finite normal-form product lacks an exact boundary comparison",
                left:=i,right:=j,target:=target,comparison:=proof); fi;
            table[i][j]:=target; proofs[i][j]:=proof;
        od;
    od;
    identity:=Position(choices,List(generators,g->0)); inverses:=[];
    for a in [1..order] do
        if table[a][identity]<>a or table[identity][a]<>a then
            Error("koFull: audited multiplication lacks the claimed identity"); fi;
        inverses[a]:=First([1..order],b->table[a][b]=identity and table[b][a]=identity);
        if inverses[a]=fail then Error("koFull: audited multiplication lacks an inverse"); fi;
        for b in [1..order] do
            if table[a][b]<>table[b][a] then Error("koFull: audited multiplication is not commutative"); fi;
            for c in [1..order] do
                if table[table[a][b]][c]<>table[a][table[b][c]] then
                    Error("koFull: audited multiplication is not associative"); fi;
            od;
        od;
    od;
    return rec(status:="computed",normalForms:=forms,normalFormCoordinates:=choices,
        table:=table,identity:=identity,inverses:=inverses,productWitnesses:=proofs,
        commutativityVerified:=true,associativityVerified:=true,
        freeGeneratorIndices:=free,scope:="finite torsion normal forms; free quotients split in the abelian abutment category");
end);
