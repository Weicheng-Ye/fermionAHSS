# Exact pullback of w2+w1^2 along the linear representation of a full space
# group. For a rank-three representation V, det(V) tensor V is oriented and
# w2(det(V) tensor V)=w2(V)+w1(V)^2. Its Spin extension therefore represents
# the requested class, including on orientation-reversing elements.

BindGlobal("KOAHSSBatchSpinMinusCocycle",function(generators)
    local identity,elements,cursor,a,b,next,gram,basis,vector,j,diagonal,
          change,inverseChange,oriented,bits,masks,multipliers,i,k,sign,
          overlap,multiply,reverse,lifts,rotation,remaining,lift,reflections,
          column,norm,reflection,cliffordVector,first,reversed,product,
          expected,table,position,scalar,result;
    identity:=IdentityMat(3);
    if not IsList(generators) or not ForAll(generators,a->
        IsMatrix(a) and Length(a)=3 and ForAll(a,row->Length(row)=3)
        and ForAll(Flat(a),IsRat) and DeterminantMat(a) in [-1,1]) then
        Error("koAHSS batch: spin cocycle needs rational rank-three linear generators");
    fi;
    # Crystallographic point groups in dimension three have order at most 48.
    # Enumerate explicitly so an invalid, infinite image cannot hang GAP.
    elements:=[identity]; cursor:=1;
    while cursor<=Length(elements) do
        for a in generators do
            next:=elements[cursor]*a;
            if not next in elements then
                Add(elements,next);
                if Length(elements)>48 then
                    Error("koAHSS batch: crystallographic linear image exceeds order 48");
                fi;
            fi;
        od;
        cursor:=cursor+1;
    od;
    Sort(elements);
    gram:=Sum(elements,a->TransposedMat(a)*a);
    if not ForAll(elements,a->TransposedMat(a)*gram*a=gram) then
        Error("koAHSS batch: linear image does not preserve the averaged metric");
    fi;
    # Rational Gram-Schmidt: no square roots or floating point approximations.
    basis:=[]; diagonal:=[];
    for i in [1..3] do
        vector:=ShallowCopy(identity[i]);
        for j in [1..Length(basis)] do
            vector:=vector-(vector*gram*basis[j])/diagonal[j]*basis[j];
        od;
        Add(basis,vector); Add(diagonal,vector*gram*vector);
    od;
    if not ForAll(diagonal,d->d>0) then
        Error("koAHSS batch: averaged metric is not positive definite");
    fi;
    change:=TransposedMat(basis); inverseChange:=change^-1;
    oriented:=Set(List(elements,a->DeterminantMat(a)*a));

    # Clifford monomials use bit masks 0,...,7 and e_i^2=-diagonal[i].
    bits:=List([0..7],a->List([0..2],i->QuoInt(a,2^i) mod 2));
    masks:=List([1..8],a->[]); multipliers:=List([1..8],a->[]);
    for a in [1..8] do
        for b in [1..8] do
            sign:=0;
            for i in [1..3] do
                for j in [1..i-1] do sign:=sign+bits[a][i]*bits[b][j]; od;
            od;
            overlap:=Filtered([1..3],i->bits[a][i]=1 and bits[b][i]=1);
            Add(masks[a],1+Sum([1..3],i->((bits[a][i]+bits[b][i]) mod 2)*2^(i-1)));
            Add(multipliers[a],(-1)^(sign+Length(overlap))*Product(overlap,i->diagonal[i]));
        od;
    od;
    multiply:=function(left,right)
        local answer,a,b;
        answer:=List([1..8],i->0);
        for a in [1..8] do
            for b in [1..8] do
                if left[a]<>0 and right[b]<>0 then
                    answer[masks[a][b]]:=answer[masks[a][b]]
                        +left[a]*right[b]*multipliers[a][b];
                fi;
            od;
        od;
        return answer;
    end;
    reverse:=q->List([1..8],i->(-1)^QuoInt(Sum(bits[i])*(Sum(bits[i])-1),2)*q[i]);
    lifts:=[];
    for rotation in oriented do
        remaining:=inverseChange*rotation*change;
        lift:=[1,0,0,0,0,0,0,0]; reflections:=0;
        # If R_v sends B e_i to e_i, replacing B by R_v B fixes e_i.
        # Thus B=R_v1 ... R_vr and v1 ... vr is its Clifford lift.
        for i in [1..3] do
            column:=List(remaining,row->row[i]); vector:=column-identity[i];
            if vector<>[0,0,0] then
                norm:=Sum([1..3],j->diagonal[j]*vector[j]^2);
                reflection:=List([1..3],j->List([1..3],k->
                    identity[j][k]-2*vector[j]*diagonal[k]*vector[k]/norm));
                remaining:=reflection*remaining; reflections:=reflections+1;
                cliffordVector:=[0,vector[1],vector[2],0,vector[3],0,0,0];
                lift:=multiply(lift,cliffordVector);
            fi;
        od;
        if remaining<>identity or IsOddInt(reflections) then
            Error("koAHSS batch: oriented image has no verified even reflection decomposition");
        fi;
        first:=PositionProperty(lift,x->x<>0); lift:=lift/lift[first];
        reversed:=reverse(lift); product:=multiply(lift,reversed); norm:=product[1];
        if norm<=0 or product<>Concatenation([norm],List([1..7],i->0)) then
            Error("koAHSS batch: Clifford lift has invalid norm");
        fi;
        for i in [1..3] do
            cliffordVector:=List([1..8],j->0); cliffordVector[2^(i-1)+1]:=1;
            product:=multiply(multiply(lift,cliffordVector),reversed)/norm;
            column:=List(inverseChange*rotation*change,row->row[i]);
            expected:=[0,column[1],column[2],0,column[3],0,0,0];
            if product<>expected then Error("koAHSS batch: Clifford lift has the wrong linear action"); fi;
        od;
        Add(lifts,lift);
    od;
    # The genuine Spin section is q/sqrt(q*reverse(q)), with the positive
    # square root. Norms multiply, so if q_a*q_b=lambda*q_ab, the normalized
    # product differs by sign(lambda). This computes its central F2 cocycle
    # exactly without adjoining square roots.
    table:=List(oriented,a->[]);
    for i in [1..Length(oriented)] do
        for j in [1..Length(oriented)] do
            position:=Position(oriented,oriented[i]*oriented[j]);
            product:=multiply(lifts[i],lifts[j]);
            first:=PositionProperty(lifts[position],x->x<>0);
            scalar:=product[first]/lifts[position][first];
            if scalar=0 or product<>scalar*lifts[position] then
                Error("koAHSS batch: Clifford lift product is not a scalar multiple");
            fi;
            if scalar<0 then Add(table[i],1); else Add(table[i],0); fi;
        od;
    od;
    result:=rec(linear_elements:=elements,oriented_elements:=oriented,
        gram_matrix:=gram,orthogonal_basis:=change,diagonal_metric:=diagonal,
        clifford_lifts:=lifts,cocycle_table:=table);
    result.value:=function(a,b)
        local i,j;
        i:=Position(oriented,DeterminantMat(a)*a);
        j:=Position(oriented,DeterminantMat(b)*b);
        if i=fail or j=fail then Error("koAHSS batch: linear cocycle input outside the point image"); fi;
        return table[i][j];
    end;
    return result;
end);

BindGlobal("KOAHSSBatchSpaceGroupPinTwist",function(model,backend,H2)
    local linear,generators,spin,bar,w,omega,class,audit;
    if model.kind<>"spacegroup" then
        Error("koAHSS batch: w2+w1^2 pullback needs a full space-group model");
    fi;
    linear:=function(g)
        if not IsMatrix(g) or Length(g)<>4 or g[4]<>[0,0,0,1] then
            Error("koAHSS batch: expected a left-action affine resolution matrix");
        fi;
        return List(g{[1..3]},row->row{[1..3]});
    end;
    generators:=List(GeneratorsOfGroup(model.group),linear);
    spin:=KOAHSSBatchSpinMinusCocycle(generators);
    bar:=KOAHSS_NaturalBarTransport(model.resolution);
    w:=simplex->spin.value(linear(simplex[1]^-1*simplex[2]),
        linear(simplex[2]^-1*simplex[3]));
    omega:=List(bar.project(2,w,0),x->x mod 2);
    if ForAny(backend.coboundary(2,omega,false),x->x mod 2<>0) then
        Error("koAHSS batch: projected space-group w2+w1^2 is not a cocycle");
    fi;
    class:=H2.class(omega);
    # JSON stores rational certificate entries as strings, without rounding.
    audit:=rec(class:="w2+w1^2",representation:="det(A)*A",
        construction:="exact rational Clifford lifts of the oriented linear image; Spin(3) pullback to the full affine group",
        point_group_order:=Length(spin.linear_elements),
        oriented_image_order:=Length(spin.oriented_elements),
        gram_matrix:=List(spin.gram_matrix,row->List(row,String)),
        orthogonal_basis:=List(spin.orthogonal_basis,row->List(row,String)),
        diagonal_metric:=List(spin.diagonal_metric,String),
        oriented_elements:=List(spin.oriented_elements,a->List(a,row->List(row,String))),
        clifford_lifts:=List(spin.clifford_lifts,row->List(row,String)),
        cocycle_table:=spin.cocycle_table,
        resolution_pullback:="homogeneous bar project followed by the actual H2 class representative");
    return rec(omega_cochain:=H2.represent(class),omega_coordinates:=Exponents(class),
        twist_audit:=audit);
end);
