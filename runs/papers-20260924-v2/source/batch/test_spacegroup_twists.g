# Bounded exact tests of the Spin pullback and its full-space-group class.
Read("load.g");;
Read("batch/koahss_groups.g");;
Read("batch/koahss_spacegroup_twists.g");;

KOAHSSSpaceGroupTwistTest:=function()
    local checks,check,verify,identity,rotation,mirror,inversion,spin,
          halfTurnX,halfTurnY,order3,order6,flip,quarterTurn,cycleAxes,model,backend,H2,twist,
          run,sg3,sg4,sg6,sg7,oddSplit,values,g,h,j;
    checks:=0;
    check:=function(condition,message)
        checks:=checks+1;
        if not condition then Error("space-group twist test: ",message); fi;
    end;
    identity:=IdentityMat(3);
    verify:=function(generators,size)
        local spin,elements;
        spin:=KOAHSSBatchSpinMinusCocycle(generators); elements:=spin.linear_elements;
        check(Length(elements)=size,"wrong finite linear image order");
        check(ForAll(elements,g->spin.value(identity,g)=0 and spin.value(g,identity)=0),
            "Spin cocycle is not normalized");
        check(ForAll(spin.oriented_elements,g->ForAll(spin.oriented_elements,h->ForAll(spin.oriented_elements,k->
            (spin.value(h,k)+spin.value(g*h,k)+spin.value(g,h*k)+spin.value(g,h)) mod 2=0))),
            "Spin signs fail the group cocycle identity");
        return spin;
    end;
    spin:=verify([],1);
    rotation:=DiagonalMat([-1,-1,1]); spin:=verify([rotation],2);
    check(spin.value(rotation,rotation)=1,"C2 rotation must have a lift squaring to -1");
    mirror:=DiagonalMat([-1,1,1]); spin:=verify([mirror],2);
    check(spin.value(mirror,mirror)=1,"mirror must detect w1 squared in the Pin-minus class");
    inversion:=-identity; spin:=verify([inversion],2);
    check(spin.value(inversion,inversion)=0,"rank-three inversion has w2+w1 squared equal to zero");
    halfTurnX:=DiagonalMat([1,-1,-1]); halfTurnY:=DiagonalMat([-1,1,-1]);
    spin:=verify([halfTurnX,halfTurnY],4);
    check(ForAll(Filtered(spin.linear_elements,g->g<>identity),g->spin.value(g,g)=1),
        "D2 spin cover must be quaternion: each involution lifts to order four");
    check((spin.value(halfTurnX,halfTurnY)+spin.value(halfTurnY,halfTurnX)) mod 2=1,
        "orthogonal half-turn lifts must anticommute");

    order3:=[[0,-1,0],[1,-1,0],[0,0,1]];
    spin:=verify([order3],3);
    check(spin.gram_matrix[1][2]<>0,"hexagonal test must use a nonorthogonal lattice metric");
    # On C3 every F2 cocycle splits. Verify this one by exhibiting a 1-cochain.
    oddSplit:=false;
    for values in Tuples([0,1],2) do
        values:=Concatenation([0],values);
        if ForAll([0..2],g->ForAll([0..2],h->
            spin.value(order3^g,order3^h)=(values[g+1]+values[h+1]+values[((g+h) mod 3)+1]) mod 2)) then
            oddSplit:=true;
        fi;
    od;
    check(oddSplit,"odd-order rotation extension did not split");
    order6:=[[0,-1,0],[1,1,0],[0,0,1]]; spin:=verify([order6],6);
    check(Sum([0..5],j->spin.value(order6^j,order6)) mod 2=1,
        "nonorthogonal C6 rotation must lift to order twelve");
    flip:=[[0,1,0],[1,0,0],[0,0,-1]]; spin:=verify([order6,flip],12);
    check(spin.value(flip,flip)=1,"hexagonal dihedral half-turn must lift to order four");
    quarterTurn:=[[0,-1,0],[1,0,0],[0,0,1]];
    cycleAxes:=[[0,1,0],[0,0,1],[1,0,0]];
    spin:=verify([quarterTurn,cycleAxes,inversion],48);
    check(Length(spin.oriented_elements)=24,"full cubic point group has oriented image of order 24");
    check(Sum([0..3],j->spin.value(quarterTurn^j,quarterTurn)) mod 2=1,
        "cubic quarter-turn must lift to order eight");

    run:=function(number)
        local model,backend,H2,twist;
        model:=KOAHSSBatchMakeGroup(rec(kind:="spacegroup",number:=number),3);
        backend:=koAHSSHAPSpace(model.resolution).koAHSS(0,0,3);
        H2:=backend.data(2,-1);
        twist:=KOAHSSBatchSpaceGroupPinTwist(model,backend,H2);
        check(ForAll(backend.coboundary(2,twist.omega_cochain,false),x->x mod 2=0),
            "canonical projected space-group class is not a cocycle");
        check(Exponents(H2.class(twist.omega_cochain))=twist.omega_coordinates,
            "space-group class and canonical representative disagree");
        return twist;
    end;
    twist:=run(1);
    check(ForAll(twist.omega_coordinates,x->x=0),"pure translations have a trivial linear twist");
    twist:=run(2);
    check(ForAll(twist.omega_coordinates,x->x=0),"inversion space group has trivial Pin-minus class");
    sg3:=run(3); sg4:=run(4);
    check(ForAny(sg3.omega_coordinates,x->x<>0),"P2 must retain the nontrivial rotation extension");
    check(ForAll(sg4.omega_coordinates,x->x=0),"P21 screw pullback must trivialize the rotation extension");
    sg6:=run(6); sg7:=run(7);
    check(ForAny(sg6.omega_coordinates,x->x<>0),"Pm must retain the nontrivial mirror extension");
    check(ForAll(sg7.omega_coordinates,x->x=0),"Pc glide pullback must trivialize the mirror extension");
    Print(checks," exact space-group twist checks passed (six full affine groups).\n");
end;;
KOAHSSSpaceGroupTwistTest();;
QUIT;
