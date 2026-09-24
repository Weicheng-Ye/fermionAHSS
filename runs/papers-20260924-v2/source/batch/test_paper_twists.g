# Bounded checks of named extensions and their actual resolution coordinates.
Read("load.g");;
Read("batch/koahss_groups.g");;
Read("batch/koahss_spacegroup_twists.g");;
Read("batch/koahss_paper_twists.g");;

KOAHSSPaperTwistTest:=function()
    local checks,check,row,run,answer,model,backend,twists,coordinates,
          elements,g,h,bar,w,pair,first,second,cyclicSum,identity;
    checks:=0;
    check:=function(condition,message)
        checks:=checks+1;
        if not condition then Error("paper twist test: ",message); fi;
    end;
    row:=function(number,factors,s,kind)
        local result;
        result:=rec(row:=number,quotient_factors_in_paper_order:=factors,
            paper_s:=s,paper_omega_kind:=kind);
        if kind="cyclic_carry" then result.paper_omega_cyclic_order:=factors[1]; fi;
        return result;
    end;
    run:=function(factors,rows,count)
        local spec,model,backend,H1,H2,twists;
        spec:=rec(kind:="abelian",factors:=factors,source_rows:=rows);
        model:=KOAHSSBatchMakeGroup(spec,3);
        backend:=koAHSSHAPSpace(model.resolution).koAHSS(0,0,3);
        H1:=backend.data(1,-1); H2:=backend.data(2,-1);
        twists:=KOAHSSBatchPaperTwists(model,backend,H1,H2,spec,fail);
        check(Length(twists)=count,"wrong number of distinct paper twists");
        check(Sum(twists,t->Length(t.paper_source_rows))=Length(rows),
            "lost a source row during deduplication");
        return rec(model:=model,backend:=backend,twists:=twists);
    end;
    answer:=run([2],[row(1,[2],"0","zero"),row(3,[2],"0","zero"),
        row(12,[2],"0","cyclic_carry"),row(19,[2],"x","zero"),
        row(20,[2],"x","cyclic_carry")],4);
    check(Length(answer.twists[1].paper_source_rows)=2,"split C2 rows not merged");
    check(Length(Set(List(answer.twists,t->[t.s_coordinates,t.omega_coordinates])))=4,
        "C2 must have four different paper pairs");

    answer:=run([2,4],[row(5,[2,4],"0","zero"),
        row(14,[2,4],"0","cyclic_carry"),
        row(15,[4,2],"0","cyclic_carry")],3);
    model:=answer.model; twists:=answer.twists;
    bar:=KOAHSS_NaturalBarTransport(model.resolution); identity:=One(model.group);
    first:=First(Elements(model.group),g->model.productCoordinates(g)=[1,0]);
    second:=First(Elements(model.group),g->model.productCoordinates(g)=[0,1]);
    cyclicSum:=function(twist,g,m)
        local value,j;
        value:=bar.lift(2,twist.omega_cochain,0);
        return Sum([0..m-1],j->value([identity,g^j,g^(j+1)])) mod 2;
    end;
    check(cyclicSum(twists[2],first,2)=1 and cyclicSum(twists[2],second,4)=0,
        "row 14 carry is not on the paper C2 factor");
    check(cyclicSum(twists[3],first,2)=0 and cyclicSum(twists[3],second,4)=1,
        "row 15 carry is not on the paper C4 factor after sorting");

    answer:=run([2,2],[row(4,[2,2],"0","zero"),
        row(13,[2,2],"0","cyclic_carry"),
        row(21,[2,2],"0","quaternion")],3);
    model:=answer.model; twists:=answer.twists;
    bar:=KOAHSS_NaturalBarTransport(model.resolution); identity:=One(model.group);
    w:=bar.lift(2,twists[3].omega_cochain,0);
    for g in Filtered(Elements(model.group),g->g<>identity) do
        check(w([identity,g,identity]) mod 2=1,
            "quaternion extension must square to fermion parity on all three involutions");
    od;
    first:=First(Elements(model.group),g->model.productCoordinates(g)=[1,0]);
    second:=First(Elements(model.group),g->model.productCoordinates(g)=[0,1]);
    w:=bar.lift(2,twists[2].omega_cochain,0);
    check(w([identity,first,identity]) mod 2=1 and w([identity,second,identity]) mod 2=0,
        "repeated C2 factor occurrence was not retained");

    for pair in [[[6],3],[[3],2],[[],0]] do
        answer:=run(pair[1],[row(pair[2],pair[1],"0","zero")],1);
        model:=answer.model; elements:=Elements(model.group);
        check(Length(Set(List(elements,model.productCoordinates)))=Size(model.group),
            "product coordinates must be bijective, including C6 and the trivial group");
        for g in elements do
            for h in elements do
                coordinates:=List([1..Length(pair[1])],j->
                    (model.productCoordinates(g)[j]+model.productCoordinates(h)[j]) mod pair[1][j]);
                check(model.productCoordinates(g*h)=coordinates,"product coordinates are not additive");
            od;
        od;
    od;
    Print(checks," exact paper-twist construction checks passed.\n");
end;;
KOAHSSPaperTwistTest();;
QUIT;
