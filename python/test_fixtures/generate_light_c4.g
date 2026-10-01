# Regenerate from the package root:
# gap -q --quitonbreak python/test_fixtures/generate_light_c4.g
# Exact C4 comparison through degree four, including the nonzero homotopy.
Read("load.g");;
CallFuncList(function()
    local R,backend,tr,setup,records,n,j,t,vertices,matrix,out;
    R:=ResolutionFiniteGroup(CyclicGroup(4),6);
    backend:=koAHSSHAPSpace(R,koAHSSNaturalOperations()).koAHSS([1],[1],5);
    tr:=KOAHSS_ExtensionNormalizedTransport(backend,2);
    Assert(0,tr.status="computed" and tr.tableMode);
    matrix:=function(n,signed)
        return List(IdentityMat(backend.dimension(n)),v->backend.coboundary(n,v,signed));
    end;
    setup:=rec(schema:=1,vertexMode:="table",gMode:="lazy",
        multiplication:=List(tr.elements,g->List(tr.elements,h->tr.vertexLabel(g*h))),
        ranks:=List([0..4],backend.dimension),s:=[1],omega:=[1],
        ordinary:=List([0..3],n->matrix(n,false)),signed:=List([0..3],n->matrix(n,true)));
    records:=[];
    for n in [0..4] do
        for j in [0..backend.dimension(n)-1] do
            Add(records,rec(kind:="g",degree:=n,vertices:=j,answer:=tr.chain(n,j)));
        od;
        for t in Tuples([0..3],n) do
            vertices:=Concatenation([0],t);
            if ForAny([1..n],i->vertices[i]=vertices[i+1]) then continue; fi;
            Add(records,rec(kind:="f",degree:=n,vertices:=vertices,answer:=tr.terms("f",vertices)));
            if n<4 then
                Add(records,rec(kind:="h",degree:=n,vertices:=vertices,answer:=tr.terms("h",vertices)));
            fi;
        od;
    od;
    out:=OutputTextFile("python/test_fixtures/light_c4.json",false);
    SetPrintFormattingStatus(out,false);
    WriteAll(out,GapToJsonString(rec(group:="C4",setup:=setup,transport:=records)));
    CloseStream(out);
end,[]);
QUIT;
