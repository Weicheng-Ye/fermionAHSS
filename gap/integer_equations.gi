# Exact integral affine systems, x * matrix = rhs.  Unlike phase equations,
# witnesses for integral coboundaries must not be reduced modulo a modulus.
# The same differential occurs with many right-hand sides in lift and gauge
# searches. Cache only its exact Smith preparation, never a solved affine
# family. Keys are immutable value snapshots, so caller mutation cannot make
# a preparation apply to a different matrix. The column count distinguishes
# empty matrices with different target ranks.
BindGlobal("KOAHSS_PrepareIntegerSystem",CallFuncList(function()
    local cache,cells,maxEntries,maxCells;
    cache:=[]; cells:=0; maxEntries:=8; maxCells:=2000000;
    return function(matrix,cols)
        local rows,index,prepared,snapshot,snf,size,old;
        rows:=Length(matrix);
        index:=PositionProperty(cache,item->item.cols=cols and item.matrix=matrix);
        if index<>fail then
            prepared:=Remove(cache,index); Add(cache,prepared); return prepared;
        fi;
        snapshot:=StructuralCopy(matrix);
        if rows=0 or cols=0 then
            snf:=fail; size:=rows;
        else
            snf:=SmithNormalFormIntegerMatTransforms(snapshot);
            # Discard the implementation's auxiliary transforms. These are
            # the only fields used by the affine solver below.
            snf:=rec(rank:=snf.rank,normal:=snf.normal,
                rowtrans:=snf.rowtrans,coltrans:=snf.coltrans);
            size:=2*rows*cols+rows^2+cols^2;
        fi;
        prepared:=rec(matrix:=snapshot,rows:=rows,cols:=cols,snf:=snf,cells:=size);
        MakeImmutable(prepared);
        if size<=maxCells then
            while Length(cache)>=maxEntries or cells+size>maxCells do
                old:=Remove(cache,1); cells:=cells-old.cells;
            od;
            Add(cache,prepared); cells:=cells+size;
        fi;
        return prepared;
    end;
end,[]));

InstallGlobalFunction(koAHSSSolveIntegerSystem, function(matrix,rhs)
    local rows, cols, prepared, snf, transformed, y, i, generators, solution;
    if not IsList(matrix) or not IsList(rhs) or not ForAll(rhs,IsInt) then
        Error("koAHSS: integer equation matrix and rhs must contain integers");
    fi;
    rows:=Length(matrix); cols:=Length(rhs);
    if not ForAll(matrix,r -> IsList(r) and Length(r)=cols and ForAll(r,IsInt)) then
        Error("koAHSS: integer equation matrix has incorrect dimensions");
    fi;
    prepared:=KOAHSS_PrepareIntegerSystem(matrix,cols);
    if rows=0 then
        if ForAny(rhs,x -> x<>0) then return fail; fi;
        return rec(particular:=[],homogeneousGenerators:=[],rank:=0,coefficientRing:=Integers);
    fi;
    if cols=0 then
        return rec(particular:=List([1..rows],i->0),homogeneousGenerators:=IdentityMat(rows),
                   rank:=0,coefficientRing:=Integers);
    fi;
    snf:=prepared.snf;
    transformed:=rhs*snf.coltrans;
    y:=List([1..rows],i->0);
    for i in [1..snf.rank] do
        if transformed[i] mod snf.normal[i][i]<>0 then return fail; fi;
        y[i]:=QuoInt(transformed[i],snf.normal[i][i]);
    od;
    if ForAny([snf.rank+1..cols],i->transformed[i]<>0) then return fail; fi;
    solution:=y*snf.rowtrans;
    generators:=List([snf.rank+1..rows],i->ShallowCopy(snf.rowtrans[i]));
    if solution*matrix<>rhs or ForAny(generators,v->ForAny(v*matrix,x->x<>0)) then
        Error("koAHSS: integer equation verification failed");
    fi;
    return rec(particular:=solution,homogeneousGenerators:=generators,
               rank:=snf.rank,coefficientRing:=Integers);
end);
