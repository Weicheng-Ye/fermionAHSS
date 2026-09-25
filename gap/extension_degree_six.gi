# Exact linear coordinates for the bounded production degree-six section.
# Copyright (c) 2026 koAHSS contributors. Distributed under the MIT license.
BindGlobal("KOAHSS_ExtensionDegreeSixData",function(model)
    local dimA,dimB,dimC,dimNext,matrix,snf,rank,order,basis,inverse,
          delta,packed,answer;
    dimA:=model.dimension(3); dimB:=model.dimension(4);
    dimC:=model.dimension(5); dimNext:=model.dimension(6);
    matrix:=model.matrix(3,true);
    if Length(matrix)<>dimA or not ForAll(matrix,
        row->IsList(row) and Length(row)=dimB and ForAll(row,IsInt)) then
        Error("koFull: degree-six integral differential has incorrect dimensions");
    fi;
    # GAP differentials act on rows. If U M V=S, the final rows of U
    # span the integral kernel; transpose them into the first basis columns.
    if dimA=0 then
        rank:=0; basis:=[]; inverse:=[];
    elif dimB=0 then
        rank:=0; basis:=IdentityMat(dimA); inverse:=IdentityMat(dimA);
    else
        snf:=SmithNormalFormIntegerMatTransforms(matrix); rank:=snf.rank;
        order:=Concatenation([rank+1..dimA],[1..rank]);
        basis:=TransposedMat(snf.rowtrans{order}); inverse:=Inverse(basis);
        if not ForAll(inverse,row->ForAll(row,IsInt)) then
            Error("koFull: degree-six kernel basis is not unimodular");
        fi;
    fi;
    delta:=List([1..dimB],j->List([1..dimA],i->matrix[i][j]));
    packed:=function(n,rows,cols)
        local m;
        m:=model.matrix(n,false);
        if Length(m)<>cols or not ForAll(m,
            row->IsList(row) and Length(row)=rows and ForAll(row,IsInt)) then
            Error("koFull: degree-six binary differential has incorrect dimensions");
        fi;
        return List([1..rows],j->Sum([1..cols],i->(m[i][j] mod 2)*2^(i-1)));
    end;
    answer:=rec(delta:=delta,basis:=basis,inverse:=inverse,
        cycleRank:=dimA-rank,dimB:=dimB,dimC:=dimC,
        deltaB:=packed(4,dimC,dimB),deltaC:=packed(5,dimNext,dimC),
        maxSectionCandidates:=4096);
    return answer;
end);
