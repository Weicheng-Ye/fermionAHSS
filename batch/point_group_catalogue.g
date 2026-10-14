# The 32 crystallographic point groups, in Table I order of
# https://arxiv.org/abs/2204.13558v2 (Zhang, Ning, Qi, and Gu).
# These are faithful three-dimensional spatial representations, not merely
# abstract groups with the same multiplication table.  Matrices act on columns.
#
# Cartesian coordinates suffice except in the trigonal/hexagonal families.
# There we use the lattice basis (e1,e2,e3), with angle(e1,e2)=60 degrees:
# r6(e1)=e2, r6(e2)=e2-e1, and r6(e3)=e3.  This keeps every entry integral.
# KOAHSSBatchSpinMinusCocycle obtains an exact invariant metric by averaging;
# no Euclidean square roots or external CrystCat catalogue are required.
#
# S4 denotes the fourfold roto-reflection point group, not SymmetricGroup(4).
# Td acts on the tetrahedron with vertices (1,1,1), (1,-1,-1),
# (-1,1,-1), (-1,-1,1); its added coordinate swap is a genuine reflection.

BindGlobal("KOAHSSPointGroupCatalogue",function()
    local identity,inversion,horizontal,vertical,c2x,r2,r4,r6,r3,
          hexVertical,hexC2,cycle,swap,entries,add;
    identity:=IdentityMat(3);
    inversion:=-identity;
    horizontal:=[[1,0,0],[0,1,0],[0,0,-1]];
    vertical:=[[1,0,0],[0,-1,0],[0,0,1]];
    c2x:=[[1,0,0],[0,-1,0],[0,0,-1]];
    r2:=[[-1,0,0],[0,-1,0],[0,0,1]];
    r4:=[[0,-1,0],[1,0,0],[0,0,1]];
    r6:=[[0,-1,0],[1,1,0],[0,0,1]];
    r3:=r6^2;
    hexVertical:=[[1,1,0],[0,-1,0],[0,0,1]];
    hexC2:=hexVertical*horizontal;
    cycle:=[[0,0,1],[1,0,0],[0,1,0]];
    swap:=[[0,1,0],[1,0,0],[0,0,1]];
    entries:=[];
    add:=function(symbol,hermannMauguin,aliases,order,generators)
        Add(entries,rec(id:=symbol,schoenflies:=symbol,
            hermannMauguin:=hermannMauguin,aliases:=aliases,order:=order,
            generators:=StructuralCopy(generators)));
    end;

    add("C1", "1", [], 1, [identity]);
    add("Ci", "-1", ["S2"], 2, [inversion]);
    add("C2", "2", [], 2, [r2]);
    add("Cs", "m", ["C1h"], 2, [horizontal]);
    add("C2h", "2/m", [], 4, [r2,horizontal]);
    add("D2", "222", ["V"], 4, [r2,c2x]);
    add("C2v", "mm2", [], 4, [r2,vertical]);
    add("D2h", "mmm", ["Vh"], 8, [r2,c2x,horizontal]);
    add("C4", "4", [], 4, [r4]);
    add("S4", "-4", [], 4, [r4*horizontal]);
    add("C4h", "4/m", [], 8, [r4,horizontal]);
    add("D4", "422", [], 8, [r4,c2x]);
    add("C4v", "4mm", [], 8, [r4,vertical]);
    add("D2d", "-42m", ["Vd","-4m2"], 8, [r4*horizontal,c2x]);
    add("D4h", "4/mmm", [], 16, [r4,c2x,horizontal]);
    add("C3", "3", [], 3, [r3]);
    add("S6", "-3", ["C3i"], 6, [r6*horizontal]);
    add("D3", "32", [], 6, [r3,hexC2]);
    add("C3v", "3m", [], 6, [r3,hexVertical]);
    add("D3d", "-3m", [], 12, [r3,hexC2,inversion]);
    add("C6", "6", [], 6, [r6]);
    add("C3h", "-6", [], 6, [r3,horizontal]);
    add("C6h", "6/m", [], 12, [r6,horizontal]);
    add("D6", "622", [], 12, [r6,hexC2]);
    add("C6v", "6mm", [], 12, [r6,hexVertical]);
    add("D3h", "-6m2", ["-62m"], 12, [r3,hexC2,horizontal]);
    add("D6h", "6/mmm", [], 24, [r6,hexC2,horizontal]);
    add("T", "23", [], 12, [cycle,c2x]);
    add("Th", "m-3", [], 24, [cycle,c2x,inversion]);
    add("Td", "-43m", [], 24, [cycle,c2x,swap]);
    add("O", "432", [], 24, [cycle,r4]);
    add("Oh", "m-3m", [], 48, [cycle,r4,inversion]);
    return entries;
end);
