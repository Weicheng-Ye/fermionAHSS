# Exact named twists, transported from group-coordinate bar cocycles. Numeric
# H2 coordinates are computed only after the extension class has been fixed.
BindGlobal("KOAHSSBatchPaperTwists",function(model,backend,H1,H2,spec,orientation)
    local result,bar,row,used,factorMap,m,j,index,sFunction,wFunction,s,w,
          sClass,wClass,sc,wc,existing,kind,distinguished,zeroS,zeroW,pinTwist;
    zeroS:=List([1..backend.dimension(1)],i->0);
    zeroW:=List([1..backend.dimension(2)],i->0);
    if model.kind="spacegroup" then
        if orientation=fail then Error("koAHSS batch: paper space-group twist needs its orientation character"); fi;
        s:=H1.represent(H1.class(orientation.cochain));
        pinTwist:=KOAHSSBatchSpaceGroupPinTwist(model,backend,H2);
        # These are two requested families, even when their classes coincide
        # after pulling back from the linear image to the full space group.
        return [rec(s_cochain:=s,
            omega_cochain:=zeroW,s_coordinates:=orientation.coordinates,
            omega_coordinates:=Exponents(One(H2.group)),paper_source_rows:=[],
            spacegroup_twist:="zero"),
            rec(s_cochain:=s,omega_cochain:=pinTwist.omega_cochain,
                s_coordinates:=orientation.coordinates,
                omega_coordinates:=pinTwist.omega_coordinates,paper_source_rows:=[],
                spacegroup_twist:="w2_plus_w1_squared",twist_audit:=pinTwist.twist_audit)];
    fi;
    if model.kind<>"abelian" or not IsBound(model.productCoordinates)
       or not IsBound(spec.source_rows) then
        Error("koAHSS batch: paper twists require explicit cyclic product coordinates and source rows");
    fi;
    bar:=KOAHSS_NaturalBarTransport(model.resolution); result:=[];
    for row in spec.source_rows do
        # Match occurrences, not just orders: repeated factors are distinct.
        # Factor 1 in paper order remains distinguished after canonical sorting.
        used:=[]; factorMap:=[];
        for m in row.quotient_factors_in_paper_order do
            if m=1 then Add(factorMap,0);
            else
                index:=First([1..Length(spec.factors)],j->
                    spec.factors[j]=m and not j in used);
                if index=fail then Error("koAHSS batch: paper factors do not match the actual group"); fi;
                Add(used,index); Add(factorMap,index);
            fi;
        od;
        if Length(used)<>Length(spec.factors) then
            Error("koAHSS batch: paper row omits an actual cyclic factor");
        fi;
        if row.paper_s="0" then s:=zeroS;
        elif row.paper_s="x" then
            if Length(factorMap)<>1 or row.quotient_factors_in_paper_order<>[2] then
                Error("koAHSS batch: time reversal x must name the paper's C2 factor");
            fi;
            distinguished:=factorMap[1];
            sFunction:=simplex->model.productCoordinates(simplex[1]^-1*simplex[2])[distinguished] mod 2;
            s:=List(bar.project(1,sFunction,0),x->x mod 2);
        else Error("koAHSS batch: unknown paper sign character"); fi;
        kind:=row.paper_omega_kind;
        if kind="zero" then w:=zeroW;
        elif kind="cyclic_carry" then
            m:=row.paper_omega_cyclic_order;
            if IsEmpty(factorMap) or factorMap[1]=0
               or row.quotient_factors_in_paper_order[1]<>m then
                Error("koAHSS batch: cyclic extension must name the first paper-order factor");
            fi;
            distinguished:=factorMap[1];
            wFunction:=function(simplex)
                local a,b;
                a:=model.productCoordinates(simplex[1]^-1*simplex[2])[distinguished];
                b:=model.productCoordinates(simplex[2]^-1*simplex[3])[distinguished];
                return QuoInt(a+b,m) mod 2;
            end;
            w:=List(bar.project(2,wFunction,0),x->x mod 2);
        elif kind="quaternion" then
            if row.quotient_factors_in_paper_order<>[2,2] then
                Error("koAHSS batch: quaternion extension needs the named C2 x C2 quotient");
            fi;
            wFunction:=function(simplex)
                local a,b,x,y;
                a:=model.productCoordinates(simplex[1]^-1*simplex[2]);
                b:=model.productCoordinates(simplex[2]^-1*simplex[3]);
                x:=factorMap[1]; y:=factorMap[2];
                return (a[x]*b[x]+a[x]*b[y]+a[y]*b[y]) mod 2;
            end;
            w:=List(bar.project(2,wFunction,0),x->x mod 2);
        else Error("koAHSS batch: unknown paper central extension"); fi;
        if ForAny(backend.coboundary(1,s,false),x->x mod 2<>0)
           or ForAny(backend.coboundary(2,w,false),x->x mod 2<>0) then
            Error("koAHSS batch: projected paper twists are not cocycles");
        fi;
        sClass:=H1.class(s); wClass:=H2.class(w);
        sc:=Exponents(sClass); wc:=Exponents(wClass);
        existing:=First(result,t->t.s_coordinates=sc and t.omega_coordinates=wc);
        if existing=fail then
            Add(result,rec(s_cochain:=H1.represent(sClass),
                omega_cochain:=H2.represent(wClass),s_coordinates:=sc,
                omega_coordinates:=wc,paper_source_rows:=[row]));
        else Add(existing.paper_source_rows,row); fi;
    od;
    return result;
end);
