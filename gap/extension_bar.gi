# Complete normalized finite-bar model for measured higher stacking relations.
CallFuncList(function()
    local path,slash;
    path:=INPUT_FILENAME();
    if path[1]<>'/' then path:=Filename(DirectoryCurrent(),path); fi;
    slash:=Last(Positions(path,'/'));
    BindGlobal("KOAHSS_EXTENSION_WORKER",Concatenation(path{[1..slash]},"../python/extension_worker.py"));
end,[]);

BindGlobal("KOAHSS_ExtensionBarModel",function(backend,maxDegree)
    local transport,group,elements,unit,nonunit,model,simplices,coordinates,
        evaluate,ctx,stream,request,setup,executable,cache,cacheCount,matrices;
    if not IsBound(backend.naturalBar) then
        return rec(status:="unresolved",reason:="higher stacking requires a normalized finite group-bar transport");
    fi;
    transport:=backend.naturalBar(); group:=transport.resolution!.group;
    if not IsFinite(group) then
        return rec(status:="unresolved",reason:="the complete bar extension model requires a finite group");
    fi;
    if (Size(group)-1)^(maxDegree+2)>8192 then
        return rec(status:="unresolved",reason:="complete bar extension model exceeds 8192 cochains in a required degree");
    fi;
    elements:=AsList(group); unit:=One(group);
    elements:=Concatenation([unit],Filtered(elements,g->g<>unit)); nonunit:=elements{[2..Length(elements)]};
    model:=rec(status:="computed",modelId:="complete-normalized-finite-bar",maxDegree:=maxDegree);
    model.supports:=k->Length(nonunit)^(k+2)<=8192
        and Length(nonunit)^(2*k+3)<=2000000;
    if not model.supports(maxDegree) then
        return rec(status:="unresolved",reason:="complete bar extension matrix exceeds the 2000000-entry budget");
    fi;
    model.dimension:=function(n)
        if n<0 then return 0; fi;
        return Length(nonunit)^n;
    end;
    simplices:=rec();
    model.simplices:=function(n)
        local key,answer,increments,sigma,g;
        if n<0 then return []; fi;
        key:=String(n);
        if not IsBound(simplices.(key)) then
            answer:=[];
            for increments in Tuples(nonunit,n) do
                sigma:=[unit];
                for g in increments do Add(sigma,Last(sigma)*g); od;
                Add(answer,sigma);
            od;
            simplices.(key):=answer;
        fi;
        return simplices.(key);
    end;
    evaluate:=function(vector,sigma)
        local index,j,g,pos;
        index:=0;
        for j in [2..Length(sigma)] do
            g:=sigma[j-1]^-1*sigma[j];
            if g=unit then return 0; fi;
            pos:=Position(nonunit,g); index:=index*Length(nonunit)+pos-1;
        od;
        return vector[index+1];
    end;
    ctx:=KOAHSS_NaturalCochains(backend);
    model.lift:=function(n,vector,signed)
        local sign,values;
        sign:=0; if signed then sign:=backend.twists.s; fi;
        if n<0 then return []; fi;
        values:=List(model.simplices(n),sigma->transport.evaluateR(n,vector,sign,sigma));
        if not signed then values:=List(values,x->x mod 2); fi;
        return values;
    end;
    model.project:=function(n,vector,signed)
        local sign;
        sign:=0; if signed then sign:=backend.twists.s; fi;
        return transport.project(n,sigma->evaluate(vector,sigma),sign);
    end;
    model.s:=List(model.lift(1,backend.twists.s,false),x->x mod 2);
    model.omega:=List(model.lift(2,backend.twists.omega,false),x->x mod 2);
    model.coboundary:=function(n,vector,signed)
        local answer,sigma,sum,j,face,sign;
        if n<0 then return List([1..model.dimension(n+1)],i->0); fi;
        KOAHSS_CC_CheckVector(vector,model.dimension(n),"bar coboundary input");
        answer:=[];
        for sigma in model.simplices(n+1) do
            sum:=0;
            for j in [1..n+2] do
                face:=sigma{Filtered([1..n+2],i->i<>j)}; sign:=(-1)^(j-1);
                if j=1 and signed then sign:=(-1)^evaluate(model.s,sigma{[1,2]}); fi;
                sum:=sum+sign*evaluate(vector,face);
            od;
            Add(answer,sum);
        od;
        return answer;
    end;
    matrices:=rec();
    model.matrix:=function(n,signed)
        local key,i,e;
        key:=Concatenation(String(n),"_",String(signed));
        if not IsBound(matrices.(key)) then
            matrices.(key):=[];
            for i in [1..model.dimension(n)] do
                e:=List([1..model.dimension(n)],j->0); e[i]:=1;
                Add(matrices.(key),model.coboundary(n,e,signed));
            od;
        fi;
        return matrices.(key);
    end;
    model.zero:=k->rec(A:=List([1..model.dimension(k-3)],i->0),
        B:=List([1..model.dimension(k-2)],i->0),C:=List([1..model.dimension(k-1)],i->0),
        D:=List([1..model.dimension(k+1)],i->0));
    executable:=Filename(DirectoriesSystemPrograms(),"python3");
    if executable=fail then Error("koFull: higher stacking requires Python 3"); fi;
    setup:=rec(operation:="setup",multiplication:=List(elements,g->List(elements,h->Position(elements,g*h)-1)),
        s:=model.s,omega:=model.omega);
    stream:=fail; cache:=NewDictionary("",true); cacheCount:=0;
    model.close:=function()
        if stream<>fail then CloseStream(stream); stream:=fail; fi;
    end;
    request:=function(input)
        local encoded,found,line,answer;
        encoded:=GapToJsonString(input); found:=LookupDictionary(cache,encoded);
        if found<>fail then return found; fi;
        if stream=fail then
            stream:=InputOutputLocalProcess(DirectoryCurrent(),executable,["-u",KOAHSS_EXTENSION_WORKER]);
            WriteLine(stream,GapToJsonString(setup)); line:=ReadLine(stream);
            if line=fail then model.close(); Error("koFull: stacking worker failed to initialize"); fi;
            answer:=JsonStringToGap(line);
            if answer.status<>"computed" then
                if answer.status="unresolved" then model.lastFailure:=answer; fi;
                model.close(); Error("koFull: stacking worker setup failed: ",answer);
            fi;
        fi;
        WriteLine(stream,encoded); line:=ReadLine(stream);
        if line=fail then model.close(); Error("koFull: stacking worker stopped before returning an exact result"); fi;
        answer:=JsonStringToGap(line);
        if not IsBound(answer.status) or answer.status<>"computed" then
            if IsBound(answer.status) and answer.status="unresolved" then model.lastFailure:=answer; fi;
            model.close(); Error("koFull: exact stacking model failed: ",answer);
        fi;
        if cacheCount>=256 then cache:=NewDictionary("",true); cacheCount:=0; fi;
        MakeImmutable(answer.state); AddDictionary(cache,encoded,answer.state); cacheCount:=cacheCount+1;
        return answer.state;
    end;
    model.ensureSix:=function()
        if not IsBound(setup.degreeSix) then
            setup.degreeSix:=CallFuncList(ValueGlobal("KOAHSS_ExtensionDegreeSixData"),[model]);
            model.close();
        fi;
    end;
    model.d:=function(k,state)
        if k>=6 then model.ensureSix(); fi;
        return request(rec(operation:="d",degree:=k,state:=state));
    end;
    model.xtimes:=function(k,x,y)
        if k>=6 then model.ensureSix(); fi;
        return request(rec(operation:="xtimes",degree:=k,state:=x,other:=y));
    end;
    return model;
end);
