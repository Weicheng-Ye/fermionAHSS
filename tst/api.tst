# Detailed group results reuse the exact backend and preserve legacy outputs.
gap> apiGroup := CyclicGroup(2);;
gap> apiExpected := [[[0],[],[2]],[[],[]],[[2]],[],[]];;
gap> apiResult := koAHSS(apiGroup,0,0,1,rec(details:=true));;
gap> Assert(0, apiResult.kind="koAHSSResult" and apiResult.status="computed");
gap> Assert(0, apiResult.maxDegree=1 and apiResult.computedThrough=6);
gap> Assert(0, apiResult.pages.kind="koAHSSPages" and apiResult.pages.pageNumbers=[6]);
gap> Assert(0, apiResult.pages.tables=[apiExpected] and Length(apiResult.pageData)=1);
gap> Assert(0, apiResult.certified_ko=false and apiResult.scope="five-row-associated-graded");
gap> Assert(0, IsIdenticalObj(apiResult._context.backend.twists,apiResult.twists));
gap> Assert(0, apiResult._context.kind="koAHSSContext");
gap> Assert(0, IsIdenticalObj(apiResult._context.group,apiGroup));
gap> Assert(0, IsIdenticalObj(apiResult._context.getCell(6,2,-4),apiResult.pageData[1][1][3]));
gap> apiCell := apiResult._context.getCell(6,2,-4);;
gap> Assert(0, ForAll(GeneratorsOfGroup(apiCell.group), g -> apiCell.project(apiCell.lift(g))=g));
gap> Assert(0, IsBound(apiCell.previous) and IsBound(apiCell.quotientMap));
gap> Assert(0, koAHSSFormat(apiResult)=koAHSSFormat(apiExpected));
gap> Assert(0, koAHSSFormat(apiResult.pages)=koAHSSFormat(apiExpected));
gap> apiPages := koAHSS(apiGroup,0,0,1,5,rec(details:=true));;
gap> Assert(0, apiPages.computedThrough=6 and apiPages.pages.pageNumbers=[2..6]);
gap> Assert(0, apiPages.pages.tables=List([1..5],i->apiExpected));
gap> Assert(0, koAHSSFormat(apiPages)=koAHSSFormat(apiPages.pages.tables));
gap> Assert(0, koAHSSFormat(apiPages,4)=koAHSSFormat(apiExpected,4));
gap> apiEarly := koAHSS(apiGroup,0,0,-1,1,rec(details:=true));;
gap> Assert(0, apiEarly.computedThrough=2 and apiEarly.pages.pageNumbers=[2]);
gap> Assert(0, koAHSS(apiGroup,0,0,-1,rec(details:=false))=[[[0]],[],[],[],[]]);
gap> Assert(0, koAHSS(apiGroup,0,0,-1,1,rec(details:=false))=[[[[0]],[],[],[],[]]]);
gap> # Low-level wrappers return a sole cell table, or an E2-first list.
gap> apiSpace := koAHSSHAPSpace(ResolutionFiniteGroup(apiGroup,3),koAHSSNaturalOperations());;
gap> apiCells := koAHSSPageData(apiSpace,0,0,-1);;
gap> Assert(0, Length(apiCells)=5 and IsRecord(apiCells[1][1]) and apiCells[1][1].page=6);
gap> apiCellPages := koAHSSPageData(apiSpace,0,0,-1,1);;
gap> Assert(0, Length(apiCellPages)=1 and Length(apiCellPages[1])=5 and apiCellPages[1][1][1].page=2);
gap> koAHSS(apiGroup,0,0,1,rec(detail:=true));
Error, koAHSS: unknown option; supported options are details
gap> koAHSS(apiGroup,0,0,1,rec(details:=1));
Error, koAHSS: details must be true or false
gap> koAHSS(apiGroup,0,0,1,rec(details:=true),5);
Error, koAHSS: n counts pages beginning with E2 and must be in [1..5]
gap> koAHSS(apiGroup,0,0,1,5,true);
Error, koAHSS: the final options argument must be a record
