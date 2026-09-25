# Formatting distinguishes pages, groups, unresolved cells, and omitted positions.
gap> displayTable := [[[0],[],[2]],[[],[]],[[2]],[],[]];;
gap> displayCopy := StructuralCopy(displayTable);;
gap> displayText := koAHSSFormat(displayTable);;
gap> Assert(0, PositionSublist(displayText,"E6") <> fail);
gap> Assert(0, PositionSublist(displayText,"Z/2") <> fail);
gap> Assert(0, PositionSublist(displayText,"-1") < PositionSublist(displayText,"-2"));
gap> Assert(0, PositionSublist(displayText,"-2") < PositionSublist(displayText,"-3"));
gap> Assert(0, PositionSublist(displayText,"-3") < PositionSublist(displayText,"-4"));
gap> Assert(0, displayTable = displayCopy);
gap> displayPages := List([1..5], i -> StructuralCopy(displayTable));;
gap> displayPages[1][1][1] := [7];;
gap> Assert(0, PositionSublist(koAHSSFormat(displayPages,2),"Z/7") <> fail);
gap> Assert(0, PositionSublist(koAHSSFormat(displayPages,6),"Z/7") = fail);
gap> Assert(0, koAHSSFormat(displayPages,6) = koAHSSFormat(displayPages[5]));
gap> Assert(0, ForAll([2..6], r -> PositionSublist(koAHSSFormat(displayPages),Concatenation("E",String(r))) <> fail));
gap> Assert(0, koAHSSFormat([displayTable]) = koAHSSFormat(displayTable,2));
gap> displayTable[1][1] := [4,2,0,2,0];;
gap> displayCopy := StructuralCopy(displayTable);;
gap> Assert(0, PositionSublist(koAHSSFormat(displayTable),"Z^2 + (Z/2)^2 + Z/4") <> fail);
gap> Assert(0, displayTable = displayCopy);
gap> displayTable[1][2] := rec(status:="unresolved",lastKnownInvariants:=[0]);;
gap> Assert(0, PositionSublist(koAHSSFormat(displayTable),"?") <> fail);
gap> Assert(0, PositionSublist(koAHSSFormat([[[]],[],[],[],[]]),"-4") <> fail);
gap> displayWide := List([8,7,6,5,4], n -> List([1..n], i -> [0,2,2,4,8,16]));;
gap> Assert(0, PositionSublist(koAHSSFormat(displayWide),"Z/16") <> fail);
gap> # Tagged page numbers select by recorded label, not list position.
gap> displayTagged := rec(kind:="koAHSSPages",pageNumbers:=[6,2],tables:=[displayPages[5],displayPages[1]]);;
gap> displayTaggedCopy := StructuralCopy(displayTagged);;
gap> Assert(0, koAHSSFormat(displayTagged,6)=koAHSSFormat(displayPages[5],6));
gap> Assert(0, koAHSSFormat(displayTagged,2)=koAHSSFormat(displayPages[1],2));
gap> Assert(0, koAHSSFormat(displayTagged)=Concatenation(koAHSSFormat(displayPages[5],6),"\n",koAHSSFormat(displayPages[1],2)));
gap> Assert(0, koAHSSFormat(rec(kind:="koAHSSResult",pages:=displayTagged))=koAHSSFormat(displayTagged));
gap> Assert(0, koAHSSFormat(rec(kind:="koFullResult",pages:=displayTagged))=koAHSSFormat(displayTagged,6));
gap> Assert(0, koAHSSFormat(rec(kind:="koFullResult",pages:=displayTagged),2)=koAHSSFormat(displayTagged,2));
gap> Assert(0, displayTagged=displayTaggedCopy);
gap> Assert(0, koAHSSFormat(rec(kind:="koAHSSPages",pageNumbers:=[6],tables:=[displayPages[5]]))=koAHSSFormat(displayPages[5]));
gap> koAHSSFormat(rec(kind:="koAHSSPages",pageNumbers:=[6],tables:=[displayPages[5]]),2);
Error, koAHSSFormat: the requested page is absent from result
gap> koAHSSFormat(rec(kind:="koAHSSPages",pageNumbers:=[6,6],tables:=[displayPages[5],displayPages[5]]));
Error, koAHSSFormat: malformed tagged pages; expected distinct page labels in [2..6] and matching valid tables
gap> koAHSSFormat(rec(kind:="koAHSSPages",pageNumbers:=[6],tables:=[]));
Error, koAHSSFormat: malformed tagged pages; expected distinct page labels in [2..6] and matching valid tables
gap> koAHSSFormat(rec(kind:="koAHSSPages",pageNumbers:=[6,2],tables:=[displayPages[5],[[[0]],[],[],[],[]]]));
Error, koAHSSFormat: all pages must have the same display window
gap> koAHSSFormat(rec(kind:="koAHSSResult"));
Error, koAHSSFormat: detailed results must contain tagged koAHSSPages
