# Human-readable view of raw koAHSS/koAHSSpages invariant tables.  Storage
# remains q=-4,...,0; only the display reverses rows to q=0,...,-4.
InstallGlobalFunction(koAHSSFormat, function(arg)
    local result, selected, isCell, isTable, cellText, centered, formatTable,
        tables, pageNumbers, width;
    if not Length(arg) in [1,2] then
        Error("usage: koAHSSFormat(result[, r])");
    fi;
    result := arg[1];
    selected := fail;
    if Length(arg)=2 then
        selected := arg[2];
        if not IsInt(selected) or not selected in [2..6] then
            Error("koAHSSFormat: page r must be an integer in [2..6]");
        fi;
    fi;
    isCell := function(cell)
        if IsRecord(cell) then
            return IsBound(cell.status) and cell.status="unresolved";
        fi;
        return IsList(cell) and IsDenseList(cell) and
            ForAll(cell, x -> IsInt(x) and (x=0 or x>=2));
    end;
    isTable := function(table)
        local columns, row;
        if not IsList(table) or not IsDenseList(table) or Length(table)<>5
            or not ForAll(table, x -> IsList(x) and IsDenseList(x)) then
            return false;
        fi;
        columns := Length(table[1]);
        # Raw tables have the physical cutoffs k=-1,...,6.  Checking the
        # complete shape also distinguishes a table from a five-page list.
        if not columns in [1..8] then return false; fi;
        for row in [1..5] do
            if Length(table[row])<>Maximum(0,columns-row+1)
                or not ForAll(table[row],isCell) then return false; fi;
        od;
        return true;
    end;
    if isTable(result) then
        tables := [result];
        if selected=fail then pageNumbers := [6];
        else pageNumbers := [selected]; fi;
    elif IsList(result) and IsDenseList(result) and Length(result) in [1..5]
        and ForAll(result,isTable) then
        width := Length(result[1][1]);
        if not ForAll(result, table -> Length(table[1])=width) then
            Error("koAHSSFormat: all pages must have the same display window");
        fi;
        if selected=fail then
            tables := result;
            pageNumbers := [2..Length(result)+1];
        else
            if selected>Length(result)+1 then
                Error("koAHSSFormat: the requested page is absent from result");
            fi;
            tables := [result[selected-1]];
            pageNumbers := [selected];
        fi;
    else
        Error("koAHSSFormat: expected a raw five-row table or E2-first list of pages with valid invariant cells");
    fi;
    cellText := function(cell)
        local factors, pair, term;
        if IsRecord(cell) then return "?"; fi;
        if IsEmpty(cell) then return "0"; fi;
        factors := [];
        # Collected sorts a copy and preserves the caller's invariants.
        for pair in Collected(cell) do
            if pair[1]=0 then term := "Z";
            else term := Concatenation("Z/",String(pair[1])); fi;
            if pair[2]>1 then
                if pair[1]<>0 then term := Concatenation("(",term,")"); fi;
                term := Concatenation(term,"^",String(pair[2]));
            fi;
            Add(factors,term);
        od;
        return JoinStringsWithSeparator(factors," + ");
    end;
    centered := function(text,width)
        local left;
        left := QuoInt(width-Length(text),2);
        return Concatenation(ListWithIdenticalEntries(left,' '),text,
            ListWithIdenticalEntries(width-Length(text)-left,' '));
    end;
    formatTable := function(table,page)
        local columns, cells, widths, row, column, labels, line, border,
            lines, unresolved;
        columns := Length(table[1]);
        cells := List([1..5], row -> List([1..columns], column -> "."));
        unresolved := false;
        for row in [1..5] do
            for column in [1..Length(table[6-row])] do
                cells[row][column] := cellText(table[6-row][column]);
                if cells[row][column]="?" then unresolved := true; fi;
            od;
        od;
        labels := List([0..columns-1],String);
        widths := List([1..columns], column -> Maximum(
            Concatenation([Length(labels[column])],
                List(cells,row -> Length(row[column])))));
        border := "+-------+";
        for column in [1..columns] do
            Append(border,ListWithIdenticalEntries(widths[column]+2,'-'));
            Append(border,"+");
        od;
        lines := [Concatenation("E",String(page)),border];
        line := "| q \\ p |";
        for column in [1..columns] do
            Append(line,Concatenation(" ",centered(labels[column],widths[column])," |"));
        od;
        Add(lines,line);
        Add(lines,border);
        for row in [1..5] do
            line := Concatenation("| ",centered(String(1-row),5)," |");
            for column in [1..columns] do
                Append(line,Concatenation(" ",centered(cells[row][column],widths[column])," |"));
            od;
            Add(lines,line);
        od;
        Add(lines,border);
        Add(lines,"0 = zero group; . = outside the displayed window.");
        if unresolved then
            Add(lines,"? = unresolved; inspect the original entry for details.");
        fi;
        return Concatenation(JoinStringsWithSeparator(lines,"\n"),"\n");
    end;
    return JoinStringsWithSeparator(List([1..Length(tables)],
        i -> formatTable(tables[i],pageNumbers[i])),"\n");
end);

InstallGlobalFunction(koAHSSDisplay, function(arg)
    local text, output;
    text := CallFuncList(koAHSSFormat,arg);
    # Use a separate unformatted stream so wide entries are not split by
    # GAP's pretty printer, and the caller's output settings stay intact.
    output := OutputTextUser();
    SetPrintFormattingStatus(output,false);
    PrintTo(output,text);
    CloseStream(output);
end);
