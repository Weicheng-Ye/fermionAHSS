# Standalone loader: Read("/absolute/path/to/fermionAHSS/load.g").
CallFuncList(function()
    local input, slash, directory, component, dependency;
    if IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then return; fi;
    if IsBoundGlobal("KOAHSS_PACKAGE_VERSION") then
        Error("fermionAHSS and koAHSS share globals; use a fresh GAP session");
    fi;
    for dependency in ["polycyclic", "hap", "json"] do
        if LoadPackage(dependency) = fail then
            Error("fermionAHSS requires GAP package ", dependency);
        fi;
    od;
    input := INPUT_FILENAME();
    slash := Positions(input, '/');
    if IsEmpty(slash) then
        directory := Directory(".");
    else
        slash := Last(slash);
        directory := Directory(input{[1..slash]});
    fi;
    if not IsBoundGlobal("koAHSSpages") then
        Read(Filename(directory, "gap/koahss.gd"));
    fi;
    Read(Filename(directory, "gap/pages.gi"));
    Read(Filename(directory, "gap/display.gi"));
    for component in ["cochains.gi", "natural_words.gi", "natural_bar.gi", "native_coherence.gi", "hap.gi", "phases.gi", "adem.gi", "universal_adem.gi", "integer_equations.gi", "secondary.gi", "natural_secondary.gi", "secondary_corrected.gi", "secondary_squarezero.gi", "normalization.gi", "secondary_coherence.gi", "defining_systems.gi", "tertiary_correction.gi", "tertiary_candidates.gi", "tertiary_equations.gi", "natural_tertiary.gi", "group_api.gi"] do
        Read(Filename(directory, Concatenation("gap/", component)));
    od;
    BindGlobal("KOAHSS_PACKAGE_VERSION", "0.1.0");
    BindGlobal("FERMION_AHSS_PACKAGE_VERSION", "0.1.0");
end, []);
