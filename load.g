# Standalone loader: Read("/absolute/path/to/fermionAHSS_new/load.g").
CallFuncList(function()
    local input, slash, directory, component, dependency;
    if IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then return; fi;
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
        directory := Directory(input{[1..Last(slash)]});
    fi;
    if not IsBoundGlobal("koAHSSpages") then
        Read(Filename(directory, "gap/api.gd"));
    fi;
    for component in ["pages.gi", "display.gi", "cochain_space.gi", "resolution_memo.gi",
        "bar_transport.gi", "binary_tensors.gi", "native_coherence.gi", "hap_space.gi",
        "linear_algebra.gi", "primary_transport.gi", "cell_transport.gi", "worker_model.gi",
        "page_operations.gi", "api.gi", "result_cache.gi",
        "extension.gi", "flat_lifts.gi", "gauge_compare.gi", "relations.gi", "light_rows.gi"] do
        Read(Filename(directory, Concatenation("gap/", component)));
    od;
    BindGlobal("FERMION_AHSS_PACKAGE_VERSION", "1.0.0");
end, []);
