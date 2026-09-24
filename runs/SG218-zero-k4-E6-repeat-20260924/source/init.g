# GAP loads package declarations before implementations.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") then
    if IsBoundGlobal("KOAHSS_PACKAGE_VERSION") then
        Error("fermionAHSS and koAHSS share globals; use a fresh GAP session");
    fi;
    if not IsBoundGlobal("koAHSSpages") then
        ReadPackage("fermionAHSS", "gap/koahss.gd");
    fi;
fi;
