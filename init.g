# GAP loads package declarations before implementations.
if not IsBoundGlobal("FERMION_AHSS_PACKAGE_VERSION") and not IsBoundGlobal("koAHSSpages") then
    ReadPackage("fermionAHSS", "gap/api.gd");
fi;
