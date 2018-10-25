@printf("# Running %s\n", basename(@__FILE__))

models = Dict(
    "NCAR" => joinpath("NCAR_CESM1-WACCM", "%s_Omon_CESM1-WACCM_piControl_r1i1p1_009601-029512.nc"),
    "GFDL" => joinpath("GFDL_GFDL-ESM2G", "SMART_%s_Omon_GFDL-ESM2G_historical_r1i1p1_186101-200512.nc"),
)


