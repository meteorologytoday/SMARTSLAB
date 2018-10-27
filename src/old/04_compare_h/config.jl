@printf("Running %s\n", basename(@__FILE__))

ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K
mon_secs = 365.0 / 12.0 * 86400.0
dt  = 1.0 * mon_secs

data_path = joinpath(dirname(@__FILE__), "..", "..", "data")
img_path = joinpath(dirname(@__FILE__), "..", "..", "img")

fn = Dict()

for varname in ["omlmax"]

    fn[varname] = joinpath(data_path, @sprintf("SMART_%s_Omon_GFDL-ESM2G_historical_r1i1p1_186101-200512.nc", varname))

end




