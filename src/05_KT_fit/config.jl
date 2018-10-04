using NetCDF

@printf("Running %s\n", basename(@__FILE__))

ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K
mon_secs = 365.0 / 12.0 * 86400.0
dt  = 1.0 * mon_secs

data_path = joinpath(dirname(@__FILE__), "..", "..", "data")
img_path = joinpath(dirname(@__FILE__), "..", "..", "img")

fn = Dict()

for varname in ["omlmax", "S", "B", "tos"]

    fn[varname] = joinpath(data_path, @sprintf("SMART_%s_Omon_GFDL-ESM2G_historical_r1i1p1_186101-200512.nc", varname))

end

rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")

SST   = ncread(fn, "tos")[:,:,:]
tp = eltype(TOT_F)

missing_value = ncgetatt(fn, "tos", "missing_value")

TOT_F[TOT_F .== ncgetatt(fn, "total_downward_heat_flux", "missing_value")] = NaN
SST[SST .== ncgetatt(fn, "tos", "missing_value")] = NaN

spatial_mask = isnan.(SST[:,:,1])
spatial_temporal_mask = isnan.(SST)

T_star = SST * ρ * c_p

mon_secs = 365.0 / 12.0 * 86400.0
nmons = length(ncread(fn, "time"))

dt = 1.0 * mon_secs

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = nmons / 12
@printf("We have %02d years of data.\n", nyrs)


dtype = eltype(T_star)
