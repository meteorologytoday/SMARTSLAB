using NetCDF
using Printf
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

rlons    = ncread(fn["omlmax"], "rlon")
rlats    = ncread(fn["omlmax"], "rlat")


rng = Colon()
rng = 1:60


SST   = ncread(fn["tos"], "tos")[:,:,rng]
S   = ncread(fn["S"], "S")[:,:,rng]
B   = ncread(fn["B"], "B")[:,:,rng]

dtype = eltype(SST)

missing_value = ncgetatt(fn["tos"], "tos", "missing_value")
missing_places = (S .== ncgetatt(fn["S"], "S", "missing_value"))


SST[missing_places] .= NaN
S[missing_places] .= NaN
B[missing_places] .= NaN

T_star = SST * ρ * c_p

mon_secs = 365.0 / 12.0 * 86400.0
nmons = size(SST)[3]

dt = 1.0 * mon_secs

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = Int(nmons / 12)
@printf("We have %02d years of data.\n", nyrs)

