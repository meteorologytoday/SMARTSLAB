using NetCDF
using Printf
@printf("Running %s\n", basename(@__FILE__))

@inline nansum(A::Array) = sum(A[isfinite.(A)])
@inline mod12(n) = mod(n-1, 12) + 1

function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%.2f ", A[i,j]) 

        end
        @printf("\n")
    end
end

dtype = Float64
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

rlons = convert(Array{dtype}, ncread(fn["omlmax"], "rlon"))
rlats = convert(Array{dtype}, ncread(fn["omlmax"], "rlat"))

#lat   = convert(Array{dtype}, ncread(fn["omlmax"], "lat"))
#=
stat_lat = collect(-90:90)
lat_mask = Array{Bool}(undef, length(stat_lat)-1, length(rlons), length(rlats))
for i = 1:size(lat_mask)[1]
    lat_mask[i, :, :] = (lat .>= stat_lat[i]) .& (lat .< stat_lat[i+1])
end
=#

rng = Colon()
#rng = 1:60

SST = convert(Array{dtype}, ncread(fn["tos"], "tos")[:,:,rng])
S   = convert(Array{dtype}, ncread(fn["S"], "S")[:,:,rng])
B   = convert(Array{dtype}, ncread(fn["B"], "B")[:,:,rng])

missing_value = convert(dtype, ncgetatt(fn["tos"], "tos", "missing_value"))
missing_places = (S .== ncgetatt(fn["S"], "S", "missing_value"))

SST[missing_places] .= NaN
S[missing_places] .= NaN
B[missing_places] .= NaN

θ = convert(Array{dtype}, SST * ρ * c_p)
F   = S + B

#println(eltype(SST))
#println(eltype(S))
#println(eltype(B))
#println(eltype(θ))

println("Data type: ", dtype)


mon_secs = 365.0 / 12.0 * 86400.0
nmons = size(SST)[3]

dt = 1.0 * mon_secs
dt2 = 2.0 * dt

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = Int(nmons / 12)
@printf("We have %02d years of data.\n", nyrs)

