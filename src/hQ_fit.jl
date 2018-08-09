include("config.jl")
include("NetCDFHelper.jl")
using NetCDF



# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_star = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) * ρ * c_p
N = size(TOT_F)[3]/12


ϕ = zeros(eltype(TOT_F), N, 2)
ϕ[:, 2] = 1.0

h = zeros(eltype(ϕ), length(rlons), length(rlats), 12)
Q = copy(h)

for m = 1:12
    println("Doing month [$m]")
    for i = 1:length(rlons), j = 1:length(rlats)

        if isnan(TOT_F[i,j,1])
            h[i,j,m] = NaN
            Q[i,j,m] = NaN
            continue
        end

        ϕ[:, 1] = TOT_F[i,j,m:12:end]
        β = ϕ \ dT_star[i,j,m:12:end]
        
        h[i,j,m] = β[1]
        Q[i,j,m] = β[2]
    end
end

# Turn a, b into h, Q
Q = Q ./ h
h = 2.0 * dt ./ h

h[isnan.(h)] = missing_value
Q[isnan.(Q)] = missing_value

time = collect(Float64, 1:12)

filename = "hQ.nc"
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])



# write h
nccreate(
    filename,
    "h",
    "rlon",
    "rlat",
    "time", time,
    atts= Dict(
        "long_name"=>"Constant Mixed-layer Thickness",
        "units"=>"m",
        "missing_value" => missing_value
    )
)
ncwrite(h, filename, "h")

# write Q
nccreate(
    filename,
    "Q",
    "rlon",
    "rlat",
    "time", time,
    atts= Dict(
        "long_name"=>"Q-flux",
        "units"=>"J / m^2",
        "missing_value" => missing_value
    )
)
ncwrite(Q, filename, "Q")

ncclose(filename)
