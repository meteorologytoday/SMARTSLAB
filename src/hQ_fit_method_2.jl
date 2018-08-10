include("config.jl")
include("NetCDFHelper.jl")
using NetCDF



# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_star_dt = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) * ρ * c_p / (2.0 * dt)
N = size(TOT_F)[3]/12


ϕ = zeros(eltype(TOT_F), N, 2)
ϕ[:, 2] = 1.0

a = zeros(eltype(ϕ), length(rlons), length(rlats), 12)
b = copy(a)

for m = 1:12
    println("Doing month [$m]")
    for i = 1:length(rlons), j = 1:length(rlats)

        if mask[i,j]
            a[i,j,m] = NaN
            b[i,j,m] = NaN
            continue
        end

        ϕ[:, 1] = dT_star_dt[i,j,m:12:end]
        β = ϕ \ TOT_F[i,j,m:12:end]
        
        a[i,j,m] = β[1]
        b[i,j,m] = β[2]
    end
end

# Turn a, b into h, Q
h = 1.0 * a
Q = - b


h[isnan.(h)] = missing_value
Q[isnan.(Q)] = missing_value

time = collect(Float64, 1:12)

filename = "hQ_method_2.nc"
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])

# write h
nccreate(
    filename,
    "h",
    "rlon",
    "rlat",
    "time", time,
    atts= Dict(
        "long_name"=>"Mixed-layer Thickness",
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
