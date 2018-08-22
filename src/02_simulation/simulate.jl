include("config_and_data.jl")
include("single_point.jl")

@printf("Running %s\n", basename(@__FILE__))

#println(rlats)
#println(rlons)

rlat_i = 130
rlon_i = 110

real_lat = lats[rlon_i, rlat_i]
real_lon = lons[rlon_i, rlat_i]

@printf("Rotated Longitude: %5.2f, Rotated Latitude: %5.2f\n", rlons[rlon_i], rlats[rlat_i])
@printf("        Longitude: %5.2f,         Latitude: %5.2f\n", lons[rlon_i, rlat_i], lats[rlon_i, rlat_i])

Qh_fn      = joinpath(data_path, "case3_hQ.jl.nc")
@printf("Will read file: %s\n", Qh_fn)


Q = convert(Array{Float64, 1}, ncread(Qh_fn,      "Q"    )[rlon_i, rlat_i, :])
h = convert(Array{Float64, 1}, ncread(Qh_fn,      "h"    )[rlon_i, rlat_i, :])

F = convert(Array{Float64, 1}, TOT_F[rlon_i, rlat_i, :])

println(Q)
println(h)


init_T_star  = T_star[rlon_i, rlat_i, 1] 
ret_len = length(F)
h_min   = 10.0

# detect h_min
for i = 1:length(h)
    h[i] = h[i] < h_min ? h_min : h[i]
end

# construct dhdt
dhdt    = copy(h)
dhdt[2:end-1] = h[3:end] - h[1:end-2]
dhdt[1]       = h[2] - h[end]
dhdt[end]     = h[1] - h[end-1]
dhdt /= dt * 2

@printf("init_T_star  : %.2f\n", init_T_star)
@printf("ret_len : %d\n", ret_len)

@printf("Doing simulation... \n")
T = simulate_single_point(
        init_T_star,
        dt,
        h,
        dhdt,
        Q,
        F,
        ret_len
) / ρ / c_p

@printf("done.\n")

using DataFrames
using CSV



df = DataFrame()
df[:SST_simulated] = T
df[:SST_real]      = SST[rlon_i, rlat_i, :]

filename = @sprintf("%s_SST.csv", basename(@__FILE__))
filename = joinpath(data_path, filename)

@printf("Gonna print data to %s\n", filename)
CSV.write(filename, df)

