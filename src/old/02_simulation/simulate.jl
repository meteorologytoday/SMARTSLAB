include("config_and_data.jl")
include("single_point.jl")

@printf("Running %s\n", basename(@__FILE__))

#println(rlats)
#println(rlons)

rlat_i = parse(Int, ARGS[1])
rlon_i = parse(Int, ARGS[2])

desc = nothing
if length(ARGS) >= 3
    desc = ARGS[3]
end

@printf("desc: %s\n", desc)

real_lat = lats[rlon_i, rlat_i]
real_lon = lons[rlon_i, rlat_i]

@printf("Rotated Longitude: %5.2f, Rotated Latitude: %5.2f\n", rlons[rlon_i], rlats[rlat_i])
@printf("        Longitude: %5.2f,         Latitude: %5.2f\n", lons[rlon_i, rlat_i], lats[rlon_i, rlat_i])

Qh_fn      = joinpath(data_path, "case3_hQ.jl.nc")
@printf("Will read file: %s\n", Qh_fn)


Q = convert(Array{Float64, 1}, ncread(Qh_fn,      "Q"    )[rlon_i, rlat_i, :])
h = convert(Array{Float64, 1}, ncread(Qh_fn,      "h"    )[rlon_i, rlat_i, :])

F_fn      = joinpath(data_path, "dT_star_dt-TOT_F.nc")

@printf("Will read file: %s\n", F_fn)
F = convert(Array{Float64, 1}, TOT_F[rlon_i, rlat_i, :])
#F = convert(Array{Float64, 1}, ncread(F_fn, "TOT_F")[rlon_i, rlat_i, :])

init_T_star  = T_star[rlon_i, rlat_i, 1] 
ret_len = 30 * 12 #length(1200)
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
T_RK4, T_EULER, used_F, used_Q, used_h, used_dhdt = simulate_single_point(
        init_T_star,
        dt,
        F,
        Q,
        h,
        dhdt,
        ret_len
)

T_RK4   /= ρ * c_p
T_EULER /= ρ * c_p

DEBUG_T_EULER_DIF = T_RK4 * 0.0
DEBUG_T_EULER_DIF[1:end-1] = T_EULER[2:end] - T_EULER[1:end-1]

DEBUG_T_RK4_DIF = T_RK4 * 0.0
DEBUG_T_RK4_DIF[1:end-1] = T_RK4[2:end] - T_RK4[1:end-1]



@printf("done.\n")

@printf("Making Plots...\n")

dat = Dict(
    :SST_RK4           => T_RK4,
    :SST_EULER         => T_EULER,
    :SST_REAL          => SST[rlon_i, rlat_i, 1:ret_len],
    :F                 => used_F,
    :Q                 => used_Q,
    :h                 => used_h,
    :dhdt              => used_dhdt,
    Symbol("-dhdt/h")  => - used_dhdt ./ used_h,
    Symbol("(F+Q)/h")  => (used_F .+ used_Q) ./ used_h,
    Symbol("-dhdt/h * SST_EULER * ρ * c_p")  => - used_dhdt ./ used_h .* T_EULER .* (ρ * c_p),
    Symbol("TOTAL / ρ / c_p")  => ( used_F .+ used_Q - T_EULER .* (ρ * c_p) .* used_dhdt) ./ used_h / (ρ * c_p),
    Symbol("TOTAL / ρ / c_p * dt")  => ( used_F .+ used_Q - T_EULER .* (ρ * c_p) .* used_dhdt) ./ used_h * dt / (ρ * c_p),
    Symbol("DEBUG_T_EULER_DIF")  => DEBUG_T_EULER_DIF, 
    Symbol("DEBUG_T_RK4_DIF")  => DEBUG_T_RK4_DIF 
)


# Print CSV file
using DataFrames
using CSV

df = DataFrame(dat)

filename = @sprintf("%s-[%.2f][%.2f].csv", basename(@__FILE__), real_lat, real_lon)
filename = joinpath(data_path, filename)

@printf("Gonna print data to %s\n", filename)
CSV.write(filename, df)





using PyPlot


t = collect(1:length(dat[:SST_REAL]))
mean_T = round(mean(dat[:SST_REAL]) / 2) * 2.0


draw_months = ret_len

fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(20, 10))

title = @sprintf("\$( \\phi, \\lambda ) = ( %.2f, %.2f )\$", real_lat, real_lon)
if desc != nothing
    title = @sprintf("%s\n", desc) * title
end

fig[:suptitle](title)
ax2_twin = ax[2][:twinx]()



ax[1][:set_xlim](1, draw_months)
ax[1][:set_ylim](mean_T - 10, mean_T + 10)

ax[3][:set_xticks](collect(1:12:draw_months))
ax[3][:set_xticklabels](collect(1:floor(Int, draw_months/12.0)))
ax[3][:set_xlabel]("Year")

ax[1][:set_ylabel]("Year")
ax[2][:set_ylabel](L"\mathrm{W} \, \mathrm{m}^{-2}")
ax2_twin[:set_ylabel](L"\mathrm{m}")
ax[3][:set_ylabel](L"\mathrm{K} \, \mathrm{s}^{-1}")




ax[1][:set_title]("SST")
ax[1][:plot](t, dat[:SST_REAL], color="r", linewidth=2.0, linestyle="-", label="Real")
ax[1][:plot](t, dat[:SST_RK4], color="g", linewidth=2.0, dashes=(7,3), label="RK4")
ax[1][:plot](t, dat[:SST_EULER], color="b", linewidth=2.0, dashes=(7,3,3,3), label="EULER")

ax[1][:legend]()

ax[2][:set_title]("F, Q, h")
ax[2][:plot](t, dat[:F], color="k", linewidth=2.0, linestyle="-", label="F")
ax[2][:plot](t, dat[:Q] ./ 100.0, color="k", linewidth=2.0, linestyle="--", label="Q/100")


ax2_twin[:plot](t, dat[:h], color="k", linewidth=2.0, linestyle="-.", label="h")

ax[2][:legend]()
ax2_twin[:legend]()


ax[3][:set_title]("F/h, Q/h, - T/h * dh/dt")
ax[3][:plot](t, dat[:F] ./ dat[:h], color="k", linewidth=2.0, linestyle="-", label="F/h")
ax[3][:plot](t, dat[:Q] ./ dat[:h], color="k", linewidth=2.0, linestyle="--", label="Q/h")
ax[3][:plot](t, - dat[:SST_REAL] ./ dat[:h] .* dat[:dhdt] .* (ρ*c_p), color="r", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (Real)")
ax[3][:plot](t, - dat[:SST_RK4] ./ dat[:h] .* dat[:dhdt] .* (ρ*c_p), color="g", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (RK4)")
ax[3][:plot](t, - dat[:SST_EULER] ./ dat[:h] .* dat[:dhdt] .* (ρ*c_p), color="b", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (EULER)")
ax[3][:legend]()


filename = @sprintf("%s-[%.2f][%.2f].png", basename(@__FILE__), real_lat, real_lon)
fig[:savefig](joinpath(img_path, filename), dpi=200)

#plt[:show]()

#=
using DataFrames
using CSV


df = DataFrame()
df[:SST_RK4]   = T_RK4
df[:SST_EULER] = T_EULER
df[:SST_REAL]  = SST[rlon_i, rlat_i, 1:ret_len]
df[:F]         = used_F
df[:Q]         = used_Q
df[:h]         = used_h
df[:dhdt]      = used_dhdt

filename = @sprintf("%s_SST.csv", basename(@__FILE__))
filename = joinpath(data_path, filename)

@printf("Gonna print data to %s\n", filename)
CSV.write(filename, df)

include("config.jl")

@printf("Running %s\n", basename(@__FILE__))

using PyPlot
using DataFrames
using CSV

df = CSV.read(joinpath(data_path, "simulate.jl_SST.csv"))


fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(20, 10))

ax[1][:set_xlim](1,240)
ax[1][:set_ylim](285,310)

t = collect(1:length(df[:SST_REAL]))

ax[1][:set_title]("SST budget")
ax[1][:plot](t, df[:SST_REAL], color="r", linewidth=2.0, linestyle="-", label="Real")
ax[1][:plot](t, df[:SST_RK4], color="g", linewidth=2.0, dashes=(7,3), label="RK4")
ax[1][:plot](t, df[:SST_EULER], color="b", linewidth=2.0, dashes=(7,3,3,3), label="EULER")

ax[1][:legend]()

ax[2][:set_title]("F, Q, h")
ax[2][:plot](t, df[:F], color="k", linewidth=2.0, linestyle="-", label="F")
ax[2][:plot](t, df[:Q] ./ 100.0, color="k", linewidth=2.0, linestyle="--", label="Q/100")

ax2_twin = ax[2][:twinx]()
ax2_twin[:plot](t, df[:h], color="k", linewidth=2.0, linestyle="-.", label="h")

ax[2][:legend]()


ax[3][:set_title]("F/h, Q/h, - T/h * dh/dt")
ax[3][:plot](t, df[:F] ./ df[:h], color="k", linewidth=2.0, linestyle="-", label="F/h")
ax[3][:plot](t, df[:Q] ./ df[:h], color="k", linewidth=2.0, linestyle="--", label="Q/h")
ax[3][:plot](t, - df[:SST_REAL] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="r", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (Real)")
ax[3][:plot](t, - df[:SST_RK4] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="g", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (RK4)")
ax[3][:plot](t, - df[:SST_EULER] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="b", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (EULER)")
ax[3][:legend]()


fig[:savefig]("")
=#

