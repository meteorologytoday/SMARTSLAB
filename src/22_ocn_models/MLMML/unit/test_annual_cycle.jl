include("../MLMML.jl")
include("../../../lib/LinearRegression.jl")
using Printf
using Statistics: mean
using .MLMML
using Formatting

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

D  = 2000.0
N  = 2001
zs = collect(Float64, range(0.0, stop=-D, length=N))

Δb_init = 0.5 * 10.0 * MLMML.α * MLMML.g  
b_ML_init = 1.0
h_init = 50.0

b_slope = 10.0 / D * MLMML.g * MLMML.α


PERIOD_N = 30

PERIOD_CNT = 360
PERIOD_TIME = PERIOD_CNT * 86400.0
TOTAL_TIME = PERIOD_N * PERIOD_TIME

SPINUP_TIME = 5 * PERIOD_TIME


ω = 2π/360.0/86400.0
t = collect(Float64, range(0.0, step=86400.0, stop=TOTAL_TIME))[1:end-1]
Δt = t[2] - t[1]
t_day = t / 86400.0

J0 = 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p

J = J0 * sin.(ω*t)
E = J0/ω * (cos.(ω*t) .- 1.0)

U10 = zeros(Float64, length(t))
U10 .= 0.0 #.+ 2.0 * rand(length(U10)) 

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = b_slope,
    b_ML    = b_ML_init,
    h_ML    = h_init,
    Δb      = Δb_init,
    K       = 1e-5
)

#if ! MLMML.checkDiffusionStability(oc, Δt=Δt)
#    throw(ErrorException("Stability criteria does not fulfill."))
#end



h_rec = [h_init]
b_rec = [b_ML_init]
hb_rec = [MLMML.getIntegratedBuoyancy(oc)]
we_rec = []
Δb_rec = [oc.b_ML - oc.bs[oc.FLDO]]
bs_rec = zeros(Float64, length(oc.bs), length(t))


bs_rec[:, 1] = oc.bs
for k = 1:length(t)-1
    println("iteration = ", k, "/", length(t)-1)
    #println("oc.h = ", oc.h)
    info = MLMML.stepOceanColumn!(
        oc=oc,
        ua=U10[k],
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h_ML)
    push!(b_rec, oc.b_ML)
    push!(Δb_rec, info[:Δb])
    push!(hb_rec, MLMML.getIntegratedBuoyancy(oc))
    push!(we_rec, ( info[:flag] == :we ) ? info[:val] : NaN)
    bs_rec[:, k+1] = oc.bs
end

for i=1:length(oc.bs)
    b_timeseries = bs_rec[i, :]
    β = LinearRegression(t, b_timeseries)
    b_timeseries -= β[1] .+ β[2] * t
    b_cyc_signal = repeat(mean( reshape( b_timeseries, PERIOD_CNT, :), dims=2)[:,1], outer=(PERIOD_N,))
    bs_rec[i, :] = b_timeseries - b_cyc_signal
    

end

#=
# Line plot figure
plt[:figure]()
for i = 1:5:size(bs_rec)[2]
    plt[:plot](bs_rec[:, i], (oc.zs[1:end-1] + oc.zs[2:end]) / 2.0, "k-")
end

=#
#=
# Diagnose Plot
fig, ax = plt[:subplots](5, 1, sharex=true)

ax[1][:plot](t_day, - h_rec, label="h")
ax[2][:plot](t_day, hb_rec, "k-", label="hb")
ax[2][:plot](t_day, hb_rec[1] .+ E, "r--", label="hb theory")
ax[3][:plot](t_day[1:end-1], we_rec, label="we")

ax[4][:plot](t_day, Δb_rec, label="Δb")
#ax[4][:plot](t_day, , label="Δb theory")
ax[5][:plot](t_day, J, label="J")

for a in ax
    a[:legend]()
end
=#

# Hovmoller diagram
gs0 = GS.GridSpec(1, 2, width_ratios=[100,5])
gs_l = GS.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], height_ratios=[1, 4])

fig = plt[:figure](figsize=(12, 6))
ax1 = plt[:subplot](gs_l[1])
ax2 = plt[:subplot](gs_l[2])
cax = plt[:subplot](gs0[2])


ax1[:plot](t_day, J, label="J")
ax1[:plot]([t[1], t[end]]/86400.0, [0, 0], "k--")

cmap = plt[:get_cmap]("jet")
clevs = (range(-1, stop=1, length=51) |> collect ) * 20
cbmapping = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec * 1e4, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, cax=cax)

cb[:set_label]("Buoyancy anomaly [\$\\times\\,10^{-3}\\,\\mathrm{m} \\, \\mathrm{s}^{-2}\$]")

ax2[:plot](t_day, - h_rec , "r--", linewidth=2, zorder=10)
ax2[:set_ylim]([-D, 0])

tlim = [t_day[1], t_day[end]]
ax1[:set_xlim](tlim)
ax2[:set_xlim](tlim)


ax1[:set_ylabel]("Insolation Flux [\$\\mathrm{m}^{2}\\, \\mathrm{s}^{-3}\$]")
ax2[:set_ylabel]("Z [m]")
ax2[:set_xlabel]("Time [year]")

xticks      = collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0
xticklabels = [format("{:d}", i-1) for i=1:length(xticks)]

ax1[:set_xticks](collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0)  
ax2[:set_xticks](collect(range(0.0, step=PERIOD_TIME, stop=TOTAL_TIME))/86400.0)

ax1[:set_xticklabels](xticklabels)
ax2[:set_xticklabels](xticklabels)

using Formatting
fig[:suptitle](format("Buoyancy anomaly (annual cycle removed) with Δt = 1 day, Δz = {:.1f}", zs[1]-zs[2]))
plt[:show]()
