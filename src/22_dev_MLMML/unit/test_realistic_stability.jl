include("../MLMML.jl")
include("../../lib/LinearRegression.jl")

using Printf
using Statistics: mean
using .MLMML
using Formatting

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

if(do_calculation)


zs = collect(Float64, range(0.0, stop=100.0 , length=21))
push!(zs, collect(Float64, range(100, step=50, stop=2000))[2:end]...)
push!(zs, collect(Float64, range(2000, step=100, stop=4000))[2:end]...)
zs *= -1.0

N = length(zs) - 1

Δb_init = 0.5 * 10.0 * MLMML.α * MLMML.g  
b_ML_init = 1.0
h_init = MLMML.h_min

b_slope = 2.0 / 4000.0 * MLMML.g * MLMML.α


PERIODS_SPINUP = 100
PERIODS_WANT = 200
PERIODS_TOTAL = PERIODS_SPINUP + PERIODS_WANT

PERIOD_CNT = 360
PERIOD_TIME = PERIOD_CNT * 86400.0
TOTAL_TIME = PERIODS_TOTAL * PERIOD_TIME





ω = 2π/360.0/86400.0
t = collect(Float64, range(0.0, step=86400.0, stop=TOTAL_TIME))[1:end-1]
Δt = t[2] - t[1]


J0 = 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p

J = J0 * sin.(ω*t)
E = J0/ω * (cos.(ω*t) .- 1.0)

U10 = zeros(Float64, length(t))
U10 .= 0.0 #.+ 2.0 * rand(length(U10)) 

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = b_slope,
    b_ML    = b_ML_init,
    h       = h_init,
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
    println("iteration = ", k)
    #println("oc.h = ", oc.h)
    info = MLMML.stepOceanColumn!(
        oc=oc,
        ua=U10[k],
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h)
    push!(b_rec, oc.b_ML)
    push!(Δb_rec, info[:Δb])
    push!(hb_rec, MLMML.getIntegratedBuoyancy(oc))
    push!(we_rec, ( info[:flag] == :we ) ? info[:val] : NaN)
    bs_rec[:, k+1] = oc.bs
end



new_rng = (PERIOD_CNT*PERIODS_SPINUP+1):length(t)
# Truncate
t = t[new_rng]
t .-= t[1]

bs_rec = bs_rec[:, new_rng]
J = J[new_rng]
h_rec = h_rec[new_rng]

using Statistics

# Doing Monthly Average
avg_t = zeros(Int(length(t)/30))
avg_bs_rec = zeros(size(bs_rec)[1], length(avg_t))
avg_h_rec  = zeros(length(avg_t))
avg_J = zeros(length(avg_t))
for i = 1:length(avg_t)
    avg_rng = (1+30*(i-1)):30*i
    avg_t[i] = mean(t[avg_rng])
    avg_h_rec[i] = mean(h_rec[avg_rng])
    avg_J[i] = mean(J[avg_rng])
    avg_bs_rec[:, i] = mean(bs_rec[:, avg_rng], dims=2)
end
avg_t .-= avg_t[1]

t=avg_t
bs_rec=avg_bs_rec
h_rec=avg_h_rec
J = avg_J


for i=1:length(oc.bs)
    b_timeseries = bs_rec[i, :]
    β = LinearRegression(t, b_timeseries)
    b_timeseries -= β[1] .+ β[2] * t
    #b_cyc_signal = repeat(mean( reshape( b_timeseries, PERIOD_CNT, :), dims=2)[:,1], outer=(PERIODS_WANT,))
    b_cyc_signal = repeat(mean( reshape( b_timeseries, 12, :), dims=2)[:,1], outer=(PERIODS_WANT,))
    bs_rec[i, :] = b_timeseries - b_cyc_signal
end

end

t_day = t / 86400.0
t_mon = t / 86400 / 30.0

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
ax1[:plot]([t_day[1], t_day[end]], [0, 0], "k--")

cmap = plt[:get_cmap]("jet")
clevs = (range(-1, stop=1, length=51) |> collect ) 
cbmapping = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec * 1e5, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = plt[:colorbar](cbmapping, cax=cax)

cb[:set_label]("Buoyancy anomaly [\$\\times\\,10^{-3}\\,\\mathrm{m} \\, \\mathrm{s}^{-2}\$]")

#ax2[:plot](t_day, - h_rec , "r--", linewidth=2, zorder=10)
ax2[:set_ylim]([-750, 0])

tlim = [t_day[1], t_day[end]]
ax1[:set_xlim](tlim)
ax2[:set_xlim](tlim)


ax1[:set_ylabel]("Insolation Flux [\$\\mathrm{m}^{2}\\, \\mathrm{s}^{-3}\$]")
ax2[:set_ylabel]("Z [m]")
ax2[:set_xlabel]("Time [year]")

CNT = PERIODS_WANT * PERIOD_CNT

xticks      = collect(0:PERIOD_CNT * 10:CNT)
xticklabels = [format("{:d}", xticks[i]/PERIOD_CNT) for i=1:length(xticks)]

ax1[:set_xticks](xticks)
ax2[:set_xticks](xticks)

ax1[:set_xticklabels](xticklabels)
ax2[:set_xticklabels](xticklabels)

using Formatting
fig[:suptitle](format("Buoyancy anomaly (annual cycle removed) with Δt = 1 day, Δz = {:.1f}", zs[1]-zs[2]))
plt[:show]()
