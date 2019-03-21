include("../MLMML.jl")
include("../../../lib/LinearRegression.jl")

using Printf
using Statistics: mean
using .MLMML
using Formatting
using StatsBase

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

if(do_calculation)


zs = [0.0, 2000.0]
zs *= -1.0

N = length(zs) - 1

Δb_0 = 20.0 * MLMML.α * MLMML.g  
b_ML_0 = 1.0
h_ML_0 = MLMML.h_ML_min

SECS_PER_DAY = 86400.0
DAYS_PER_MON = 30
MONS_PER_YEAR= 12
DAYS_PER_YEAR = DAYS_PER_MON * MONS_PER_YEAR
SECS_PER_YEAR = DAYS_PER_YEAR * SECS_PER_DAY

SPINUP_YEARS = 10
YEARS_WANTED = 50
TOTAL_YEARS = SPINUP_YEARS + YEARS_WANTED

TOTAL_DAYS = TOTAL_YEARS * DAYS_PER_YEAR
TOTAL_SECS = TOTAL_DAYS * SECS_PER_DAY

SPINUP_DAYS = SPINUP_YEARS * DAYS_PER_YEAR
DAYS_WANTED = YEARS_WANTED * DAYS_PER_YEAR

ω = 2π/360.0/86400.0
t_sim = collect(Float64, range(0.0, step=SECS_PER_DAY, stop=TOTAL_SECS))[1:end-1]
Δt = t_sim[2] - t_sim[1]

t = t_sim[SPINUP_DAYS+1:end]
t .-= t[1]


J0 = 125.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p

J = J0 * cos.(ω*t_sim)

U10 = zeros(Float64, length(t_sim))
U10 .= 8.0 #.+ 2.0 * cos.(ω*t_sim)

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = 0.0,
    b_ML    = b_ML_0,
    h_ML    = h_ML_0,
    Δb      = Δb_0,
    K       = 0.0,
)

mon_bs_detrend = zeros(Float64, length(oc.bs), MONS_PER_YEAR)
mon_bs = zeros(Float64, length(oc.bs), MONS_PER_YEAR)
trends = zeros(Float64, length(oc.bs))
means = zeros(Float64, length(oc.bs))

h_rec = zeros(Float64, DAYS_WANTED)
b_rec = copy(h_rec)
hb_rec = copy(h_rec)
J_rec = copy(h_rec)
we_rec = copy(h_rec)
convadj_rec = copy(h_rec)
Δb_rec = copy(h_rec)
bs_rec = zeros(Float64, length(oc.bs), DAYS_WANTED)

for k = 1:length(t_sim)
    println("iteration = ", k, "/", length(t_sim))
    if k != 1
        info = MLMML.stepOceanColumn!(
            oc=oc,
            ua=U10[k],
            B0=0.0,
            J0=J[k],
            Δt=Δt,
        )
    end

    if k > SPINUP_DAYS
        h_rec[i]       = oc.h_ML
        b_rec[i]       = oc.b_ML
        Δb_rec[i]      = info[:Δb]
        hb_rec[i]      = MLMML.getIntegratedBuoyancy(oc)
        we_rec[i]      = ( info[:flag] == :we ) ? info[:val] : NaN
        convadj_rec[i] = info[:convective_adjustment]
        bs_rec[:, i]   = oc.bs
    end
end

using Statistics

# Doing Monthly Average
avg_t = zeros(Int(length(t)/DAYS_PER_MON))
avg_bs_rec = zeros(size(bs_rec)[1], length(avg_t))
avg_h_rec  = zeros(length(avg_t))
avg_J_rec = zeros(length(avg_t))
avg_convadj_rec = zeros(length(avg_t))
for i = 1:length(avg_t)
    avg_rng = (1+DAYS_PER_MON*(i-1)):DAYS_PER_MON*i
    avg_t[i] = mean(t[avg_rng])
    avg_h_rec[i] = mean(h_rec[avg_rng])
    avg_J_rec[i] = mean(J[avg_rng])
    avg_bs_rec[:, i] = mean(bs_rec[:, avg_rng], dims=2)
    avg_convadj_rec[i] = (sum(convadj_rec[avg_rng]) != 0) ? 1.0 : 0.0
end
avg_t .-= avg_t[1]

t=avg_t
bs_rec=avg_bs_rec
h_rec=avg_h_rec
J_rec = avg_J_rec
convadj_rec = avg_convadj_rec

# Remove seasonal cycle
for i=1:length(oc.bs)
    b_timeseries = bs_rec[i, :]
    mon_bs[i, :] = mean( reshape( b_timeseries, MONS_PER_YEAR, :), dims=2)[:,1]

    # Detrend
    β = LinearRegression(t, b_timeseries)
    b_timeseries -= β[1] .+ β[2] * t
    mon_bs_detrend[i, :] = mean( reshape( b_timeseries, MONS_PER_YEAR, :), dims=2)[:,1]
    
    means[i]  = β[1]
    trends[i] = β[2]

    
    b_cyc_signal = repeat(mon_bs_detrend[i, :], outer=(YEARS_WANTED,))
    bs_rec[i, :] = b_timeseries - b_cyc_signal
end

end

t_day = t / 86400.0
t_mon = t / 86400 / 30.0
mid_zs = (zs[1:end-1] + zs[2:end]) / 2.0
b2T = 1.0 / (MLMML.α * MLMML.g)

#SST_rec = mean( reshape(bs_rec[1,:] * b2T, 12*5, :), dims=1)[1, :]
#t_SST_rec = mean( reshape(t_day, 12*5, :), dims=1)[1, :]

SST_rec = bs_rec[1, :] * b2T
t_SST_rec = t_day

T_rec = bs_rec * b2T
mon_SST = mon_bs * b2T
mon_SST_detrend = mon_bs_detrend * b2T
T_means = means * b2T
T_trends = trends * b2T



auto_SST = autocor(SST_rec, collect(0:120); demean=true) 
auto_SST_yr = collect(Float64, 0:length(auto_SST)-1) / 12.0

plt[:figure]()
plt[:plot](auto_SST_yr, auto_SST)







# Hovmoller diagram
fig, (ax1, ax2, ax3) = plt[:subplots](3, 1, figsize=(12, 6), sharex=true)

# SST record
ax1[:plot]([t_day[1], t_day[end]], [0, 0], "k--")
ax1[:plot](t_SST_rec, SST_rec, label="SST")

ax2[:plot](t_day, T_rec[1, :])

ax3[:plot](t_day, -h_rec)

convadj_rec[convadj_rec .== 0.0] .= NaN
ax3[:set_ylim]([-1500, 0])





tlim = [t_day[1], t_day[end]]
ax1[:set_xlim](tlim)

ax1[:set_ylabel]("SST\n[\$\\mathrm{K}\$]")
ax2[:set_ylabel]("Deep Ocean Temperature\n[\$\\mathrm{K}\$]")
ax3[:set_ylabel]("Z [m]")
ax3[:set_xlabel]("Time [year]")

ticks      = collect(0:DAYS_PER_YEAR * 50:DAYS_WANTED)
ticklabels = [format("{:d}", ticks[i]/DAYS_PER_YEAR) for i=1:length(ticks)]

ax1[:set_xticks](ticks)
ax2[:set_xticks](ticks)

ax1[:set_xticklabels](ticklabels)
ax2[:set_xticklabels](ticklabels)

using Formatting
fig[:suptitle](
    format(
        "Temperature anomaly (annual cycle removed)\nΔt = 1 day, spin up time = {:d} years", SPINUP_YEARS)
    )
plt[:show]()


# Monthly structure
fig, ax = plt[:subplots](1, 4, figsize=(20,6), sharey=true)

fig[:suptitle]("Vertical profile of monthly mean buoyancy")

for i = 1:size(mon_bs)[2]

    offset =  i

    ax[1][:plot](mon_SST[:, i], mid_zs, label="$i")
    ax[1][:text](mon_SST[1, i], zs[1]+offset, "$i", va="bottom", ha="center")

    ax[2][:plot](mon_SST_detrend[:, i], mid_zs, label="$i")
    ax[2][:text](mon_SST_detrend[1, i], zs[1]+offset, "$i", va="bottom", ha="center")
end

ax[1][:set_title]("Original")
ax[2][:set_title]("After Detrend")

ax[1][:set_ylim]([-200, 20])
ax[1][:legend]()
ax[2][:legend]()

ax[1][:set_xlabel]("Temperature anomaly [\$\\mathrm{K}\$]")
ax[2][:set_xlabel]("Temperature anomaly [\$\\mathrm{K}\$]")


ax[3][:plot]([0, 0], [-200, 20], "--", color="#888888")
ax[3][:plot](T_trends * SECS_PER_YEAR * 1e4, mid_zs, "r-", label="trend")
ax[3][:set_title]("Trend for each layer")

ax[4][:plot](T_means, mid_zs, "r-", label="mean")
ax[4][:set_title]("Mean for each layer")

ax[3][:set_xlabel]("Trend of Temperature [\$\\times\\,10^{-4}\\,\\mathrm{K}\\,\\mathrm{yr}^{-1} \$]")
ax[4][:set_xlabel]("Mean Temperature [\$\\mathrm{K}\$]")
plt[:show]()
