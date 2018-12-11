include("config.jl")
@printf("Running %s\n", basename(@__FILE__))

include("../lib/AnalyizeTimeseries.jl")
using .AnalyzeTimeseries
using PyPlot
using MultivariateStats
using Formatting

idx = (43, 21)


obs_SST = readModelVar("tos", (idx..., init_time:init_time+sim_len-1))


print(format("Simulation file: {}\n", sim_nc_filename))

ds = Dataset(hQ_nc_filename, "r")
h = nomissing(ds["h_mean"][:][idx..., :], NaN)
Q = nomissing(ds["Q_mean"][:][idx..., :], NaN)
close(ds)

ds = Dataset(sim_nc_filename, "r")
sim_SST = nomissing(ds["SST"][:][idx..., :], NaN)
close(ds)

println(size(sim_SST))
println(size(obs_SST))

#sim_SST = AnalyzeTimeseries.detrend(sim_SST)
#obs_SST = AnalyzeTimeseries.detrend(obs_SST)

c_sim = AnalyzeTimeseries.SpectralVariance(sim_SST)
c_obs = AnalyzeTimeseries.SpectralVariance(obs_SST)

# Let the 1-year peak be the reference intensity
norm_period = 12.0
norm_idx = floor(Int, length(sim_SST) / norm_period)
c_sim = c_sim / c_sim[norm_idx]
c_obs = c_obs / c_obs[norm_idx]

cutoff_period = 12.0 * [2.0, 50.0]
cutoff_coe = sort(length(sim_SST) ./ cutoff_period)
period_ticks = 12.0 * [.25, 0.5, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
freq_ticks   = length(sim_SST) ./ period_ticks
c_span = collect(Float64, 1:length(c_sim))

t = (collect(1:sim_len) .- 1) / 12.0

fig, ax = plt[:subplots](3, 1, figsize=(20, 10))

ax[1][:loglog](c_span, c_sim, "r-", label="sim_SST")
ax[1][:loglog](c_span, c_obs, "k--", label="obs_SST")

#ax[1][:set_xlim](cutoff_coe...)
ax[1][:set_xticks](freq_ticks)
ax[1][:set_xticklabels](["3m", "6m", "1y", "5y", "10y", "20y", "30y", "40y", "50y"])
#ax[1][:set_ylim]([0, 0.001])

#ax[1][:plot](t, obs_SST, "k--", label="obs_SST")
#ax[1][:plot](t, sim_SST, "r-", label="sim_SST")

#ax[1][:plot](c_span, c_obs, "k--", label="obs_SST")
#ax[1][:plot](c_span, c_sim, "r-", label="sim_SST")



ax[1][:set_ylim]([ 1e-1, 1e3])

#=
ax[2][:plot](1:12, h)
ax[3][:plot](1:12, Q)

ax[2][:invert_yaxis]()

ax[2][:set_title]("h")
ax[3][:set_title]("Q")


ax[2][:set_ylabel]("MLD [\$\\mathrm{m}\$]")
ax[3][:set_ylabel]("Q flux [\$\\mathrm{W} \\, \\mathrm{m}^{-2} \$]")
=#
fig[:suptitle](format("{}\n(lat, lon) = ({:.2f}, {:.2f})", sim_nc_filename, lat[idx[2]], lon[idx[1]]))




plt[:show]()

