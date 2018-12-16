using Printf
@printf("Running %s\n", basename(@__FILE__))

include("../lib/AnalyizeTimeseries.jl")
include("../01_config/paths.jl")
using .AnalyzeTimeseries
using PyPlot
using MultivariateStats
using Formatting

using JLD
jld_fn = joinpath(data_path, "simulated_and_cpld_PDO.jld")
data = load(jld_fn)

obs_PDO = data["obs_PDO"]
sim_PDO = data["sim_PDO"]
nmons = length(sim_PDO)


c_sim = AnalyzeTimeseries.SpectralVariance(sim_PDO)
c_obs = AnalyzeTimeseries.SpectralVariance(obs_PDO)

# Let the 1-year peak be the reference intensity
norm_period = 12.0
norm_idx = floor(Int, length(sim_PDO) / norm_period)
#c_sim = c_sim / c_sim[norm_idx]
#c_obs = c_obs / c_obs[norm_idx]


cutoff_period = 12.0 * [2.0, 50.0]
cutoff_coe = sort(nmons ./ cutoff_period)
period_ticks = 12.0 * [.25, 0.5, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
freq_ticks   = nmons ./ period_ticks
c_span = collect(Float64, 1:length(c_sim))

t = (collect(1:nmons) .- 1) / 12.0

fig, ax = plt[:subplots](1, 1, figsize=(8,6))

ax[:loglog](c_span, c_sim, "r-", label="sim_PDO")
ax[:loglog](c_span, c_obs, "k--", label="obs_PDO")
ax[:legend]()
#ax[1][:set_xlim](cutoff_coe...)
ax[:set_xticks](freq_ticks)
ax[:set_xticklabels](["3m", "6m", "1y", "5y", "10y", "20y", "30y", "40y", "50y"])
#ax[1][:set_ylim]([0, 0.001])

#ax[1][:plot](t, obs_SST, "k--", label="obs_SST")
#ax[1][:plot](t, sim_SST, "r-", label="sim_SST")

#ax[1][:plot](c_span, c_obs, "k--", label="obs_SST")
#ax[1][:plot](c_span, c_sim, "r-", label="sim_SST")



ax[:set_ylim]([ 1e-1, 1e10])

#=
ax[2][:plot](1:12, h)
ax[3][:plot](1:12, Q)

ax[2][:invert_yaxis]()

ax[2][:set_title]("h")
ax[3][:set_title]("Q")


ax[2][:set_ylabel]("MLD [\$\\mathrm{m}\$]")
ax[3][:set_ylabel]("Q flux [\$\\mathrm{W} \\, \\mathrm{m}^{-2} \$]")
=#
#fig[:suptitle](format("{}\n(lat, lon) = ({:.2f}, {:.2f})", sim_nc_filename, lat[idx[2]], lon[idx[1]]))

plt[:show]()

