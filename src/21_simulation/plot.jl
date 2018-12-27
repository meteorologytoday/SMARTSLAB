using Printf
@printf("Running %s\n", basename(@__FILE__))

include("../lib/AnalyizeTimeseries.jl")
include("../01_config/paths.jl")
using .AnalyzeTimeseries
using PyPlot
using MultivariateStats
using Formatting
using Statistics: mean, std

using JLD
jld_fn = joinpath(data_path, "simulated_and_cpld_PDO.jld")
data = load(jld_fn)

obs_F_PDO = data["obs_F_PDO"]
obs_PDO = data["obs_PDO"]
sim_PDO = data["sim_PDO"]
nmons = length(sim_PDO)

c_sim = AnalyzeTimeseries.SpectralVariance(sim_PDO)
c_obs = AnalyzeTimeseries.SpectralVariance(obs_PDO)
c_obs_F = AnalyzeTimeseries.SpectralVariance(obs_F_PDO)

# Let the 1-year peak be the reference intensity
#norm_period = 12.0
#norm_idx = floor(Int, length(sim_PDO) / norm_period)
#c_sim = c_sim / c_sim[norm_idx]
#c_obs = c_obs / c_obs[norm_idx]


cutoff_period = 12.0 * [2.0, 50.0]
cutoff_coe = sort(nmons ./ cutoff_period)
period_ticks = 12.0 * [.25, 0.5, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0]
freq_ticks   = nmons ./ period_ticks
c_span = collect(Float64, 1:length(c_sim))

t = (collect(1:nmons) .- 1) / 12.0

# Figure 1: Frequency analysis
fig, ax = plt[:subplots](1, 1, figsize=(8,6))

ax[:loglog](c_span, c_sim, "r-", label="sim_PDO")
ax[:loglog](c_span, c_obs, "k-.", label="obs_PDO")
ax[:loglog](c_span, c_obs_F, "b--", label="obs_F_PDO")
ax[:legend]()

ax[:set_xticks](freq_ticks)
ax[:set_xticklabels](["3m", "6m", "1y", "5y", "10y", "20y", "30y", "40y", "50y"])

ax[:set_ylim]([ 1e-1, 1e10])

plt[:show]()


# Figure 2: PDO index timeseries

nmons  = length(obs_PDO)
nyrs   = Int(nmons / 12)
t_mon  = collect(Float64, 0:nmons-1) / 12.0
t_yr   = collect(Float64, 0:nyrs-1)


fig, ax = plt[:subplots](3, 1,sharex=true)

obs_PDO_annual_avg = mean(reshape(obs_PDO, 12, :), dims=1)[1,:]
obs_F_PDO_annual_avg = mean(reshape(obs_F_PDO, 12, :), dims=1)[1,:]
sim_PDO_annual_avg = mean(reshape(sim_PDO, 12, :), dims=1)[1,:]

ax[1][:set_title]("obs_PDO")
ax[1][:plot](t_mon, obs_PDO)
ax[1][:plot](t_yr,  obs_PDO_annual_avg)

ax[2][:set_title]("sim_PDO")
ax[2][:plot](t_mon, sim_PDO)
ax[2][:plot](t_yr,  sim_PDO_annual_avg)

ax[3][:set_title]("obs_F_PDO")
#ax[3][:plot](t_mon, obs_F_PDO)
ax[3][:plot](t_yr,  obs_F_PDO_annual_avg)


plt[:show]()
