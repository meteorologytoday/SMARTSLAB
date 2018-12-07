include("config.jl")
@printf("Running %s\n", basename(@__FILE__))
using PyPlot
using MultivariateStats

idx = (38, 29)


obs_SST = readModelVar("tos", (idx..., init_time:init_time+sim_len-1))

ds = Dataset(hQ_nc_filename, "r")
h = nomissing(ds["h_mean"][:][idx..., :], NaN)
Q = nomissing(ds["Q_mean"][:][idx..., :], NaN)
close(ds)

ds = Dataset(sim_nc_filename, "r")
sim_SST = nomissing(ds["SST"][:][idx..., :], NaN)
close(ds)

println(size(sim_SST))
println(size(obs_SST))


t = (collect(1:sim_len) .- 1) / 12.0

fig, ax = plt[:subplots](3, 1, figsize=(20, 10))

#ax[1][:set_xlim](1,240)
#ax[1][:set_ylim](285,310)

ax[1][:plot](t, obs_SST, "k--", label="obs_SST")
ax[1][:plot](t, sim_SST, "r-", label="sim_SST")

ax[2][:plot](1:12, h)
ax[3][:plot](1:12, Q)

plt[:show]()

