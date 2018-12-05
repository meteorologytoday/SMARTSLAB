
include("config.jl")
include("../01_config/general_config.jl")

using PyCall
#@pyimport matplotlib as mpl
#mpl.use("Agg")

using PyPlot

using NCDatasets
using Formatting
using Statistics: mean, std

# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
img_dir  = normpath(joinpath( dirname(@__FILE__), "..", "..", "img"))
data_filename = joinpath(data_dir, "HMC_NCAR_5deg_init-omlmax_c4_s1000_w200_pdf_explore_lat035.jld")

using JLD
β = load(data_filename)["β"]

h = β[:, :, 1:12, :]

fig, axes = plt[:subplots](3, 4, sharex=true, sharey=true, figsize=(12, 6))
for t = 1:12
    ax = axes[t]
    d = h[:, :, t, :]
    ax[:hist](d[isfinite.(d)], bins=collect(range(0, stop=500, step=5)), density=1)
    ax[:text](0.9, 0.8, format("{:02d}", t), va="top", ha="right", transform=ax[:transAxes], size=15)
end

fig[:text](0.5,  0.05, "Mixed-Layer Depth [\$\\mathrm{m}\$]", va="top", ha="center", size=15)
fig[:text](0.05, 0.5,  "PDF", rotation=90, va="center", ha="right", size=15)
fig[:suptitle](format("File: {}\nLatitude: {:.2f}", basename(data_filename), lat[lat_i]))

axes[1][:set_xlim]([-20, 500])
axes[1][:set_xticks]([0,1,2,3,4,5] *100)

filename = joinpath(img_dir, format("{}-lat{:.2f}-all_lon.png", splitext(basename(data_filename))[1], lat[lat_i]))
#filename = joinpath(img_dir, format("{}-lat{:.2f}-lon{:.2f}.png", splitext(basename(data_filename))[1], lat[lat_i], lon[lon_i]))
fig[:savefig](filename)
plt[:show]()
