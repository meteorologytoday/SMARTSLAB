include("config.jl")

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
data_filename = joinpath(data_dir, "HMC_NCAR_5deg_init-omlmax_c4_s100_w20_pdf_explore_lat035.jld")

using JLD
β = load(data_filename)["β"]

h = β[:, :, 1:12, :]

fig, axes = plt[:subplots](3, 4, figsize=(12, 6))
for t = 1:12
    ax = axes[t]
    d = h[:, :, t, :]
    ax[:hist](d[isfinite.(d)], bins=collect(range(0, stop=300, step=10)))
    ax[:set_title](format("{:02d}", t))
end

fig[:suptitle](format("File: {}\nLatitude: {:.2f}", basename(data_filename), lat[lat_i]))
plt[:show]()
