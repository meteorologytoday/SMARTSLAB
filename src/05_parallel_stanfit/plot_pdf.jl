using PyCall
#@pyimport matplotlib as mpl
#mpl.use("Agg")

using PyPlot

@pyimport mpl_toolkits.basemap as basemap
@pyimport mpl_toolkits.basemap.cm as cm

using NCDatasets
using Formatting
using Statistics: mean, std

include("config_plot.jl")


function nanmean(x; dims)
    n = isnan.(x)

    y = copy(x)
    y[n] .= 0.0

    return sum(y; dims=dims) ./ sum(1 .- n; dims=dims)
end

lat_i = 36

h_mean  = nanmean(data_h; dims=(1,))[1,:, :]


lab_lats = collect(-90:30: 90)

# setting maps
filename = normpath(joinpath(
    img_dir,
    format("{}-pdf.png", splitext(basename(nc_filename))[1])
))

fig, axes = plt[:subplots](3, 4, figsize=(12, 6))
println(size(lat))
println(lat)
for t = 1:12
    ax = axes[t]
    d = data_h[:, lat_i, t]
    ax[:hist](d[isfinite.(d)], bins=collect(range(0, stop=150, step=10)))
    ax[:set_title](format("{:02d}", t))
end

fig[:suptitle](format("File: {}\nLatitude: {:.2f}", basename(nc_filename), lat[lat_i]))
plt[:show]()
