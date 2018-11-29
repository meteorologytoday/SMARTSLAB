using PyCall
@pyimport matplotlib as mpl
mpl.use("Agg")

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

#=
function nanmean(x; dims)
    xdev_sq = x -
    return sum(x; dims=dims) ./ sum(m; dims=dims)
end

=#

#complex_model_filename = joinpath(data_dir, "NCAR_CESM1-WACCM", "SMART_5deg_omlmax_Omon_CESM1-WACCM_piControl_r1i1p1_009601-029512.nc")

lon_rng = range(1, step=1, stop=length(lon)) |> collect

h_rng = [0, 200]
Q_rng = [-100, 100]
Td_rng = [-15, 25]

lat_rng = [-90, 90]

h_mean  = nanmean(data_h; dims=(1,))[1,:, :]
Q_mean  = nanmean(data_Q; dims=(1,))[1,:, :]
#Td_mean = nanmean(data_Td; dims=(1,))[1,:]

lab_lats = collect(-90:30: 90)

#=
### plot Td ###
println("Plotting Td")
# setting maps
filename = normpath(joinpath(
    img_dir,
    format("{}-slice-Td.png", splitext(basename(nc_filename))[1])
))

fig, ax = plt[:subplots](1, 1, figsize=(12, 6))

ax[:set_xticks](lab_lats)
ax[:set_xlabel]("Latitude [\$\\mathrm{deg}\$]")
ax[:set_ylabel]("Deep Ocean Temperature \$T_\\mathrm{d}\$ [\$ {}^\\circ\\mathrm{C}\$]")
ax[:set_xlim](lat_rng)
ax[:set_ylim](Td_rng)

for j = 1:length(lon_rng)
    ax[:plot](lat, data_Td[lon_rng[j], :])
end
ax[:plot](lat, Td_mean, color="k", linestyle="--", linewidth=2, label="mean")

ax[:legend]()
ax[:grid]()

fig[:suptitle](format("File: {}", basename(nc_filename)))


#plt[:show]()
fig[:savefig](filename, dpi=100)
println("Output file:", filename)
plt[:close](fig)
=#

### plot h and Q ###
for i = 1:size(data_h)[3]
    println("Plotting month ", i)
    # setting maps
    filename = normpath(joinpath(
        img_dir,
        format("{}-slice-{:02d}.png", splitext(basename(nc_filename))[1], i)
    ))

    fig, axes = plt[:subplots](2, 1, figsize=(10, 6), sharex=true)

    axes[2][:set_xticks](lab_lats)
    axes[2][:set_xlabel]("Latitude [\$\\mathrm{deg}\$]")
    
    axes[1][:set_ylabel]("MLD [\$\\mathrm{m}\$]")
    axes[2][:set_ylabel]("Q flux [\$\\mathrm{W} \\, \\mathrm{m}^{-2} \$]")

   
    axes[1][:set_xlim](lat_rng)
    axes[1][:set_ylim](h_rng)
    axes[2][:set_ylim](Q_rng)

    axes[1][:invert_yaxis]()

    println(size(data_h))
    for j = 1:length(lon_rng)
        axes[1][:plot](lat, data_h[lon_rng[j], :, i])
        axes[2][:plot](lat, data_Q[lon_rng[j], :, i])
    end
 
    axes[1][:plot](lat, h_mean[:, i], color="k", linestyle="--", linewidth=2, label="mean")
    axes[2][:plot](lat, Q_mean[:, i], color="k", linestyle="--", linewidth=2, label="mean")

    axes[1][:legend]()
    axes[2][:legend]()

    axes[1][:grid]()
    axes[2][:grid]()
    
    fig[:suptitle](format("File: {} \n Month: {:02d}", basename(nc_filename), i))


    #plt[:show]()
    fig[:savefig](filename, dpi=100)
    println("Output file:", filename)
    plt[:close](fig)
end
