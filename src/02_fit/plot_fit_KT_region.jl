using PyPlot
using JLD
import Statistics: mean, std

include("../01_config/general_config.jl")

filename = joinpath(data_path, "NCAR-fit_KT_region.jl.jld")
data = load(filename)

for region_name in keys(regions)

    # Reading data
    #=
    filename = joinpath(data_path, "NCAR-fit_KT_region.jl.nc")
    ds = Dataset(filename, "r")
    d["h"] = reshape(ds[format("{}_{}", region_name, "h")][:], 12, :)
    d["Q"] = reshape(ds[format("{}_{}", region_name, "Q")][:], 12, :)
    close(ds)
    =#

    d = data[region_name]

    filename = joinpath(data_path, "NCAR-fit_KT_region.jl.nc")

    omlmax = reshape(readModelRegionVar(region_name, "omlmax"), 12, :)

    omlmax_mean = mean(omlmax; dims=(2,))
    omlmax_std  = std(omlmax; dims=(2,))

    println(omlmax_std)

    t = collect(1:12)

    fig, ax = plt[:subplots](2, 1, figsize=(8,10), sharex=true)

    fig[:suptitle](format("Ocean: {} \n Best fit \$(\\alpha, \\beta) = ({:.2f}, {:.2f})\$", region_name, d["final_ab"]...))




    ax[2][:set_xlabel]("Month")
    ax[2][:set_xlim](0, 13)
    ax[2][:set_xticks](t)

    ax[1][:set_ylim](-10, 200)
    ax[2][:set_ylim](-3000, 3000)


    ax[1][:set_yticks](collect(0:20:200))
    ax[2][:set_yticks](collect(-3000:500:3000))

    ax[1][:grid](which="major", alpha=0.5)
    ax[2][:grid](which="major", alpha=0.5)
    
    ax[1][:invert_yaxis]()

    # Plot MLD
    ax[1][:set_ylabel]("MLD [m]")
    ax[1][:plot](t, d["h"], color="k", label="calculated")
    ax[1][:scatter](t, d["h"], color="k")

    ax[1][:plot](t, omlmax_mean, color="k", dashes=(5,2), label="model_output")
    ax[1][:errorbar](t, omlmax_mean, yerr=omlmax_std, color="k", fmt="o")
    ax[1][:legend]()

    # Plot Q flux
    ax[2][:set_ylabel]("Q [ \$ \\mathrm{W} \\, \\mathrm{m}^{-2} \$ ]")
    ax[2][:plot](t, d["Q"], color="k", label="calculated")

    ax[2][:legend]()

    imgname = joinpath(img_path, format("region_fit-{}-{}.png", model_name, region_name))
    @printf("Save image: %s  ...", imgname)
    fig[:savefig](imgname, dpi=100)
    println("done.")

    plt[:close](fig)

end

