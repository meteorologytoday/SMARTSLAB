using PyPlot
using JLD
import Statistics: mean, std

include("../01_config/general_config.jl")

for region_name in keys(regions)

    # Reading data
    #=
    filename = joinpath(data_path, "NCAR-fit_KT_region.jl.nc")
    ds = Dataset(filename, "r")
    d["h"] = reshape(ds[format("{}_{}", region_name, "h")][:], 12, :)
    d["Q"] = reshape(ds[format("{}_{}", region_name, "Q")][:], 12, :)
    close(ds)
    =#
    data = Dict()
    for scenario in ["degenerate", "init_omlmax", "init_zero"]
        filename = joinpath(data_path, format("NCAR-fit_KT_region_analytic.jl-{}.jld", scenario))
        data[scenario] = load(filename)
    end

    model_d = Dict()
    for varname in ["omlmax", "tos", "hfds"]
        tmp = reshape(readModelRegionVar(region_name, varname), 12, :)
        model_d[varname] = Dict(
            "mean" => mean(tmp; dims=(2,)),
            "std"  => std(tmp; dims=(2,))
        )
    end

    model_d["tos"]["mean"] .-= 273.15

    init_omlmax_converge = data["init_omlmax"][region_name]["h"][1] != -999999.0

    t = collect(1:12)

    fig, ax = plt[:subplots](3, 1, figsize=(8,12), sharex=true)

    fig[:suptitle](
        format("Ocean: {} \n Best fit \$\\alpha = \$\ndegenerate: {:.2e}\ninit_zero: {:.2e}\ninit_omlmax: {}",
            region_name,
            data["degenerate"][region_name]["final_a"]...,
            data["init_zero"][region_name]["final_a"]...,
            (init_omlmax_converge) ? format("{:.2e}", data["init_omlmax"][region_name]["final_a"]...) : "not converge"
        )
    )


    ax1_twin = ax[1][:twinx]()

    ax[3][:set_xlabel]("Month")
    ax[3][:set_xlim](0, 13)
    ax[3][:set_xticks](t)

    ax[1][:set_ylim](-5, 35)
    ax1_twin[:set_ylim](-200, 200)
    ax[2][:set_ylim](-10, 200)
    ax[3][:set_ylim](-3000, 3000)


    ax[1][:set_yticks](collect(-5:5:35))
    ax1_twin[:set_yticks](collect(-200:50:200))
    ax[2][:set_yticks](collect(0:20:200))
    ax[3][:set_yticks](collect(-5000:1000:5000))

    ax[1][:grid](which="major", alpha=0.5)
    ax[2][:grid](which="major", alpha=0.5)
    ax[3][:grid](which="major", alpha=0.5)
    
    ax[2][:invert_yaxis]()


    ax[1][:set_ylabel]("SST [K]")
    ax1_twin[:set_ylabel]("Total Downward Heat Flux [\$ \\mathrm{W} \\, \\mathrm{m}^{-2} \$ ]")
    ax[2][:set_ylabel]("MLD [m]")
    ax[3][:set_ylabel]("Q [ \$ \\mathrm{W} \\, \\mathrm{m}^{-2} \$ ]")

    ms=5

    # Plot tos and hfds
    ln1 = ax[1][:errorbar](t, model_d["tos"]["mean"], yerr=model_d["tos"]["std"], color="k", linestyle="-", fmt="o", label="SST")
    ln2 = ax1_twin[:errorbar](t, model_d["hfds"]["mean"], yerr=model_d["hfds"]["std"], color="r", dashes=(5,2), fmt="s", label="F")
    
    ax[1][:legend]((ln1, ln2), ("SST", "F"))
     

    # Plot MLD
    #ax[2][:plot](t, degenerate["h"], color="b", marker="o", markersize=ms, dashes=(5,2,2,2), label="degenerate")
    ax[2][:plot](t, data["degenerate"][region_name]["h"], color="k", marker="o", markersize=ms, label="degenerate")
    ax[2][:plot](t, data["init_zero"][region_name]["h"], color="b", marker="s", markersize=ms, label="init_zero")
    ax[2][:plot](t, data["init_omlmax"][region_name]["h"], color="r", marker="d", markersize=ms, label=format("init_omlmax{}", (init_omlmax_converge) ? "" : "(not converge)"))


    ax[2][:errorbar](t, model_d["omlmax"]["mean"], yerr=model_d["omlmax"]["std"], color="k", dashes=(5,2), fmt="s", label=format("{}-omlmax", model_name))

    ax[2][:legend]()

    # Plot Q flux
    #ax[3][:plot](t, degenerate["Q"], color="b", marker="o", markersize=ms, dashes=(5,2,2,2), label="degenerate")
    ax[3][:plot](t, data["degenerate"][region_name]["Q"], color="k", marker="o", markersize=ms, label="degenerate")
    ax[3][:plot](t, data["init_zero"][region_name]["Q"], color="b", marker="s", markersize=ms, label="init_zero")
    ax[3][:plot](t, data["init_omlmax"][region_name]["Q"], color="r", marker="d", markersize=ms, label=format("init_omlmax{}", (init_omlmax_converge) ? "" : "(not converge)"))


    ax[3][:legend]()

    imgname = joinpath(img_path, format("region_fit-{}-{}_analytic.png", model_name, region_name))
    @printf("Save image: %s  ...", imgname)
    fig[:savefig](imgname, dpi=100)
    println("done.")

    plt[:close](fig)

end

