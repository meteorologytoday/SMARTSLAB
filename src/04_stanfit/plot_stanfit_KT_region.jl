using Printf
@printf("Importing libraries...")
using PyPlot
using JLD
import Statistics: mean, std
println(" done.")
include("../01_config/general_config.jl")

for region_name in keys(regions)
    if region_name != "SPAC-1"
        continue
    end
    # Reading data

    filename = joinpath(data_path, format("NCAR-{}-stanfit_KT_region.jl.jld", region_name))
    data = load(filename)

    model_d = Dict()
    for varname in ["omlmax", "tos", "hfds"]
        tmp = reshape(readModelRegionVar(region_name, varname), 12, :)
        model_d[varname] = Dict(
            "mean" => mean(tmp; dims=(2,)),
            "std"  => std(tmp; dims=(2,))
        )
    end

    model_d["tos"]["mean"] .-= 273.15

    t = collect(1:12)

    fig, ax = plt[:subplots](3, 1, figsize=(8,12), sharex=true)

    fig[:suptitle](
        format("Ocean: {} \n \$T_d = {:.2f} \\pm {:.2f} \$",
            region_name,
            data["Td_mean"],
            data["Td_std"]
        )
    )


    ax1_twin = ax[1][:twinx]()

    ax[3][:set_xlabel]("Month")
    ax[3][:set_xlim](0, 13)
    ax[3][:set_xticks](t)

    ax[1][:set_ylim](-5, 35)
    ax1_twin[:set_ylim](-200, 200)
    ax[2][:set_ylim](-10, 200)
    ax[3][:set_ylim](-50, 50)


    ax[1][:set_yticks](collect(-5:5:35))
    ax1_twin[:set_yticks](collect(-200:50:200))
    ax[2][:set_yticks](collect(0:20:200))
    ax[3][:set_yticks](collect(-50:10:50))

    ax[1][:grid](which="major", alpha=0.5)
    ax[2][:grid](which="major", alpha=0.5)
    ax[3][:grid](which="major", alpha=0.5)
    
    ax[2][:invert_yaxis]()


    ax[1][:set_ylabel]("SST [deg C]")
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
    ax[2][:errorbar](t, data["h_mean"], yerr=data["h_std"], linestyle="-", color="r", fmt="s", label="HMC")
    ax[2][:errorbar](t, model_d["omlmax"]["mean"], yerr=model_d["omlmax"]["std"], color="k", dashes=(5,2), fmt="o", label=format("{}-omlmax", model_name))
    ax[2][:legend]()

    # Plot Q flux
    ax[3][:errorbar](t, data["Q_s_mean"], yerr=data["Q_s_std"], color="r", fmt="s", label="HMC")
    ax[3][:legend]()

    imgname = joinpath(img_path, format("stanfit_region_fit-{}-{}.png", model_name, region_name))
#    @printf("Save image: %s  ...", imgname)
#    fig[:savefig](imgname, dpi=100)
#    println("done.")

    plt[:show]()
    plt[:close](fig)

end

