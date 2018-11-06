include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Formatting
import Statistics.mean

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan, Mamba
@printf("done\n")

model_script = read(joinpath(dirname(@__FILE__), "KT.stan"), String)

@printf("Now we are going to build stan model...\n")
stanmodel = Stanmodel(name="KT", nchains=4, model=model_script)
display(stanmodel)



data = Dict()
for region_name in keys(regions)
    # Read model omlmax
    omlmax = reshape(readModelRegionVar(region_name, "omlmax"), 12, :)
    omlmax = mean(omlmax; dims=(2,))

    # Generate region data
    θ = readModelRegionVar(region_name, "tos") * (ρ * c_p)
    F = readModelRegionVar(region_name, "hfds")
    N = length(θ) 


    data = Dict(
       "N"       => N, 
       "period"  => 12, 
       "dt"      => Δt, 
       "theta"   => θ, 
       "F"       => F,
    )
    
    rc, sim1 = stan(stanmodel, [data], "tmp", CmdStanDir=ENV["CMDSTAN_HOME"])

    if rc == 0
        println("Subset Sampler Output")
        sim = sim1[1:1000, :, :]
        describe(sim)
    else
        println("There are errors!!")
    end
    exit()
end
#=
using JLD

filename = format("{}-{}-{}.jld", model_name, basename(@__FILE__), scenario)
filename = joinpath(data_path, filename)
println("Output filename: ", filename)

rm(filename, force=true)
save(filename, data)
=#

