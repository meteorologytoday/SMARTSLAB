include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Printf
using Formatting
import Statistics.mean

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan, Mamba
@printf("done\n")

model_script = read(joinpath(dirname(@__FILE__), "KT.stan"), String)

@printf("Now we are going to build stan model...\n")
nchains     = 10
num_samples = 1000
stanmodel = Stanmodel(name="KT", nchains=nchains, num_samples=num_samples, model=model_script)
display(stanmodel)


h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q_s.{}", i) for i = 1:12 ]

data = Dict()

beg_time = Base.time()
println("Start counting time...")
for region_name in keys(regions)

    println("### Doing region: ", region_name)

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

    println("Doing STAN...")
    rc, sim1 = stan(stanmodel, [data], "tmp", CmdStanDir=ENV["CMDSTAN_HOME"])
    if rc != 0
        println("There are errors!!")
        continue
    end

    println("Extracting result...")
    h_mean = zeros(12)
    h_std  = zeros(12)
    Q_mean = zeros(12)
    Q_std  = zeros(12)


    data_h = sim1[:, h_key, :].value
    data_Q = sim1[:, Q_key, :].value

    for i = 1:12
        h_mean[i] = mean(data_h[:, i, :])
        h_std[i]  = std(data_h[:, i, :])

        Q_mean[i] = mean(data_Q[:, i, :])
        Q_std[i]  = std(data_Q[:, i, :])
    end


    using JLD
    filename = format("{}-{}-{}.jld", model_name, region_name, basename(@__FILE__))
    filename = joinpath(data_path, filename)
    println("Output filename: ", filename)
    rm(filename, force=true)
    save(filename, Dict(
        "h_mean" => h_mean,
        "h_std"  => h_std,
        "Q_s_mean" => Q_mean,
        "Q_s_std"  => Q_std,
    ))

    break
end

end_time = Base.time()

@printf("Total time used: %.2f min for %d regions, with nchains = %d, num_samples = %d\n",
    (end_time-beg_time)/ 60.0,
    length(keys(regions)),
    nchains,
    num_samples
)
