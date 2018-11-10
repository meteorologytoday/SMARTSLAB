program_beg_time = Base.time()

include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Printf
using Formatting
import Statistics.mean

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan
@printf("done\n")

model_script = read(joinpath(dirname(@__FILE__), "KT.stan"), String)

@printf("Now we are going to build stan model...\n")
nchains     = 1
num_samples = 1000
stanmodel = Stanmodel(
    Stan.Optimize();
    name="KT",
    nchains=nchains,
    num_samples=num_samples,
    model=model_script
)
display(stanmodel)


h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q_s.{}", i) for i = 1:12 ]

data = Dict()
time_stat = Dict()

for region_name in keys(regions)

    if region_name != "SPAC-1"
        continue
    end
    beg_time = Base.time()
    println("### Doing region: ", region_name)

    # Read model omlmax
    tmp = reshape(readModelRegionVar(region_name, "omlmax"), 12, :)
    omlmax = mean(tmp; dims=(2,)),

    # Generate region data
    θ = readModelRegionVar(region_name, "tos") * (ρ * c_p)
    F = readModelRegionVar(region_name, "hfds")
    N = length(θ) 

    data = Dict(
        "N"           => N, 
        "period"      => 12, 
        "dt"          => Δt, 
        "theta"       => θ, 
        "F"           => F,
        "epsilon_std" => 10.0,
        "Q_std"       => 100.0,
    )

    println("Doing STAN...")
    rc, sim1 = stan(
        stanmodel,
        [data],
        "tmp";
        CmdStanDir = ENV["CMDSTAN_HOME"],
        init = Dict(
            "h" => omlmax
        )
    )

    if rc != 0
        println("There are errors!!")
        continue
    end



    println("Extracting result...")

    v = sim1[1]["optimize"]

    h_est = zeros(12)
    Q_est = zeros(12)

    Td_est = v["theta_d"][1] / ρ / c_p

    for i = 1:12
        h_est[i] = v[h_key[i]][1] 
        Q_est[i] = v[Q_key[i]][1] 
    end

    time_stat[region_name] = Base.time() - beg_time
    println(format("Region {} took {:.2f} min to do STAN fit.", region_name, time_stat[region_name] / 60.0 ))

    using JLD
    filename = format("{}-{}-{}.jld", model_name, region_name, basename(@__FILE__))
    filename = joinpath(data_path, filename)
    println("Output filename: ", filename)
    rm(filename, force=true)
    save(filename, Dict(
        "h_est"   => h_est,
        "Q_s_est" => Q_est,
        "Td_est"  => Td_est,
    ))

end

program_end_time = Base.time()

@printf("Total time used: %.2f min for %d regions, with nchains = %d, num_samples = %d\n",
    (program_end_time-program_beg_time)/ 60.0,
    length(keys(regions)),
    nchains,
    num_samples
)
