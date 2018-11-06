include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Formatting
import Statistics.mean

newton_fail_max = 10
newton_η = 1e-5

fail_count_max = 5

init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)
β = zeros(dtype, 24)

for scenario in keys(ab_pairs)
    println("# Scenario: ", scenario)

    data = Dict()
    for region_name in keys(regions)
        # Read model omlmax
        omlmax = reshape(readModelRegionVar(region_name, "omlmax"), 12, :)
        omlmax = mean(omlmax; dims=(2,))

        # Generate region data
        θ = readModelRegionVar(region_name, "tos") * (ρ * c_p)
        F = readModelRegionVar(region_name, "hfds")
        S = F 
        B = S * 0.0

        init_Q .=       0.0
        β      .= -999999.0
        
        if scenario == "init_omlmax"
            init_h[:] = omlmax
        else
            init_h .= 0
        end

        iter_obj = ab_pair_iterate_struct(ab_pairs[scenario], fail_count_max)

        continue_flag = true
        while continue_flag
            
            try
                h, Q = NewtonApproach.fit(
                    bundle       = bundle,
                    init_h       = init_h,
                    init_Q       = init_Q,
                    θ            = θ,
                    S            = S,
                    B            = B,
                    θd           = 273.15,
                    ab_pairs     = [iter_obj.ab_pair,],
                    max          = newton_fail_max,
                    η            = newton_η,
                    verbose      = false
                )

                # Converge successfully
                β[ 1:12] = h
                β[13:24] = Q
                println("Converge for pair: ", iter_obj.ab_pair)
                continue_flag = ab_pair_iterate!(iter_obj, true)
                
                init_h[:] = h
                init_Q[:] = Q

            catch err
                if isa(err, NewtonMethod.NotConvergeException)
                    #println("Does not converge.")
                else
                    println("Unknown Error")
                    println(err)
                    exit()
                end
                continue_flag = ab_pair_iterate!(iter_obj, false)
            end
            
        end
        data[region_name] = Dict(
            "h" => copy(β[ 1:12]),
            "Q" => copy(β[13:24]),
            "final_ab" => copy(iter_obj.ab_pair),
        )

    end

    using JLD

    filename = format("{}-{}-{}.jld", model_name, basename(@__FILE__), scenario)
    filename = joinpath(data_path, filename)
    println("Output filename: ", filename)

    rm(filename, force=true)
    save(filename, data)

end

