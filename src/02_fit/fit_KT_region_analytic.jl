include("../01_config/general_config.jl")
include("../01_config/regions.jl")
include("./core/Newton_approach_analytic_core.jl")
include("./core/Param_Control.jl")

using .NewtonApproach
using .NewtonMethod
using Formatting
import Statistics.mean
#using NCDatasets

#filename = format("{}-{}.nc", model_name, basename(@__FILE__))

newton_fail_max = 20
newton_η = 1e-5

fail_count_max = 5

test_scenarios = Dict(
    "degenerate" => [
         1.00;
    ],

    "init_zero" => [
        1.00;
        0.75;
        0.50;
        0.25;
        0.00
    ],
    "init_omlmax" => [
        1e-5;
        1e-7;
         0.0;
    ]
)

bundle = NewtonApproach.Bundle(dtype; N=(years-2)*12, period=12, Δt=Δt)

init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)
β = zeros(dtype, 24)

for scenario in keys(test_scenarios)
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
        @printf("Initializing parameter controller... ")
        p_ctl = param_controller(test_scenarios[scenario], fail_count_max)
        @printf("done\n")

        continue_flag = true
        while continue_flag
            
            try
                #println(typeof(p_ctl.test_param))
                #println(p_ctl.test_param)
                h, Q = NewtonApproach.fit(;
                    bundle       = bundle,
                    init_h       = init_h,
                    init_Q       = init_Q,
                    θ            = θ,
                    S            = S,
                    B            = B,
                    θd           = 273.15,
                    a            = p_ctl.test_param[1],
                    max          = newton_fail_max,
                    η            = newton_η,
                    verbose      = false
                )

                # Converge successfully
                β[ 1:12] = h
                β[13:24] = Q
                println("Converge for a = ", p_ctl.test_param)
                continue_flag = iterate_and_adjust!(p_ctl, true)
                
                init_h[:] = h
                init_Q[:] = Q

            catch err
                if isa(err, NewtonMethod.NotConvergeException)
                    println("Does not converge.")
                else
                    throw(err)
                end
                continue_flag = iterate_and_adjust!(p_ctl, false)
            end
            
        end
        data[region_name] = Dict(
            "h" => copy(β[ 1:12]),
            "Q" => copy(β[13:24]),
            "final_a" => copy(p_ctl.test_param),
        )

    end

    using JLD

    filename = format("{}-{}-{}.jld", model_name, basename(@__FILE__), scenario)
    filename = joinpath(data_path, filename)
    println("Output filename: ", filename)

    rm(filename, force=true)
    save(filename, data)

end
