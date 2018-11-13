include("../01_config/general_config.jl")
include("../01_config/regions.jl")
include("./core/Newton_approach_core.jl")
include("./core/ab_training_toolkit.jl")

using .NewtonApproach
using .NewtonMethod
using Formatting
import Statistics.mean
#using NCDatasets

#filename = format("{}-{}.nc", model_name, basename(@__FILE__))

newton_fail_max = 10
newton_η = 1e-5

fail_count_max = 5

bundle = NewtonApproach.Bundle(dtype; N=(years-2)*12, period=12, Δt=Δt)


ab_pairs = Dict(
    "degenerate" => [
         1.0 0.00;
         1.0 0.00;
    ],

    "init_zero" => [
        1e-7 0.00;
        1e-7 0.50;
        1e-7 1.00;
        0.0 1.00;
    ],

    "init_omlmax" => [
        1e-7 0.50;
        1e-7 1.00;
         0.0 1.00;
    ]
)





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



#=
ds = Dataset(filename, "c")
defDim(ds, "time", 12)
defVar(ds, "time", Float64, ("time",))[:] = convert(Array{Float64, 1}, collect(1:12))

for region_name in keys(regions)
    for varname in ["h", "Q"]
        d = data[region_name][varname]
        new_varname = format("{}_{}", region_name, varname)
        var = defVar(ds, new_varname, eltype(d), ("time",) ; fillvalue=missing_value)
        var[:] = d
    end
end
close(ds)
=#
