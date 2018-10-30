include("../01_config/general_config.jl")
include("../01_config/regions.jl")
include("./core/Newton_approach_core.jl")

using .NewtonApproach
using .NewtonMethod
using Formatting
#using NCDatasets

#filename = format("{}-{}.nc", model_name, basename(@__FILE__))
filename = format("{}-{}.jld", model_name, basename(@__FILE__))
filename = joinpath(data_path, filename)
println("Output filename: ", filename)

data = Dict()

newton_fail_max = 10
newton_η = 1e-5

fail_count_max = 5

bundle = NewtonApproach.Bundle(dtype; N=(years-2)*12, period=12, Δt=Δt)

ab_pairs = [
    1e-7 0.00;
    1e-7 0.50;
    1e-7 1.00;
     0.0 1.00;
]



mutable struct ab_pair_iterate_struct{T <: AbstractFloat}
    beg_i           :: Integer
    end_i           :: Integer
    ab_pairs        :: Array{T, 2}
    ab_pair         :: Array{T, 1}
    ab_pair_tmp     :: Array{T, 1}
    fail_count      :: Integer
    fail_count_max  :: Integer

    function ab_pair_iterate_struct(
        ab_pairs::Array{T, 2},
        fail_count_max::Integer
    ) where T <: AbstractFloat

        return new{T}(
            1,
            1,
            ab_pairs,
            copy(ab_pairs[1,:]),
            copy(ab_pairs[1,:]),
            0,
            fail_count_max
        )
    end
end

function ab_pair_iterate!(
    o::ab_pair_iterate_struct,
    converge_flag::Bool
)

    if converge_flag == true
        o.fail_count = 0
        if o.ab_pair == o.ab_pairs[o.end_i, :]
            #println("Touch endpoint")
            # If now touches the end point
            
            if o.end_i == size(o.ab_pairs)[1]
                # If it is entirely complete
                return false
            else
                # If not then move on to next interval
                o.beg_i, o.end_i = o.end_i, o.end_i + 1
                o.ab_pair[:] = o.ab_pairs[o.end_i, :]
                return true
            end
        else
            #println("Does not touch endpoint")
            #println(o.ab_pair)
            # If not touching end point
            o.ab_pair_tmp[:] = o.ab_pair
            o.ab_pair[:] = o.ab_pairs[o.end_i, :]

            return true
        end
    else
        # If not converge
        o.fail_count += 1
        if o.fail_count >= o.fail_count_max
            return false
        else
            o.ab_pair[:] = (o.ab_pair_tmp[:] + o.ab_pair) / 2.0
            #println("Fail count: ", fail_count, "; ab_pair change to :", ab_pair)
            return true
        end
    end
end


init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)
β = zeros(dtype, 24)

for region_name in keys(regions)
    global fail_count 
    # Generate region data
    θ = readModelRegionVar(region_name, "tos") * (ρ * c_p)
    F = readModelRegionVar(region_name, "hfds")
    S = F 
    B = S * 0.0

    init_h .=       0.0 
    init_Q .=       0.0
    β      .= -999999.0
    
    iter_obj = ab_pair_iterate_struct(ab_pairs, fail_count_max)

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

    println(β[1:12])

    data[region_name] = Dict(
        "h" => copy(β[ 1:12]),
        "Q" => copy(β[13:24]),
        "final_ab" => copy(iter_obj.ab_pair),
    )

end


rm(filename, force=true)

using JLD
save(filename, data)

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
