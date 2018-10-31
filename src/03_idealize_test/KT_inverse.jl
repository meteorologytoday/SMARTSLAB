include("KT_simulate.jl")

# Fitting 

years = Int(length(t_year) / 12.0)

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
        1e-7 1.00;
         0.0 1.00;
    ]
)


init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)
β = zeros(dtype, 24)

data = Dict()
for scenario in keys(ab_pairs)
    println("# Scenario: ", scenario)

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
                θ            = θs,
                S            = S,
                B            = B,
                θd           = Td,
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

    data[scenario] = Dict(
        "h" => copy(β[ 1:12]),
        "Q" => copy(β[13:24]),
    )
end

