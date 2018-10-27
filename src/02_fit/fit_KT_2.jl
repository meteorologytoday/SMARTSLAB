include("../01_config/general_config.jl")
include("config.jl")

include("./core/Newton_approach_core.jl")

using .NewtonApproach
using .NewtonMethod
using Formatting

ab_pairs = [
    [1e-5, 0.00],
    [1e-5, 0.25],
    [1e-5, 0.30],
    [1e-5, 0.35],
    [1e-5, 0.40],
    [1e-5, 0.45],
    [1e-5, 0.50],
    [1e-5, 0.75],
    [1e-5, 1.00],
    [.8e-5, 1.00],
    [.6e-5, 1.00],
    [.4e-5, 1.00],
    [.2e-5, 1.00],
    [1e-6, 1.00],
    [.8e-6, 1.00],
    [.6e-6, 1.00],
    [.5e-6, 1.00],
    [.2e-6, 1.00],
]

β = zeros(dtype, lon_len, lat_len, 24)

S = F 
B = S * 0.0


bundle = NewtonApproach.Bundle(dtype; N=(years-2)*12, period=12, Δt=Δt)

selected_index = (125, 137)
for i = 1:lon_len, j = 1:lat_len
    global bundle

    #=
    if (i, j) != selected_index
        β[i, j, :] .= NaN
        #continue
    end
    =#
    
    if lat[j] > 40.0 || lat[j] < -40.0
        β[i, j, :] .= NaN
        continue
    end

    if isnan(θ[i, j, 1])
        β[i, j, :] .= NaN
        continue
    end

    @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", i, lon[i], j, lat[j])
    init_h = zeros(dtype, 12) 
    init_Q = zeros(dtype, 12) 
    try
        h, Q = NewtonApproach.fit(
            bundle       = bundle,
            init_h       = init_h,
            init_Q       = init_Q,
            θ            = θ[i, j, :],
            S            = S[i, j, :],
            B            = B[i, j, :],
            θd           = 273.15,
            ab_pairs     = ab_pairs,
            max          = 10,
            η            = 1e-5
        )

        β[i, j,  1:12] = h
        β[i, j, 13:24] = Q
        println("Converge.")
    catch err
        β[i, j, :] .= unconverge_value
        if isa(err, NewtonMethod.NotConvergeException)
            println("Does not converge.")
        else
            println("Unknown Error")
            println(err)
        end
    end

end

println("# h")
println(β[selected_index...,  1:12])

println("# Q")
println(β[selected_index..., 12:24])

nan2missing!(β)

time = collect(Float64, 1:12)
filename = format("{}-{}.nc", model_name, basename(@__FILE__))
filename = joinpath(data_path, filename)

rm(filename, force=true)

for obj in [
    [
        β[:, :,  1:12], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "_FillValue" => missing_value,
            "missing_value" => missing_value,
            "unconverge_value" => unconverge_value

        )
    ], [
        β[:, :, 13:24], "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "_FillValue" => missing_value,
            "missing_value" => missing_value,
            "unconverge_value" => unconverge_value
        )
    ]
]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    println("Writing ", varname, " with size: ", size(var))    

    nccreate(
        filename,
        varname,
        "lon", lon,
        "lat", lat,
        "time", 12,
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)

@printf("Output file: %s\n", filename)
