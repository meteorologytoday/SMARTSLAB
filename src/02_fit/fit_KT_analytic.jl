include("../01_config/general_config.jl")
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
    "init_zero" => collect(range(1.0, stop=0.0, length=50))
)

bundle = NewtonApproach.Bundle(dtype; N=(years-2)*12, period=12, Δt=Δt)

init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)
β = zeros(dtype, length(lon), length(lat), 24)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds") * ρ * c_p
B = zeros(dtype, length(time))

β .= missing_value
for scenario in keys(test_scenarios)
    println("# Scenario: ", scenario)
    for j = 1:length(lat), i = 1:length(lon)

        if isnan(θ[i, j, 1]) || rand() < 0.9999
            continue
        end

        @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", i, lon[i], j, lat[j])

        init_Q .= 0.0
        init_h .= 0.0

        p_ctl = param_controller(test_scenarios[scenario], fail_count_max)

        continue_flag = true
        #println("θ", θ[i, j, :])
        #println("F", F[i, j, :])
        while continue_flag
            
            try
                h, Q = NewtonApproach.fit(;
                    bundle       = bundle,
                    init_h       = init_h,
                    init_Q       = init_Q,
                    θ            = θ[i, j, :],
                    S            = F[i, j, :],
                    B            = B,
                    θd           = 273.15,
                    a            = p_ctl.test_param[1],
                    max          = newton_fail_max,
                    η            = newton_η,
                    verbose      = false
                )

                println(p_ctl.test_param[1], h)

                # Converge successfully
                if p_ctl.test_param[1] < 1.5
                    println("Converge for a = ", p_ctl.test_param[1])
                    println(h)
                    β[i, j, 1:12] = h
                    β[i, j, 13:24] = Q
                    println(β[i, j, :])
                end

                continue_flag = iterate_and_adjust!(p_ctl, true)
                
                init_h[:] = h
                init_Q[:] = Q

            catch err
                if isa(err, NewtonMethod.NotConvergeException)
                    println("Does not converge for a = ", p_ctl.test_param[1])
                else
                    throw(err)
                end
                continue_flag = iterate_and_adjust!(p_ctl, false)
            end
            
        end
    end

    # output data
    filename = format("{}-{}-{}.nc", model_name, basename(@__FILE__), scenario)
    filename = joinpath(data_path, filename)

    ds = Dataset(filename, "c")

    println("Creating dimension...")
    defDim(ds, "time", 12)
    defDim(ds, "lat", length(lat))
    defDim(ds, "lon", length(lon))

    defVar(ds, "time", dtype, ("time",))[:] = collect(1:12)
    defVar(ds, "lat",  dtype, ("lat",))[:]  = lat
    defVar(ds, "lon",  dtype, ("lon",))[:]  = lon
     
    println("done")
    for o in (
        [
            "h", β[:, :, 1:12], Dict(
            "long_name"=>"Mixed-layer Depth",
            "units"=>"m",
            )
        ], [
            "Q", β[:, :,13:24], Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            )
        ],
    )
        varname, vardata, varatts = o
        println("Writing ", varname, " with size: ", size(vardata))

        ncvar = defVar(ds, varname, eltype(vardata), ("lon", "lat", "time"))
        ncvar.attrib["_FillValue"] = missing_value
        for key in keys(varatts)
            ncvar.attrib[key] = varatts[key]
        end

        ncvar[:] = vardata
        println("done.")
    end

    close(ds)
    println("Output file: ", filename)
end
