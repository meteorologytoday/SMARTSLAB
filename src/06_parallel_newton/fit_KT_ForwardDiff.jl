include("./core/Newton_approach_ForwardDiff_core.jl")
include("./core/Param_Control.jl")

include("../01_config/general_config.jl")
include("config_a_lon.jl")

using .NewtonApproach
using .NewtonMethod
using Formatting
import Statistics.mean

if length(ARGS) != 2
    throw(ErrorException("Length of ARGS must be 2. The first is the name of the run and the second is the longitude index."))
end

exp_name = ARGS[1]
lon_i = parse(Int, ARGS[2])

# construct data folder
main_dir = joinpath(data_path, splitext(basename(@__FILE__))[1], exp_name)
tmp_dir = joinpath(main_dir, "tmp", format("{:03d}", lon_i))
mkpath(main_dir)
mkpath(tmp_dir)

filename = format("{:03d}.jld", lon_i)
filename = joinpath(main_dir, filename)

println("This program is going to fit lon[", lon_i, "] = ", lon[lon_i])
if isfile(filename)
    println("File ", filename, " already exists. End program.")
    exit()
end

# Newton setting
newton_fail_max = 100
newton_η = 1e-2

fail_count_max = 5

test_scenarios = Dict(
    "init_zero" => collect(range(1.0, stop=0.0, length=10)),
)

N = (years-2)*12
period = 12
beg_t  = period + 1

init_θd = 273.15 * ρ * c_p
init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)

β = zeros(dtype, length(lon), length(lat), 25)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds")
B = zeros(dtype, length(time))

β .= missing_value

for scenario in keys(test_scenarios)
    println("# Scenario: ", scenario)
    for j = 1:length(lat), i = 1:length(lon)

        if isnan(θ[i, j, 1]) || i != 156 || j != 110
            continue
        end

        @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", i, lon[i], j, lat[j])

        test_Q  = copy(init_Q)
        test_h  = copy(init_h)
        test_θd = init_θd

        p_ctl = param_controller(test_scenarios[scenario], fail_count_max)

        continue_flag = true
        while continue_flag
            
            try
                fit_β = NewtonApproach.fit(;
                    N            = N,
                    period       = period,
                    beg_t        = beg_t,
                    Δt           = Δt,
                    init_h       = test_h,
                    init_Q       = test_Q,
                    init_θd      = test_θd,
                    θ            = θ[i, j, :],
                    S            = F[i, j, :],
                    B            = B,
                    a            = p_ctl.test_param[1],
                    max          = newton_fail_max,
                    η            = newton_η,
                    verbose      = true
                )

                # Converge successfully
                if p_ctl.test_param[1] <= 1.0
                    println("Converge and reach acceptable a = ", p_ctl.test_param[1])
                    β[i, j, 1:24] = fit_β
                    #println(β[i, j, 1:period])
                    #println(β[i, j, period+1:2period])
                    #println(β[i, j, end]/ρ/c_p)
                end

                continue_flag = iterate_and_adjust!(p_ctl, true)
                
                test_h[:] = fit_β[1:period]
                test_Q[:] = fit_β[period+1: 2*period]
                test_θd   = fit_β[end]

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


    println("β[156, 110, :] = ", β[156, 110,:])
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
