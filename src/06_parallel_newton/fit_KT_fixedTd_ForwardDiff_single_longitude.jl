include("./core/Newton_approach_fixedTd_ForwardDiff_core.jl")
include("./core/Param_Control.jl")
include("config.jl")

using .ParamControl
using .NewtonApproach
using .NewtonMethod
using Formatting
using Statistics: mean, std

if length(ARGS) != 1
    throw(ErrorException("The first argument should be the longitude index."))
end

lon_i = parse(Int, ARGS[1])


N = (years-2)*12
period = 12
beg_t  = period + 1

θd = 273.15 * ρ * c_p
init_h = zeros(dtype, 12) 
init_Q = zeros(dtype, 12)

β = zeros(dtype, length(lat), 24)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds")
B = zeros(dtype, length(time))

β .= missing_value

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

program_beg_time = Base.time()
for j = 1:length(lat)
    if isnan(θ[lon_i, j, 1])
        continue
    end

    @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", lon_i, lon[lon_i], j, lat[j])
    writelog(lon_i, "Doing lon[{:d}]: {:.2f}, lat[{:d}], {:.2f}\n", lon_i, lon[lon_i], j, lat[j])

    test_β  = zeros(dtype, 2*period)

    p_ctl = ParamControl.param_controller(param_path, fail_count_max)

    continue_flag = ParamControl.STATUS_CONTINUE

    while continue_flag == ParamControl.STATUS_CONTINUE
        
        try
            fit_β = NewtonApproach.fit(;
                N            = N,
                period       = period,
                beg_t        = beg_t,
                Δt           = Δt,
                init_h       = test_β[1:period],
                init_Q       = test_β[period+1:2*period],
                θd           = θd,
                θ            = θ[lon_i, j, :],
                S            = F[lon_i, j, :],
                B            = B,
                a            = p_ctl.test_param[1],
                max          = newton_fail_max,
                η            = newton_η,
                σ_ϵ          = 10.0,
                σ_Q          = 100.0,
                σ_h          = 1.0,
                verbose      = true
            )

            continue_flag = ParamControl.iterate_and_adjust!(p_ctl, true)
            
            test_β[:] = fit_β

        catch err
            if isa(err, NewtonMethod.NotConvergeException)
                println("Does not converge for a = ", p_ctl.test_param[1])
            else
                throw(err)
            end
            continue_flag = ParamControl.iterate_and_adjust!(p_ctl, false)
        end
                    
    end

    # Converge successfully
    if continue_flag == ParamControl.STATUS_ACCEPT
        println("Converge and reach acceptable a = ", p_ctl.test_param[end])
        β[j, 1:24] = test_β
    else
        println("Fail to converge to a = ", p_ctl.test_param[end])
        β[j, 1:24] .= NaN
    end
end

using JLD

println("Output filename: ", filename)
save(filename, Dict("β_bestfit" => β))

program_end_time = Base.time()

@printf("Total time used: %.2f min for %d points.\n",
    (program_end_time-program_beg_time)/ 60.0,
    length(lat),
)
