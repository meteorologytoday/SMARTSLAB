model_name = "NCAR_5deg"

include("../lib/fit_cores/linear_regression/SLAB.jl")
include("../lib/fit_cores/bayesian/Newton_SLAB.jl")
include("../01_config/general_config.jl")

using Printf
using Formatting
using .BayesianNewtonSLAB
using .NewtonMethod

period = 12

# Newton setting
newton_fail_max = 100
newton_η = 1e-2
fail_count_max = 5

# Posterior control
σ_ϵ          = 10.0
σ_Q          = 100.0
σ_h          = 1.0
h_rng        = [0, 5000.0]
verbose      = false

β = zeros(dtype, lon_len, lat_len, 24)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds")

time_beg = Base.time()
for i = 1:lon_len, j = 1:lat_len

    if j == 1
        @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", i, lon[i], j, lat[j])
    end
    
    if isnan(θ[i, j, 1])
        β[i, j, :] .= NaN
        continue
    end
    
    β[i, j, :] = SLAB!(
        period = period,
        F  = F[i, j, :],
        θ  = θ[i, j, :],
        Δt = Δt,
        reinterpolate = false
    )

    try 
        _β = BayesianNewtonSLAB.fit(;
                    N            = N,
                    period       = period,
                    beg_t        = beg_t,
                    Δt           = Δt,
                    init_h       = β[i, j,        1 :   period],
                    init_Q       = β[i, j, period+1 : 2*period],
                    θ            = θ[i, j, :],
                    F            = F[i, j, :],
                    max          = newton_fail_max,
                    η            = newton_η,
                    σ_ϵ          = σ_ϵ,
                    σ_Q          = σ_Q,
                    σ_h          = σ_h,
                    h_rng        = h_rng,
                    verbose      = verbose
        )

        # since now the last month (12.5) is actually 0.5, we
        # need to shift β

        β[i,j,       1 :   period] = circshift(_β[       1:  period], 1)
        β[i,j,period+1 : 2*period] = circshift(_β[period+1:2*period], 1)
    
    catch err
        if isa(err, NewtonMethod.NotConvergeException)
            println("Does not converge.")
        else
            throw(err)
        end
    end



end
time_end = Base.time()

print(format("The fitting uses {:.2f} min.\n", (time_end - time_beg) / 60.0))


nan2missing!(β)

output_vars = [
    [
        "h", β[:, :, 1:12], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Mixed-layer Depth",
        "units"=>"m",
        )
    ], [
        "Q", β[:, :,13:24], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Q-flux",
        "units"=>"W / m^2",
        )
    ],
]

filename = normpath(joinpath(
    data_path,
    format("{}.nc", splitext(basename(@__FILE__))[1])
))

ds = Dataset(filename,"c")
defDim(ds,"time", 12)
defDim(ds,"lat", length(lat))
defDim(ds,"lon", length(lon))

defVar(ds, "time", Float64, ("time",))[:] = collect(1:12) .- 0.5
defVar(ds, "lat", Float64, ("lat",))[:] = lat
defVar(ds, "lon", Float64, ("lon",))[:] = lon

for o in output_vars
    varname, vardata, vardims, varatts = o
    println("Writing ", varname, " with size: ", size(vardata), " ; dim: ", vardims)

    ncvar = defVar(ds, varname, eltype(vardata), vardims)
    ncvar.attrib["_FillValue"] = missing_value
    for key in keys(varatts)
        ncvar.attrib[key] = varatts[key]
    end

    ncvar[:] = vardata
    println("done.")
end

close(ds)
println("Output file: ", filename)
