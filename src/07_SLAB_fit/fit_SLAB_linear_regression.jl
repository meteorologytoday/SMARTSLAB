model_name = "NCAR_5deg"

include("../lib/fit_cores/linear_regression/SLAB.jl")
include("../01_config/general_config.jl")

using Formatting

period = 12

β = zeros(dtype, lon_len, lat_len, 24)

θ = readModelVar("tos") * ρ * c_p
F = readModelVar("hfds")

for i = 1:lon_len, j = 1:lat_len

    if j == 1
        @printf("Doing lon[%d]: %.2f, lat[%d], %.2f\n", i, lon[i], j, lat[j])
    end
    
    if isnan(θ[i, j, 1])
        β[i, j, :] .= NaN
        continue
    end

    _β = SLAB!(
        period = period,
        F  = F[i, j, :],
        θ  = θ[i, j, :],
        Δt = Δt,
        reinterpolate = false
    )

    β[i,j,       1 :   period] = circshift(_β[       1:  period], 1)
    β[i,j,period+1 : 2*period] = circshift(_β[period+1:2*period], 1)
 
end

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
