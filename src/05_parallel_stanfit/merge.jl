using JLD
using Formatting
using NCDatasets

include("config.jl")
include("../01_config/general_config.jl")

main_dir = joinpath(data_path, "stanfit_KT_single_longitude", exp_name)

β_mean = zeros(dtype, length(lon), length(lat), 25)
β_std  = copy(β_mean)

β_mean .= NaN
β_std  .= NaN

for i = 1:length(lon)
    filename = normpath(joinpath(main_dir, format("{:03d}.jld", i)))
    if isfile(filename)
        println("Doing file: ", filename)
        d = JLD.load(filename)
        β_mean[i, :, :] = d["β_mean"][:, :]
        β_std[i, :, :]  = d["β_std"][:, :]
    else
        println("File: ", filename, " does not exist. Skip this one.")
    end
end

nan2missing!(β_mean)
nan2missing!(β_std)

filename = joinpath(data_path, format("{}.nc", exp_name))
ds = Dataset(filename,"c")
defDim(ds,"time", 12)
defDim(ds,"lat", length(lat))
defDim(ds,"lon", length(lon))

defVar(ds, "time", Float64, ("time",))[:] = collect(1:12)
defVar(ds, "lat", Float64, ("lat",))[:] = lat
defVar(ds, "lon", Float64, ("lon",))[:] = lon

for o in (
    [
        "h_mean", β_mean[:, :, 1:12], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Mixed-layer Depth",
        "units"=>"m",
        )
    ], [
        "Q_mean", β_mean[:, :,13:24], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Q-flux",
        "units"=>"W / m^2",
        )
    ], [
        "Td_mean", β_mean[:, :, 25], ("lon", "lat"), Dict(
        "long_name"=>"Mean of Deep Ocean Temperature",
        "units"=>"K",
        )
    ], [
        "h_std", β_std[:, :, 1:12], ("lon", "lat", "time"), Dict(
        "long_name"=>"Standard Deviation of Mixed-layer Depth",
        "units"=>"m",
        )
    ], [
        "Q_std", β_std[:, :,13:24], ("lon", "lat", "time"), Dict(
        "long_name"=>"Standard Deviation of Q-flux",
        "units"=>"W / m^2",
        )
    ], [
        "Td_std", β_std[:, :, 25], ("lon", "lat"), Dict(
        "long_name"=>"Standard Deviation of Deep Ocean Temperature",
        "units"=>"K",
        )
    ]
)
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

