using JLD
using Formatting
using NCDatasets

include("config.jl")

main_dir = joinpath(data_path, "fit_KT_fixedTd_ForwardDiff_single_longitude", exp_name)

β_bestfit = zeros(dtype, length(lon), length(lat), 24)

β_bestfit .= NaN

for i = 1:length(lon)
    filename = normpath(joinpath(main_dir, format("{:03d}.jld", i)))
    if isfile(filename)
        println("Doing file: ", filename)
        d = JLD.load(filename)
        β_bestfit[i, :, :] = d["β_bestfit"][:, :]
    else
        println("File: ", filename, " does not exist. Skip this one.")
    end
end

nan2missing!(β_bestfit)

output_vars = [
    [
        "h_bestfit", β_bestfit[:, :, 1:12], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Mixed-layer Depth",
        "units"=>"m",
        )
    ], [
        "Q_bestfit", β_bestfit[:, :,13:24], ("lon", "lat", "time"), Dict(
        "long_name"=>"Mean of Q-flux",
        "units"=>"W / m^2",
        )
    ],
]

if size(β_bestfit)[3] == 25
    push!(output_vars, [
        "Td_mean", β_bestfit[:, :, 25], ("lon", "lat"), Dict(
        "long_name"=>"Mean of Deep Ocean Temperature",
        "units"=>"K",
        )
    ])
end
 




filename = joinpath(data_path, format("{}.nc", exp_name))
ds = Dataset(filename,"c")
defDim(ds,"time", 12)
defDim(ds,"lat", length(lat))
defDim(ds,"lon", length(lon))

defVar(ds, "time", Float64, ("time",))[:] = collect(1:12)
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

