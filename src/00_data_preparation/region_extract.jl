include("../01_config/general_config.jl")
include("../01_config/regions.jl")

using Formatting

function nanmean(x)
    mask = isfinite.(x)
    return sum(x[mask]) / sum(mask)
end


data = Dict()

# load data
for region_name in keys(regions)
    println("Region:", region_name)
    data[region_name] = Dict()
    mask = region_mask(lon, lat, region_name)
    local time_len = 0
    for varname in varnames
        println("Processing var: ", varname)
        x = readModelVar(varname)
        time_len = size(x)[3]
        new_x = zeros(eltype(x), time_len)
        
        for i in 1:length(new_x)
            new_x[i] = nanmean(x[:, :, i][mask])
        end

        data[region_name][varname] = new_x
    end


    filename = format("{}-region-{}.nc", model_name, region_name)
    filename = joinpath(data_path, filename)
    rm(filename, force=true)

    println("Output file: ", filename)
    
    ds = Dataset(filename,"c")
    defDim(ds,"time", time_len)

    time = defVar(ds, "time", Float64, ("time",))
    time[:] = convert(Array{Float64, 1}, collect(1:time_len))

    for varname in keys(data[region_name])

        var = data[region_name][varname]
        println("Writing ", varname, " with size: ", size(var))

        v = defVar(ds, varname, eltype(var), ("time",))
        v.attrib["_FillValue"] = missing_value
        v[:] = var
        println("Missing_value:", missing_value)


    end
    close(ds)

end





