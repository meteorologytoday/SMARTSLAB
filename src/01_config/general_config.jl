using NCDatasets
using Printf
using Formatting
@printf("# Running %s\n", basename(@__FILE__))

include("models.jl")
include("constants.jl")
include("regions.jl")

dtype = Float64
data_path = joinpath(dirname(@__FILE__), "..", "..", "data")
img_path = joinpath(dirname(@__FILE__), "..", "..", "img")

fn = Dict()
fn_region = Dict()

try
    model_name
    models[model_name]
catch e

    if isa(e, UndefVarError)
        global model_name = "NCAR"
        println("model_name not set. use the deafult: ", model_name)
    end
end



model_fname_pattern = models[model_name]
println("# Files on the list:")
varnames = ["omlmax", "tos", "hfds"]
for varname in varnames
    fn[varname] = joinpath(data_path, (@eval @sprintf($model_fname_pattern, $varname)))
    println("   ", fn[varname], " -- ", (isfile(fn[varname]) ? "EXISTS" : "MISSING"))
end
for region_name in keys(regions)
    fn_region[region_name] = joinpath(data_path, format("{}-region-{}.nc", model_name, region_name))
    println("   ", fn_region[region_name], " -- ", (isfile(fn_region[region_name]) ? "EXISTS" : "MISSING"))
end


println()

function missing2nan!(x)
    x[x .== missing_value] .= NaN
end

function nan2missing!(x)
    x[isnan.(x)] .= missing_value
end


function readModelVar(varname, range_tuple=())
    local var
    @printf("# Reading var [%s]...", varname)
    ds = Dataset(fn[varname],"r")
    v = ds[varname]

    if range_tuple == ()
        range_tuple = repeat([Colon(),], outer=(length(size(v)),))
    end

    v = v[range_tuple...]
    v = nomissing(v, NaN)
    var = convert(Array{dtype}, v)
    close(ds)

    println("done.")

    return var
end

function readModelRegionVar(region_name, varname, range_tuple=())
    local var
    @printf("# Reading [%s] [ %s]...", region_name, varname)

    ds = Dataset(fn_region[region_name], "r")
    v = ds[varname]

    if range_tuple == ()
        range_tuple = repeat([Colon(),], outer=(length(size(v)),))
    end

    v = v[range_tuple...]
    v = nomissing(v, NaN)
    var = convert(Array{dtype}, v)
    close(ds)

    println("done. Size: ", size(v))

    return var
end


function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%.2f ", A[i,j]) 

        end
        @printf("\n")
    end
end

ds = Dataset(fn["tos"],"r")
missing_value = convert(dtype, ds["tos"].attrib["_FillValue"])
lon  = ds["lon"][:]
lat  = ds["lat"][:]
time = ds["time"][:]

close(ds)



mon_secs = 365.0 / 12.0 * 86400.0
Δt  = 1.0 * mon_secs

lon_len = length(lon)
lat_len = length(lat)
time_len = length(time)


if time_len % 12 != 0
    error("There are $months months, not a multiple of 12")
end

years = Int(time_len / 12)
N     = (years-2) * 12 # Discard the first and final year
beg_t = 13


println("# Data info")
@printf("  Model  : %s\n", model_name)
@printf("  Length : %02d years.\n", years)
println("  Data type: ", dtype)
println("  Resolution: lon=", length(lon), "; lat=", length(lat), "; time=", length(time))
println()

println("# Constants")
@printf("ρ   : %.2f kg / m^3\n", ρ)
@printf("c_p : %.2f J / kg / K\n", c_p)
@printf("Δt  : %.2f s\n", Δt)
println()
