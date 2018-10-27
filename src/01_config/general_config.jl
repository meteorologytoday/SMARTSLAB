using NetCDF
using Printf

@printf("# Running %s\n", basename(@__FILE__))

include("models.jl")
include("constants.jl")
include("regions.jl")

dtype = Float64
data_path = joinpath(dirname(@__FILE__), "..", "..", "data")
img_path = joinpath(dirname(@__FILE__), "..", "..", "img")

fn = Dict()

model_name = "NCAR"
model_fname_pattern = models[model_name]
println("# Files on the list:")
for varname in ["omlmax", "S", "B", "tos", "hfds"]
    fn[varname] = joinpath(data_path, (@eval @sprintf($model_fname_pattern, $varname)))
    println("   ", fn[varname], " -- ", (isfile(fn[varname]) ? "EXISTS" : "MISSING"))
end
println()

missing_value = convert(dtype, ncgetatt(fn["tos"], "tos", "_FillValue"))

function missing2nan!(x)
    x[x .== missing_value] .= NaN
end

function nan2missing!(x)
    x[isnan.(x)] .= missing_value
end


function readModelVar(varname, range_tuple=())
    local var

    @printf("# Reading var [%s]...", varname)

    if range_tuple == ()
        var = convert(Array{dtype}, ncread(fn[varname], varname))
    else
        var = convert(Array{dtype}, ncread(fn[varname], varname)[range_tuple...])
    end

    var[var .== missing_value] .= NaN

    println("done.")

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


mon_secs = 365.0 / 12.0 * 86400.0
Δt  = 1.0 * mon_secs

lon = ncread(fn["tos"], "lon")
lat = ncread(fn["tos"], "lat")

lon_len = length(lon)
lat_len = length(lat)


months = length(ncread(fn["tos"], "time"))
if months % 12 != 0
    error("There are $months months, not a multiple of 12")
end

years = Int(months / 12)
N     = (years-2) * 12 # Discard the first and final year
beg_t = 13


println("# Data info")
@printf("  Model  : %s\n", model_name)
@printf("  Length : %02d years.\n", years)
println("  Data type: ", dtype)
println()

println("# Constants")
@printf("ρ   : %.2f kg / m^3\n", ρ)
@printf("c_p : %.2f J / kg / K\n", c_p)
@printf("Δt  : %.2f s\n", Δt)
println()
