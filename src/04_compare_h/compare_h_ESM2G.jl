include("config.jl")
include("../lib/NetCDFHelper.jl")

using NetCDF

fn["fithQ"] = joinpath(data_path, "case3_hQ.jl.nc")
output_fn = joinpath(data_path, @sprintf("%s.nc", basename(@__FILE__)) )

h_fit  = ncread(fn["fithQ"],  "h")[:, :, :]
omlmax = ncread(fn["omlmax"], "omlmax")[:, :, :]

missing_value = ncgetatt(fn["fithQ"], "h", "missing_value")

h_fit[h_fit .== ncgetatt(fn["fithQ"], "h", "missing_value")] = NaN
omlmax[omlmax .== ncgetatt(fn["omlmax"], "omlmax", "missing_value")] = NaN

h_fit  = convert(Array{Float32}, h_fit)
omlmax = convert(Array{Float32}, omlmax)


h_model = Array{eltype(omlmax), 3}(size(h_fit)...)

for i = 1:12
    h_model[:,:,i] = mean(omlmax[:, :, i:end:12], 3)
end

NetCDFHelper.specialCopyNCFile(fn["omlmax"], output_fn, ["lat", "lon", "lat_vertices", "lon_vertices"])

deviation = h_fit - h_model

println(eltype(h_fit))
println(eltype(h_model))
println(eltype(deviation))
println(typeof(missing_value))

h_fit[isnan.(h_fit)] = missing_value
omlmax[isnan.(omlmax)] = missing_value
deviation[isnan.(deviation)] = missing_value

time = collect(Float32, 1:12)

for obj in [
    [
        h_fit, "h_fit", Dict(
            "long_name"=>"Fitted Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        h_model, "h_model", Dict(
            "long_name" => "Model estimated monthly maximum MLD thickness",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ], [
        deviation, "deviation", Dict(
            "long_name" => "Difference between fitted and model estimated MLD thickness",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ]

]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    @printf("Dumping variable: %s\n", varname)

    nccreate(
        output_fn,
        varname,
        "rlon",
        "rlat",
        "time", time,
        atts=varatts
    )
    
    ncwrite(var, output_fn, varname)

end

ncclose(output_fn)





#=

using PyPlot

deviation = h_fit - omlmax


=#
