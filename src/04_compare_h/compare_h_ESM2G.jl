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
h_model_std = copy(h_model)


for i = 1:12
    month_collect      = omlmax[:, :, i:12:end]
    h_model[:,:,i]     = mean(month_collect, 3)
    h_model_std[:,:,i] = std(month_collect, 3)
end



NetCDFHelper.specialCopyNCFile(fn["omlmax"], output_fn, ["lat", "lon", "lat_vertices", "lon_vertices"])

deviation = h_fit - h_model

h_fit[isnan.(h_fit)] = missing_value
h_model[isnan.(h_model)] = missing_value
h_model_std[isnan.(h_model_std)] = 
omlmax[isnan.(omlmax)] = missing_value
deviation[isnan.(deviation)] = missing_value


time = collect(Float32, 1:12)

for obj in [
    [
        h_fit, "h_fit", Dict{String, Any}(
            "long_name"=>"Fitted MLD",
            "units"=>"m"
        )
    ], [
        h_model, "h_model", Dict{String, Any}(
            "long_name" => "Model estimated monthly maximum MLD",
            "units"=>"m"
        )
    ], [
        deviation, "deviation", Dict{String, Any}(
            "long_name" => "Difference between fitted and model estimated MLD",
            "units"=>"m"
        )
    ], [
        h_model_std, "h_model_std", Dict{String, Any}(
            "long_name" => "Standard deviation of model estimated monthly maximum MLD",
            "units"=>"m"
        )

    ]

]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    varatts["_FillValue"] = missing_value
    varatts["missing_value"] = missing_value

    @printf("Dumping variable: %s whose type is %s\n", varname, typeof(var))

    nccreate(
        output_fn,
        varname,
        "rlon",
        "rlat",
        "time", time,
        atts=varatts,
        t=NC_FLOAT
    )
    
    ncwrite(var, output_fn, varname)

end

ncclose(output_fn)




