include("../01_config/general_config.jl")
include("config.jl")

include("../lib/NetCDFHelper.jl")
include("./core/LR_flexible_slab_core.jl")


β     = zeros(dtype, lon_len, lat_len, 24)
β_std = copy(β)

for i = 1:lon_len, j = 1:lat_len

    global β, ϕ

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, lon[i])
    end

    if isnan(θ[i, j, 1])
        β[i,j,:] .= NaN
        continue
    end
    
    h, Q = LR_shallow_water!(
        pts_per_year = 12,
        F            = F[i, j, :],
        θ            = θ[i, j, :],
        Δt           = Δt,
    )

    β[i, j,  1:12] = h
    β[i, j, 13:24] = Q
    

end

nan2missing!(β)

time = collect(Float64, 1:12)

filename = @sprintf("%s.nc", basename(@__FILE__))
filename = joinpath(data_path, filename)

rm(filename)
for obj in [
    [
        β[:, :,  1:12], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "_FillValue" => missing_value,
            "missing_value" => missing_value,

        )
    ], [
        β[:, :, 13:24], "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "_FillValue" => missing_value,
            "missing_value" => missing_value
        )
    ]
]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    println("Writing ", varname, " with size: ", size(var))    

    nccreate(
        filename,
        varname,
        "lon", lon,
        "lat", lat,
        "time", 12,
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)

@printf("Output file: %s\n", filename)
