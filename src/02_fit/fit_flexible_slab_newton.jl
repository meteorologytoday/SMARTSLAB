include("../01_config/general_config.jl")

include("../lib/NetCDFHelper.jl")
include("./core/LR_flexible_slab_core.jl")


θ = readModelVar("tos")
F = readModelVar("hfds") * ρ * c_p

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

println("β[156, 110, :] = ", β[156, 110,:])

# output data
filename = format("{}-{}.nc", model_name, basename(@__FILE__))
filename = joinpath(data_path, filename)

ds = Dataset(filename, "c")

println("Creating dimension...")
defDim(ds, "time", 12)
defDim(ds, "lat", length(lat))
defDim(ds, "lon", length(lon))

defVar(ds, "time", dtype, ("time",))[:] = collect(1:12)
defVar(ds, "lat",  dtype, ("lat",))[:]  = lat
defVar(ds, "lon",  dtype, ("lon",))[:]  = lon
 
println("done")
for o in (
    [
        "h", β[:, :, 1:12], Dict(
        "long_name"=>"Mixed-layer Depth",
        "units"=>"m",
        )
    ], [
        "Q", β[:, :,13:24], Dict(
        "long_name"=>"Q-flux",
        "units"=>"W / m^2",
        )
    ],
)
    varname, vardata, varatts = o
    println("Writing ", varname, " with size: ", size(vardata))

    ncvar = defVar(ds, varname, eltype(vardata), ("lon", "lat", "time"))
    ncvar.attrib["_FillValue"] = missing_value
    for key in keys(varatts)
        ncvar.attrib[key] = varatts[key]
    end

    ncvar[:] = vardata
    println("done.")
end

close(ds)
println("Output file: ", filename)
