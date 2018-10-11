include("config.jl")
include("NetCDFHelper.jl")

using NetCDF

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

ϕ = zeros(dtype, N, 12 * 2)   # for h and Q

β     = zeros(dtype, length(rlons), length(rlats), 24)
β_std = copy(β)

dh_dt = zeros(dtype, length(rlons), length(rlats), 12)

println("dt2: ", dt2)
@inline mod12(n) = mod(n-1, 12) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

F = S + B

for i = 1:length(rlons), j = 1:length(rlats)
    global β, ϕ
    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end

    if i != 68 || j != 124
    #    continue
    end


    if missing_places[i, j, 1]
        β[i,j,:] .= NaN
        continue
    end

    ϕ .= 0.0

    for t = 1:N

        # ϕ_h
        ϕ[t, mod12(t  )] = - θ[i, j, (beg_t + t - 1)    ] / dt
        ϕ[t, mod12(t+1)] =   θ[i, j, (beg_t + t - 1) + 1] / dt

        # ϕ_Q
        ϕ[t, 12 + mod12(t  )] =  - .5
        ϕ[t, 12 + mod12(t+1)] =  - .5

    end

    _F =  (F[i, j, rng1] + F[i, j, rng2]) / 2.0

    # Solve normal equation
    # ϕ β = F => β = ϕ \ F
    β[i, j, :] = ϕ \ _F

    #=
    # Estimate standard deviation
    ϵ = _F - ϕ * β[i, j, :]
    var = inv(ϕ'*ϕ) * (ϵ' * ϵ) / (1 + N)
    for k = 1:24
        β_std[i, j, k] = (var[k, k])^0.5
    end
    =#  
end

println("## h")
println(β[68, 124, :])


# Derive dh_dt
for i = 1:length(rlons), j = 1:length(rlats)
    if missing_places[i,j,1]
        dh_dt[i,j,:] .= NaN
        continue
    end

    for t = 1:12
        dh_dt[i, j, t] = (β[i, j, mod12(t+1)] - β[i, j, mod12(t-1)]) / dt2
    end
end


mask = isnan.(β)

β[mask] .= missing_value
β_std[mask] .= missing_value

dh_dt[isnan.(dh_dt)] .= missing_value

time = collect(Float64, 1:12)

filename = @sprintf("%s.nc", basename(@__FILE__))
filename = joinpath(data_path, filename)

NetCDFHelper.specialCopyNCFile(fn["tos"], filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


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
    ], [
        β_std[:, :,  1:12], "h_std", Dict(
            "long_name"=>"Mixed-layer Thickness Standard Deviation",
            "units"=>"m",
            "_FillValue" => missing_value,
            "missing_value" => missing_value
        )
    ], [
        β_std[:, :, 13:24], "Q_std", Dict(
            "long_name"=>"Q-flux Standard Deviation",
            "units"=>"W / m^2",
            "_FillValue" => missing_value,
            "missing_value" => missing_value
        )
    ], [
        dh_dt, "dh_dt", Dict(
            "long_name" => "Mixed-layer Thickness Changing Rate",
            "units"=>"m / s",
            "_FillValue" => missing_value,
            "missing_value" => missing_value
        )
    ]
]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    println("Writing ", varname)

    nccreate(
        filename,
        varname,
        "rlon",
        "rlat",
        "time", 12,
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)

@printf("Output file: %s\n", filename)
