include("config.jl")
include("NetCDFHelper.jl")

using NetCDF
using LinearAlgebra

# Equation:  ∂hθ/∂t = F + Q
# h and Q are functions of time and space

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

ϕ = zeros(dtype, N, 12 * 2)   # for h and Q

β     = zeros(dtype, length(rlons), length(rlats), 24)
β_std = copy(β)

dh_dt = zeros(dtype, length(rlons), length(rlats), 12)

println("dt2: ", dt2)


rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

ϵ2sum = 0.0
ϵ2count = 0
unreal_h_count = 0
for i = 1:length(rlons), j = 1:length(rlats)
    global ϵ2sum, ϵ2count, unreal_h_count
    if j == 1
#        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end

    if i != 68 || j != 124
        continue
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
    #β[i, j, :] = convert(Array{Float32}, ϕ) \ convert(Array{Float32}, _F)

   
 
    println("sum ϕ: ", sum(abs.(ϕ - convert(Array{Float32,2}, ϕ))))
    Φ = convert(Array{Float32,2}, ϕ)
    println("det( ϕTϕ ): ", det(ϕ' * ϕ))
    println("det( ΦTΦ ): ", det(Φ' * Φ))


    println(β[i, j, :])


    # Estimate standard deviation
    ϵ = _F - ϕ * β[i, j, :]
    ϵ2sum += sum(ϵ' * ϵ)
    ϵ2count += length(ϵ)
 
    if sum(β[i, j, 1:12] .< 0) != 0
        unreal_h_count += 1
    end

end

Q2 = β[:,:,13:24].^2.0
Q2count = sum(isfinite.(Q2), dims=(1,3))
Q2[isnan.(Q2)] .= 0.0
Q2sum = sum(Q2, dims=(1,3))
q = ((Q2sum ./ Q2count).^.5)[1,:,1]
println("length(q): ", length(q))

println("######")
println("Mean Q magnitude: ", q)


println("######")
@printf("Mean residue: %.e\n", (ϵ2sum / ϵ2count)^0.5)
@printf("Unrealistic h grid points: %d\n", unreal_h_count)

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


    var[isnan.(var)] .= missing_value

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

for obj in (
    [
        q, "q", Dict(
            "long_name" => "Zonal mean of q",
            "units"=>"W / m^2",
            "_FillValue" => missing_value,
            "missing_value" => missing_value
        )
    ],
)
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    var[isnan.(var)] .= missing_value
    println("Writing ", varname)

    nccreate(
        filename,
        varname,
        "rlat",
        atts=varatts
    )
    ncwrite(var, filename, varname)

end


ncclose(filename)

@printf("Output file: %s\n", filename)
