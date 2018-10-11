include("config.jl")
include("NetCDFHelper.jl")
include("BacktrackingLineSearchStruct.jl")

using NetCDF
using LinearAlgebra

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year
#rng = beg_t : beg_t + (N-1)

ϕ = zeros(dtype, N, 12)
@inline mod12(n) = mod(n-1, 12) + 1

β = zeros(dtype, length(rlons), length(rlats), 12)
println("dt2: ", dt2)
# For each grid point
for i = 1:length(rlons), j = 1:length(rlats)
    global ϕ, β

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end

    if i != 68 || j != 124
    #    continue
    end


    if missing_places[i, j, 1]
        continue
    end

    # Construct normal matrix
    ϕ .= 0.0
    for t = 1:N
        # ϕ_h
        ϕ[t, mod12(t+1)] =   θ[i, j, (beg_t + t - 1) + 1] / dt2
        ϕ[t, mod12(t-1)] = - θ[i, j, (beg_t + t - 1) - 1] / dt2

        # ϕ_Q
        #ϕ[t, 12 + mod12(t)] =  - 1.0
    end

    #prtArr(ϕ)

    F =  S[i, j, beg_t:(beg_t + N - 1)] + B[i, j, beg_t:(beg_t + N - 1)]
    
    #println(F)
    
    # Solve normal equation
    # ϕ β = F ⇒ β = ϕ \ F
    β[i, j, :] = ϕ \ F
    #println(β[i, j, :])
    #println(i, j)

end


println("## h")
println(β[68, 124, :])


β[isnan.(β)] .= missing_value 

time = collect(Float32, 1:12)

filename = @sprintf("%s.nc", basename(@__FILE__))
filename = joinpath(data_path, filename)

NetCDFHelper.specialCopyNCFile(fn["tos"], filename, ["lat", "lon", "lat_vertices", "lon_vertices"])

println("Missing Value: ", missing_value)

println("Output h")
for obj in (
    (
        β[:, :, 1:12], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value,
        )
    ),
)

    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]

    nccreate(
        filename,
        varname,
        "rlon", rlons, 
        "rlat", rlats,
        "time", 12,

        atts=varatts
    )
    ncwrite(var, filename, varname)

end


ncclose(filename)

@printf("Output file: %s\n", filename)
