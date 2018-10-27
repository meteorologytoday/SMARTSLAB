include("config.jl")
include("NetCDFHelper.jl")

# Equation:  h ∂θ/∂t = F + Q
# h and Q are only functions of space. No time dependence.

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

F    = (F[:, :, beg_t:beg_t+N-1] + F[:, :, beg_t+1:beg_t+N]) / 2.0 
∂θ∂t = (θ[:, :, beg_t+1:beg_t+N] - θ[:, :, beg_t:beg_t+N-1]) / dt

β = zeros(dtype, length(rlons), length(rlats), 2)

ϕ = zeros(dtype, N, 2)

ϵ2sum = 0.0
ϵ2count = 0
println("Doing calculation... ")
for i = 1:length(rlons), j = 1:length(rlats)
    global ϵ2sum, ϵ2count

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end
    
    if missing_places[i,j,1]
        β[i, j, :] .= NaN
        continue
    end

    ϕ[:, 1] .= ∂θ∂t[i, j, :]
    ϕ[:, 2] .= -1.0
    
    # Solve normal equation
    # ϕ β = _F => β = ϕ \ _F
    _F = F[i, j, :]
    β[i, j, :] = ϕ \ _F

    ϵ = _F - ϕ * β[i, j, :]
    ϵ2sum += sum(ϵ' * ϵ)
    ϵ2count += length(ϵ)
 
end

Q2 = β[:,:,2].^2.0
Q2count = sum(isfinite.(Q2), dims=(1,))
Q2[isnan.(Q2)] .= 0.0
Q2sum = sum(Q2, dims=(1,))
q = ((Q2sum ./ Q2count).^.5)[1,:]

println("######")
println("Mean Q magnitude: ", q)


println("######")
h = β[:, :, 1]
h = h[isfinite.(h)]
unreal_h_count = sum(h .< 0)
@printf("Mean residue: %.e\n", (ϵ2sum / ϵ2count)^0.5)
@printf("Unrealistic h grid points: %d / %d \n", unreal_h_count, length(h))



println("done.")

h = β[:, :, 1]
h = h[isfinite.(h)]
unreal_h_count = sum(h .< 0)
@printf("Unrealistic h grid points: %d / %d \n", unreal_h_count, length(h))




filename = @sprintf("data/%s.nc", basename(@__FILE__))
NetCDFHelper.specialCopyNCFile(fn["tos"], filename, ["lat", "lon", "lat_vertices", "lon_vertices"])

for obj in [
    [
        β[:,:,1], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        β[:,:,2], "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ]
]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]
    
    var[isnan.(var)] .= missing_value
    nccreate(
        filename,
        varname,
        "rlon",
        "rlat",
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
