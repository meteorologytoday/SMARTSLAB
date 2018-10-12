include("config.jl")
include("NetCDFHelper.jl")

assumption_break_value = -9999.0
singular_matrix_value = -999.0

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

TOT_F = TOT_F[:, :, beg_t:beg_t+N-1]  # Discard the first and last year
dT_star_dt = (T_star[:, :, beg_t+1:beg_t+N] - T_star[:, :, beg_t-1:beg_t+N-2]) / (2.0 * mon_secs)

β = zeros(tp, length(rlons), length(rlats), 1)
β_std = copy(β)

ϕ = zeros(tp, N, 1)



println("Doing calculation... ")

# Equation:  h ∂θ/∂t = F  
# h is only function of space. No time dependence.


ϵ2sum = 0.0
count = 0
for i = 1:length(rlons), j = 1:length(rlats)
    global ϵ2sum, count

    local var

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end
    
    if spatial_mask[i,j]

        β[i, j, :] .= NaN
        β_std[i, j, :] .= NaN
        continue
    end

    ϕ[:, 1] = dT_star_dt[i, j, :]
    
    # Solve normal equation
    # ϕ β = F => β = ϕ \ TOT_F
    β[i, j, :] = ϕ \ TOT_F[i, j, :]

    # Estimate standard deviation
    ϵ = TOT_F[i, j, :] - ϕ * β[i, j, :]
    ϵ2sum += sum(ϵ' * ϵ)
    count += length(ϵ)
    try
        var = inv(ϕ'*ϕ) * (ϵ' * ϵ) / (1 + N)
    catch e
        #println("Exception happened")
        #if isa(e, Base.LinAlg.SingularException)
            #println("Singular detected")
            β_std[i, j, :] = singular_matrix_value
            continue
        #else
        #    println(e)
        #    throw(e)
        #end
    end
    
    for k = 1:1
        β_std[i, j, k] = var[k, k] < 0 ? assumption_break_value : sqrt(var[k, k])
    end

end
println("done.")

h = β[:, :, 1]
h = h[isfinite.(h)]
unreal_h_count = sum(h .< 0)
@printf("Total Residue: %.e\n", ϵ2sum)
@printf("Valid count: %d\n", count)
@printf("Mean residue: %.e\n", (ϵ2sum / count)^0.5)
@printf("Unrealistic h grid points: %d / %d \n", unreal_h_count, length(h))




β[isnan.(β)] .= missing_value
β_std[isnan.(β_std)] .= missing_value

filename = @sprintf("data/%s.nc", basename(@__FILE__))
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])

for obj in [
    [
        β[:,:,1], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        β_std[:,:,1], "h_std", Dict(
            "long_name"=>"Mixed-layer Thickness Standard Deviation",
            "units"=>"m",
            "missing_value" => missing_value,
            "assumption_break_value" => assumption_break_value,
            "singular_matrix_value" => singular_matrix_value
        )
    ]


]
    var     = obj[1]
    varname = obj[2]
    varatts = obj[3]
    nccreate(
        filename,
        varname,
        "rlon",
        "rlat",
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)
