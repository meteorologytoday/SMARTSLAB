include("config.jl")
include("NetCDFHelper.jl")
function nansum(A::Array)
    s, n = 0.0, 0.0
    for val in A
        if !isnan(val)
            s += val
            n += 1.0
        end
    end
    return s / n
end

# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 

N = Int((nyrs-2) * 12)



dT_star_dt = (T_star[:, :, 14:end-11] - T_star[:, :, 12:end-13]) / (2.0 * mon_secs)
h = Array{tp, 2}(length(rlons), length(rlats))
h_std = copy(h)

ϕ = Array{tp, 2}(N, 2)

ϕ[:, 2] = 0.0

for i = 1:length(rlons), j = 1:length(rlats)

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end
    
    if spatial_mask[i,j]
        β[i, j, :] = NaN
        continue
    end

    ϕ[:, 1] = T_star_dt[i, j]
    
    # Solve normal equation
    # ϕ v = F => v = ϕ \ F
    β = ϕ \ F
    h[i, j]     = β[1]

    # Estimate standard deviation
    ϵ = F - ϕ * β[i, j, :]
    var = inv(ϕ'*ϕ) * (ϵ' * ϵ) / (1 + N)
    for k = 1:24
        β_std[i, j, k] = sqrt(var[k, k])
    end
 

end



print("Doing calculation... ")
#h = (sum(TOT_F .* dT_star_dt, 3) ./ sum(dT_star_dt .^ 2.0, 3))[:,:,1]





println("done.")


h[spatial_mask] = missing_value



filename = @sprintf("data/%s.nc", basename(@__FILE__))
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


for obj in [
    [
        h, "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
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
