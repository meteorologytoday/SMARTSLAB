include("config.jl")
include("NetCDFHelper.jl")

function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%d ", A[i,j]) 

        end
        @printf("\n")
    end
end



using NetCDF

T_star = SST * ρ * c_p

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

ϕ = Base.SparseArrays.spzeros(eltype(T_star), N, 12 * 2)   # for h and Q
ϕ = zeros(eltype(T_star), N, 12 * 2)   # for h and Q

v     = zeros(eltype(T_star), length(rlons), length(rlats), 24)
v_std = copy(v)

dt2 = 2.0 * dt
@inline mod12(n) = mod(n-1, 12) + 1
for i = 1:length(rlons), j = 1:length(rlats)

    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end
    
    if mask[i,j]
        v[i,j,:] = NaN
        continue
    end

    ϕ[:,:] = 0.0

    for t = 1:N

        # ϕ_h
        #@printf("t = %d, center: %d\n", t, beg_t + t - 1)
        #@printf("%d, %d\n", t, t_mod_12 + 1)
        ϕ[t, mod12(t+1)] =   T_star[i, j, (beg_t + t - 1) + 1] / dt2
        ϕ[t, mod12(t-1)] = - T_star[i, j, (beg_t + t - 1) - 1] / dt2

        # ϕ_Q
        ϕ[t, 12 + mod12(t)] =  - 1.0

    end

    #prtArr(ϕ[1:36, :])
    #exit()

    F =  TOT_F[i, j, beg_t:(beg_t + N - 1)]

    # Solve normal equation
    # ϕ v = F => v = ϕ \ F
    v[i, j, :] = ϕ \ F
    # Estimate standard deviation
    ϵ = F - ϕ * v[i, j, :]
    var = inv(ϕ'*ϕ) * (ϵ' * ϵ) / (1 + N)
    for k = 1:24
        v_std[i, j, k] = sqrt(var[k, k])
    end
        
end

mask = isnan.(v)

v[mask] = missing_value
v_std[mask] = missing_value

time = collect(Float64, 1:12)

filename = @sprintf("data/%s.nc", basename(@__FILE__))

NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


for obj in [
    [
        v[:, :,  1:12], "h", Dict(
            "long_name"=>"Constant Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        v[:, :, 13:24], "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ], [
        v_std[:, :,  1:12], "h_std", Dict(
            "long_name"=>"Constant Mixed-layer Thickness Standard Deviation",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        v_std[:, :, 13:24], "Q_std", Dict(
            "long_name"=>"Q-flux Standard Deviation",
            "units"=>"W / m^2",
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
        "time", time,
        atts=varatts
    )
    ncwrite(var, filename, varname)

end

ncclose(filename)

@printf("Output file: %s\n", filename)
