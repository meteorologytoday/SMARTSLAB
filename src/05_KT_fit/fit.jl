include("config.jl")
include("NetCDFHelper.jl")

using NetCDF
using LinearAlgebra

function extend(a::AbstractArray, len::Int)
    if length(a) < len
        a = repeat(a, outer=ceil(Int, len / length(a)))
    end

    return a[1:len]
end

function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end



N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(h)
output_Td_star = 273.0 * ρ * c_p

dt2 = 2.0 * dt

# Gradient Descent method parameters
iter_max = 1000000
#update_bounds = [0.16, 0.16]
ϵ_converge_ratio_threshold = 5e-3
converge_count_target = 10
ϵ_mem = Inf
η = 10.0

@inline mod12(n) = mod(n-1, 12) + 1

# Assign h and Q initial condition
output_h *= 0.0
output_h += 30.0

# In the length of N
h = zeros(dtype, N)
Q = copy(h)
∂h∂t = copy(h)


# Construct delta function
I12 = zeros(dtype, 12, 12) + I
δ    = repeat(I12, outer=(1, nyrs-2))
∂δ∂t = repeat(
    (circshift(I12,(0, 1)) - circshift(I12,(0, -1))) / dt2,
    outer=(nyrs-2,1)
)

# For each grid point
for i = 1:length(rlons), j = 1:length(rlats)
    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
    end

    ∂T_star∂t = (T_star[i, j, beg_t + 1 : beg_t + 1 + (N - 1)]
               - T_star[i, j, beg_t - 1 : beg_t - 1 + (N - 1)]) / dt2

   
    h[1:12]  = output_h[i, j]
    Q[1:12] = output_Q[i, j]

    S_term  = S[i, j, beg_t : beg_t + (N-1)]
    B_term  = B[i, j, beg_t : beg_t + (N-1)]
    ΔT_star = T_star[i, j] .- Td_star

    # Gradient Descent Init
    converge_count = 0

    # For each iteration
    for k = 1 : iter_max

        repeat_fill!(h,    h[1:12])
        repeat_fill!(Q,   Q[1:12])
        repeat_fill!(∂h∂t, (h[14:25] - h[12:23]) / dt2)  # Assume the length is long enough

        Λ = convert(Array{dtype}, ∂h∂t .> 0.0)
        
        # Calculate ϵ and ϵ^2
        ϵ =  h .* ∂T_star∂t + ΔT_star .* ∂h∂t .* Λ .- S_term .- B_term .- Q
        ϵ2_sum = ϵ' * ϵ

        # Determine if iteration should stop or not
        Δϵ_ratio = (ϵ2_sum^0.5 - ϵ_mem) / ϵ_mem


        if Δϵ_ratio < 0.0 && abs(Δϵ) < ϵ_converge_ratio_threshold  # if it is converging
            converge_count +=1
        else # if it is diverging
            converge_count = 0
        end

        if converge_count >= converge_count_threshold
            output_h[i, j, :] = h[:]
            output_Q[i, j, :] = Q[:]
            break
        end
       
        # Calculate ∂Post∂h, ∂Post∂Q
        # Assume flat prior for now
        ∂Post∂h  = - (δ * (ϵ .* ∂T_star∂t) + ∂δ∂t * (Λ .* ΔT_star))
        ∂Post∂Q = δ * ϵ

        # Update h and Q
        h += ∂Post∂h .* η
        Q += ∂Post∂Q .* η 

    end
        
end










mask = isnan.(β)

β[mask] = missing_value
β_std[mask] = missing_value

dh_dt[isnan.(dh_dt)] = missing_value

time = collect(Float64, 1:12)

filename = @sprintf("%s.nc", basename(@__FILE__))
filename = joinpath(data_path, filename)

NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


for obj in [
    [
        β[:, :,  1:12], "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        β[:, :, 13:24], "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ], [
        β_std[:, :,  1:12], "h_std", Dict(
            "long_name"=>"Mixed-layer Thickness Standard Deviation",
            "units"=>"m",
            "missing_value" => missing_value
        )
    ], [
        β_std[:, :, 13:24], "Q_std", Dict(
            "long_name"=>"Q-flux Standard Deviation",
            "units"=>"W / m^2",
            "missing_value" => missing_value
        )
    ], [
        dh_dt, "dh_dt", Dict(
            "long_name" => "Mixed-layer Thickness Changing Ratei",
            "units"=>"m / s",
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
