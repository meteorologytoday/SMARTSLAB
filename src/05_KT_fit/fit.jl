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
output_Qf = copy(h)
output_Td_star = 273.0 * ρ * c_p

dt2 = 2.0 * dt
@inline mod12(n) = mod(n-1, 12) + 1

# Assign h and Qf initial condition
output_h *= 0.0
output_h += 30.0

# In the length of N
h = zeros(dtype, N)
Qf = copy(h)
∂h∂t = copy(h)


#
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
    Qf[1:12] = output_Qf[i, j]

    S_term  = S[i, j, beg_t : beg_t + (N-1)]
    B_term  = B[i, j, beg_t : beg_t + (N-1)]
    ΔT_star = T_star[i, j] .- Td_star
    # For each iteration
    for k = 1 : iter_max

        repeat_fill!(h,    h[1:12])
        repeat_fill!(Qf,   Qf[1:12])
        repeat_fill!(∂h∂t, (h[14:25] - h[12:23]) / dt2)  # Assume the length is long enough

        Λ = convert(Array{dtype}, ∂h∂t .> 0.0)
        
        # Calculate ϵ
        ϵ =  h .* ∂T_star∂t + ΔT_star .* ∂h∂t .* Λ .- S_term .- B_term .- Qf
        ϵ2 = ϵ' * ϵ
       
        # Calculate ∂Post∂h, ∂Post∂Qf
        # Assume flat prior for now
        ∂Post∂h  = - δ * (ϵ .* ∂T_star∂t) + ∂δ∂t * (Λ .* (ΔT_star))
        ∂Post∂Qf = δ * ϵ


        # Update h and Qf
        
        # Calculate ϵ^2

        # Test if ϵ^2 decreases and converges
        # If it does then do the next grid point
    end
        
end

# Derive dh_dt
for i = 1:length(rlons), j = 1:length(rlats)
    if spatial_mask[i,j]
        dh_dt[i,j,:] = NaN
        continue
    end

    for t = 1:12
        dh_dt[i, j, t] = (β[i, j, mod12(t+1)] - β[i, j, mod12(t-1)]) / dt2
    end
end





function Λ(x)
    return (x >= 0) ? 1.0 : 0.0
end





   
        


function cal_ϵ(h, Ts,)
    global dt2
    

end

function cal_∂P∂h()

end

function cal_∂P∂Qf()
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
