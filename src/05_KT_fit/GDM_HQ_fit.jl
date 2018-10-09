include("config.jl")
include("NetCDFHelper.jl")
include("BacktrackingLineSearchStruct.jl")

using NetCDF
using LinearAlgebra

function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end

@inline eucLen    = x -> (sum(x.^2.0))^(0.5)
@inline normalize = x -> x / encLen(x)

function getϵandϵ2sum(h, ∂h∂t, θ, ∂θ∂t, S, B, Q)
    ϵ = h .* ∂θ∂t + θ .* ∂h∂t .- S .- B .- Q
    return ϵ, ϵ' * ϵ
end

function GDM(
    iter_max,
    h_init,
    Q_init,
    ∂θ∂t,
    θ,
    S,
    B)
 
    global dt2, δ, ∂δ∂t, BLSS
  
    N = length(S)
    h = S * 0.0
    Q = coph(h)
    ∂h∂t = copy(h)

    dtype = eltype(S)
    ∇Post = zeros(dtype, 24)  # 1:12 = h1~h12  13:24 = Q1~Q12
 
    for k = 1 : iter_max

        repeat_fill!(h,     h[1:12])
        repeat_fill!(Q,     Q[1:12])
        repeat_fill!(∂h∂t, (h[14:25] - h[12:23]) / dt2)  # Assume the length is long enough

        ϵ, ϵ2sum = getϵandϵ2sum(h, ∂h∂t, θ, ∂θ∂t, S, B, Q)

        ∇Post[ 1:12] = - (δ * (ϵ .* ∂θ∂t) + ∂δ∂t * (ϵ .* θ)) / length(ϵ)
        ∇Post[13:24] = δ * ϵ / length(ϵ)

        ∇Post_len  = eucLen(∇Post)
        ∇Post_unit = normalize(∇Post)


        # 1: Test if this iteration should stop
        if ∇Post_len < BLSS.η
            @printf("The stop condition is met. End loop immediately.\n")
            break
        end

        # 2: Compare values of log(Posterior) of changing (h,Q) and Taylor extrapolation

        new_h[1:12] = h[1:12] + BLSS.t * ∇Post_unit[ 1:12]
        new_Q[1:12] = Q[1:12] + BLSS.t * ∇Post_unit[13:24]
 
        repeat_fill!(new_h, new_h[1:12])
        repeat_fill!(new_Q, new_Q[1:12])
        repeat_fill!(new_∂h∂t, (new_h[14:25] - new_h[12:23]) / dt2)
        
        new_ϵ_from_new_hQ, new_Post_from_new_hQ = getϵandϵ2sum(h, ∂h∂t, θ, ∂θ∂t, S, B, Q)
        new_Post_from_new_hQ *= -1.0

        new_Post_from_taylor = -ϵ2sum

        

         

        # Update h and Q
        h[1:12] += ∂Post∂h * η
        Q[1:12] += ∂Post∂Q * η

        ϵ_mem = ϵ_now 
    end






N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(output_h)
output_converge = zeros(dtype, length(rlons), length(rlats))

dt2 = 2.0 * dt

# Gradient Descent method parameters
iter_max = 10

BLSS = BacktrackingLineSearchStruct(
    N = iter_max,
    α = 0.1,
    β = 0.8,
    t = 10.0,
    η = 1e-3
)


# Assign h and Q initial condition
fn_HQ_REAL = joinpath(data_path, "case3_hQ.jl.nc")
output_h = ncread(fn_HQ_REAL, "h")
output_Q = ncread(fn_HQ_REAL, "Q")

# In the length of N
h = zeros(dtype, N)
Q = copy(h)
∂h∂t = copy(h)

# Construct delta function
I12 = zeros(dtype, 12, 12) + I
δ   = repeat(I12, outer=(1, nyrs-2))
∂δ∂t = repeat(
    (circshift(I12,(0, 1)) - circshift(I12,(0, -1))) / dt2,
    outer=(1, nyrs-2)
)


# For each grid point
for i = 1:length(rlons), j = 1:length(rlats)

#    if j == 1
#        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
#    end

    if i != 68 || j != 124
        continue
    end


    if missing_places[i, j, 1]
        continue
    end


    ∂T_star∂t = (T_star[i, j, beg_t + 1 : beg_t + 1 + (N - 1)]
               - T_star[i, j, beg_t - 1 : beg_t - 1 + (N - 1)]) / dt2

   
    h[1:12] = output_h[i, j, 1:12]
    Q[1:12] = output_Q[i, j, 1:12]


    #@printf("%d : %d\n", beg_t, beg_t + (N-1))
    rng = beg_t : beg_t + (N-1)
    S_term  = S[i, j, rng]
    B_term  = B[i, j, rng]
    ΔT_star = T_star[i, j, rng] .- Td_star

    # Gradient Descent Init
    converge_count = 0


    GDM(h, Q)
    output_h[i, j, :], output_Q[i, j, :] = GDM(h, Q, )


    
    output_h[i, j, :] = h[1:12]
    output_Q[i, j, :] = Q[1:12]

    if converge_count < converge_count_threshold
        output_converge[i, j] = -999.0
    end

end


missing_places_year   = missing_places[:, :, 1:12]
missing_places_single = missing_places[:, :, 1]

output_h[missing_places_year] .= missing_value 
output_Q[missing_places_year] .= missing_value 
output_converge[missing_places_single] .= missing_value

time = collect(Float32, 1:12)

filename = @sprintf("%s.nc", basename(@__FILE__))
filename = joinpath(data_path, filename)



NetCDFHelper.specialCopyNCFile(fn["tos"], filename, ["lat", "lon", "lat_vertices", "lon_vertices"])

println("Missing Value: ", missing_value)

println("Output converge count")
nccreate(
    filename,
    "converge",
    "rlon", rlons,
    "rlat", rlats,
    atts=Dict(
        "long_name" => "Converge count",
        "units"=>"(Scalar)",
        "missing_value" => missing_value,
    )

)
ncwrite(output_converge, filename, "converge")


println("Output h and Q")
for obj in [
    [
        output_h, "h", Dict(
            "long_name"=>"Mixed-layer Thickness",
            "units"=>"m",
            "missing_value" => missing_value,
        )
    ], [
        output_Q, "Q", Dict(
            "long_name"=>"Q-flux",
            "units"=>"W / m^2",
            "missing_value" => missing_value,
        )
    ]

]

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
