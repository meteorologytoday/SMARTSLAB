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

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)

function getϵandϵ2avg(h, ∂h∂t, θ, ∂θ∂t, S, B)
    ϵ = h .* ∂θ∂t + θ .* ∂h∂t .- S .- B
    return ϵ, ϵ' * ϵ / length(ϵ)
end

function GDM(;
    h_init,
    ∂θ∂t,
    θ,
    S,
    B)
 
    global dt2, δ, ∂δ∂t, BLSS
 
    N = length(S)
    h = S * 0.0
    ∂h∂t = copy(h)

    new_h = copy(h)
    new_∂h∂t = copy(h)


    dtype = eltype(S)
    ∇Post = zeros(dtype, 12)  # 1:12 = h1~h12  13:24 = Q1~Q12

    iter_count = 0

    h[1:12] = h_init

    for k = 1 : BLSS.N

        iter_count = k

        if_update = false

        repeat_fill!(h,     h[1:12])
        repeat_fill!(∂h∂t, (h[14:25] - h[12:23]) / dt2)  # Assume the length is long enough

        ϵ, ϵ2avg = getϵandϵ2avg(h, ∂h∂t, θ, ∂θ∂t, S, B)

        Post = - ϵ2avg

        ∇Post[:] = - (δ * (ϵ .* ∂θ∂t) + ∂δ∂t * (ϵ .* θ))
 

        ∇Post /= length(ϵ)

        ∇Post_len  = eucLen(∇Post)
        ∇Post_unit = normalize(∇Post)

        #println(∇Post_unit)
        #println(∇Post_unit)

        # 1: Test if this iteration should stop
        if ∇Post_len < BLSS.η
            @printf("The stop condition is met. End loop immediately.\n")
            break
        end

        # 2: Compare values of log(Posterior) of changing (h,Q) and Taylor extrapolation

        # Add the gradient to move to larger posterior
        new_h[1:12] = h[1:12] + BLSS.t * ∇Post_unit
 
        repeat_fill!(new_h, new_h[1:12])
        repeat_fill!(new_∂h∂t, (new_h[14:25] - new_h[12:23]) / dt2)
        
        new_ϵ_from_new_hQ, new_Post_from_new_hQ = getϵandϵ2avg(new_h, new_∂h∂t, θ, ∂θ∂t, S, B)
        new_Post_from_new_hQ *= -1.0 

        new_Post_from_taylor = Post + BLSS.α * BLSS.t * ∇Post_len


        if new_Post_from_new_hQ < new_Post_from_taylor
            BLSS.t *= BLSS.β
        else
            h[1:12] = new_h[1:12]
            if_update = true
        end
    
        @printf("Iter: %d. Update? %s. Post: %f, ΔPost_from_new_hQ: %f, ΔPost_from_taylor: %f, BLSS.t: %.2e, ∇Post_len: %f\n",
            k,
            (if_update) ? "YES" : "NO",
            Post,
            new_Post_from_new_hQ - Post,
            new_Post_from_taylor - Post,
            BLSS.t,
            ∇Post_len
	) 

    end

    return h[1:12], iter_count
end





N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_converge = zeros(dtype, length(rlons), length(rlats))

dt2 = 2.0 * dt

# Gradient Descent method parameters
BLSS = BacktrackingLineSearchStruct(
    N = 20,
    α = 0.1,
    β = 0.5,
    t = 10.0,
    η = 1e-3
)


# Assign h and Q initial condition
fn_HQ_REAL = joinpath(data_path, "case3_hQ.jl.nc")
output_h = ncread(fn_HQ_REAL, "h") * 0 .+ 100
output_Q = ncread(fn_HQ_REAL, "Q") * 0

# Construct delta function
I12 = zeros(dtype, 12, 12) + I
δ   = repeat(I12, outer=(1, nyrs-2))
∂δ∂t = repeat(
    (circshift(I12,(0, 1)) - circshift(I12,(0, -1))) / dt2,
    outer=(1, nyrs-2)
)
println("## h")
println(output_h[68, 124, :])
println("## Q")
println(output_Q[68, 124, :])



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


    h_init = output_h[i, j, 1:12]

    rng = beg_t : beg_t + (N-1)
    S_now    = S[i, j, rng]
    B_now    = B[i, j, rng]
    θ_now    = T_star[i, j, rng]
    ∂θ∂t_now = (T_star[i, j, beg_t + 1 : beg_t + 1 + (N - 1)]
               - T_star[i, j, beg_t - 1 : beg_t - 1 + (N - 1)]) / dt2




    output_h[i, j, :], iter_count = GDM( 
        h_init   = h_init,
        ∂θ∂t     = ∂θ∂t_now,
        θ        = θ_now,
        S        = S_now,
        B        = B_now
    )

end


println("## h")
println(output_h[68, 124, :])
println("## Q")
println(output_Q[68, 124, :])


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
