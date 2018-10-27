include("config.jl")
include("NetCDFHelper.jl")
include("BacktrackingLineSearchStruct.jl")

using NetCDF
using LinearAlgebra

Q_scale = 1.0

function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)
@inline mod12(n) = mod(n-1, 12) + 1





function getϵandϵ2avg(h, θ_p1, θ_m1, S, B, Q)
    global dt2
    h_p1 = circshift(h, -1)
    h_m1 = circshift(h,  1)
    ϵ = (h_p1 .* θ_p1 - h_m1 .* θ_m1) / dt2 - S - B - Q * Q_scale

    #= 
    ϵ = zeros(dtype, length(S))
    for i = 1:length(S)
        ϵ[i] = (h[mod12(i+1)] * θ_p1[i] - h[mod12(i-1)] * θ_m1[i]) / dt2 - S[i] - B[i] - Q[i]
        #ϵ[i] = (h_p1[i] * θ_p1[i] - h_m1[i] * θ_m1[i]) / dt2 - S[i] - B[i] - Q[i]
    end
    =#

    return ϵ, transpose(ϵ) * ϵ / length(ϵ)
end

function GDM(;
    h_init,
    Q_init,
    θ_p1,
    θ_m1,
    S,
    B)
 
    global dt2, δ, δ_p1, δ_m1, BLSS
 
    N = length(S)
    h = S * 0.0
    Q = copy(h)

    new_h = copy(h)
    new_Q = copy(h)

    dtype = eltype(S)
    ∇Post = zeros(dtype, 24)  # 1:12 = h1~h12  13:24 = Q1~Q12

    iter_count = 0

    h[1:12] = h_init
    Q[1:12] = Q_init

    for k = 1 : BLSS.N
        iter_count = k

        if_update = false

        repeat_fill!(h,     h[1:12])
        repeat_fill!(Q,     Q[1:12])
        ϵ, ϵ2avg = getϵandϵ2avg(h, θ_p1, θ_m1, S, B, Q)
        Post = - ϵ2avg

        ∇Post[ 1:12] = - (δ_p1 * (ϵ .* θ_p1) - δ_m1 * (ϵ .* θ_m1)) / dt2
        ∇Post[13:24] = δ * ϵ * Q_scale
        #∇Post[ 1:12] .= 0
        ∇Post /= length(ϵ)

        
        # Complex derivative check

        complex_check = ∇Post * 0.0

        η = 1e-50
        for d = 1:12
            h_tmp = h * 1.0
            h_tmp[d] += η * im 
            repeat_fill!(h_tmp, h_tmp[1:12])
            _, new_ϵ2avg = getϵandϵ2avg(h_tmp, θ_p1, θ_m1, S, B, Q)
            complex_check[d] = - imag(new_ϵ2avg) / η
        end

        for d = 1:12
            Q_tmp = Q * 1.0
            Q_tmp[d] += η * im 
            repeat_fill!(Q_tmp, Q_tmp[1:12])
            _, new_ϵ2avg = getϵandϵ2avg(h, θ_p1, θ_m1, S, B, Q_tmp)
            println(new_ϵ2avg)
            complex_check[d+12] = - imag(new_ϵ2avg) / η
        end

        println("∇Post:")
        println(∇Post)
        println("Complex Check:")
        println(complex_check)

        println(∇Post ./ complex_check)

        exit()



        ∇Post_len  = eucLen(∇Post)
        ∇Post_unit = normalize(∇Post)

        # 1: Test if this iteration should stop
        if ∇Post_len < BLSS.η
            @printf("The stop condition is met. End loop immediately.\n")
            break
        end

        # 2: Compare values of log(Posterior) of changing (h,Q) and Taylor extrapolation
        ## Add the gradient to move to larger posterior
        new_h[1:12] = h[1:12] + BLSS.t * ∇Post_unit[ 1:12]
        new_Q[1:12] = Q[1:12] + BLSS.t * ∇Post_unit[13:24]
 
        repeat_fill!(new_h, new_h[1:12])
        repeat_fill!(new_Q, new_Q[1:12])
        new_ϵ_from_new_hQ, new_Post_from_new_hQ = getϵandϵ2avg(new_h, θ_p1, θ_m1, S, B, new_Q)
        new_Post_from_new_hQ *= -1.0 

        new_Post_from_taylor = Post + BLSS.α * BLSS.t * ∇Post_len


        if new_Post_from_new_hQ < new_Post_from_taylor
            BLSS.t *= BLSS.β
        else
            #BLSS.t /= BLSS.β
            h[1:12] = new_h[1:12]
            Q[1:12] = new_Q[1:12]
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
        #println(Q[1:12]) 
    end

    return h[1:12], Q[1:12], iter_count
end





N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(output_h)
output_converge = zeros(dtype, length(rlons), length(rlats))

dt2 = 2.0 * dt

# Construct δ* function
Δ = (i, j) -> (i == (mod(j-1, 12) + 1)) ? 1 : 0
δ    = zeros(dtype, 12, N)
δ_p1 = copy(δ)
δ_m1 = copy(δ)

for k=1:12, i=1:N
    δ[k, i]    = Δ(k, i  )
    δ_p1[k, i] = Δ(k, i+1)
    δ_m1[k, i] = Δ(k, i-1)
end


# Gradient Descent method parameters
BLSS = BacktrackingLineSearchStruct(
    N = 100000,
    α = 0.1,
    β = 0.8,
    t = 1000.0,
    η = 1e-3
)


# Assign h and Q initial condition
#fn_HQ_REAL = joinpath(data_path, "case3_hQ.jl.nc")
#output_h = ncread(fn_HQ_REAL, "h") * 0 .+ 30.0
#output_Q = ncread(fn_HQ_REAL, "Q") * 0

output_h = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(output_h)




# For each grid point
for i = 1:length(rlons), j = 1:length(rlats)

#    if j == 1
#        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
#    end

    if i != 125 || j != 137
        continue
    end


    if missing_places[i, j, 1]
        continue
    end


    h_init = output_h[i, j, 1:12]
    Q_init = output_Q[i, j, 1:12]

    rng = collect(beg_t : beg_t + (N-1))
    S_now    = S[i, j, rng]
    B_now    = B[i, j, rng]
    θ_p1_now    = θ[i, j, rng .+ 1]
    θ_m1_now    = θ[i, j, rng .- 1]

    output_h[i, j, :], output_Q[i, j, :], iter_count = GDM( 
        h_init   = h_init,
        Q_init   = Q_init,
        θ_p1     = θ_p1_now,
        θ_m1     = θ_m1_now,
        S        = S_now,
        B        = B_now
    )

end


println(rlons[125], "  ", rlats[137])
println("## h")
println(output_h[125, 137, :])
println("## Q")
println(output_Q[125, 137, :] * Q_scale)


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
