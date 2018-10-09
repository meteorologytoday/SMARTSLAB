include("config.jl")
include("NetCDFHelper.jl")

using NetCDF
using LinearAlgebra
NC_VERBOSE = true
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

Td_star = ( 273.15 + 4.0 ) * ρ * c_p

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(output_h)
output_converge = zeros(dtype, length(rlons), length(rlats))

#output_Td_star = 273.0 * ρ * c_p

dt2 = 2.0 * dt

# Gradient Descent method parameters
iter_max = 10
#update_bounds = [0.16, 0.16]
ϵ_converge_ratio_threshold = 5e-3
converge_count_threshold = 10
ϵ_mem = Inf
η = 1e-5

@inline mod12(n) = mod(n-1, 12) + 1

# Assign h and Q initial condition
output_h .= 30.0

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
    if j == 1
        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
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
    
    # For each iteration
    for k = 1 : iter_max
        global h, Q, ϵ_mem
        #@printf("Iterating: %d. ", k)
        repeat_fill!(h,     h[1:12])
        repeat_fill!(Q,     Q[1:12])
        repeat_fill!(∂h∂t, (h[14:25] - h[12:23]) / dt2)  # Assume the length is long enough

        #print(h)
        #@printf("Repeat_fill done.\n")
        Λ = convert(Array{dtype}, ∂h∂t .> 0.0)
       
         
        #@printf("Doing terrible matrix multiplication.\n")

        #=
        println(size(h))
        println(size(∂T_star∂t))

        println(size(ΔT_star))
        println(size(∂h∂t))
        println(size(Λ))
        
        println(size(S_term))
        println(size(B_term))
        println(size(Q))
        =#


        # Calculate ϵ and ϵ^2
        ϵ =  h .* ∂T_star∂t + ΔT_star .* ∂h∂t .* Λ .- S_term .- B_term .- Q
        ϵ2_sum = ϵ' * ϵ
        ϵ_now = (ϵ2_sum / N)^0.5

        #@printf("Determine convergence condition.\n")
        # Determine if iteration should stop or not
        Δϵ_ratio = (ϵ_now - ϵ_mem) / ϵ_mem

        #@printf("ϵ_now: %f; Δϵ_ratio * 100: %f\n", ϵ_now, Δϵ_ratio * 100.0)

        #@printf("Converging?\n")
        if Δϵ_ratio < 0.0 && abs(Δϵ_ratio) < ϵ_converge_ratio_threshold  # if it is converging
            converge_count +=1
        else # if it is diverging
            converge_count = 0
        end

        #@printf("Converge count threshold reached?\n")
        if converge_count >= converge_count_threshold
            break
        end
       
        #@printf("Gradient descent.\n")
        # Calculate ∂Post∂h, ∂Post∂Q
        # Assume flat prior for now
        #=
        println(size(δ))
        println(size(ϵ))
        println(size(∂T_star∂t))
        println(size(∂δ∂t))
        println(size(Λ))
        println(size(ΔT_star))
        =#
        #∂Post∂h  = ϵ .* ∂T_star∂t
        #∂Post∂h  = ∂δ∂t * (Λ .* ΔT_star)
        
        ∂Post∂h = - (δ * (ϵ .* ∂T_star∂t) + ∂δ∂t * (ϵ .* Λ .* ΔT_star))
        ∂Post∂Q = δ * ϵ

        #@printf("Update h and Q.\n")
        # Update h and Q
        h[1:12] += ∂Post∂h * η
        Q[1:12] += ∂Post∂Q * 0

        ϵ_mem = ϵ_now 
    end


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
