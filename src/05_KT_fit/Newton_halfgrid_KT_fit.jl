include("config.jl")
include("NetCDFHelper.jl")
include("Newton.jl")

using NetCDF
using LinearAlgebra

iii = 168
jjj = 46

iii, jjj = 125, 137

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

output_h  = zeros(dtype, length(rlons), length(rlats), 12)
output_Q = copy(output_h)
output_converge = zeros(dtype, length(rlons), length(rlats))

dt2 = 2.0 * dt

η = 1e-5

# Construct δ* function
Δ = (i, j) -> (i == (mod(j-1, 12) + 1)) ? 1 : 0
δ    = zeros(dtype, 12, N)
δ_p1 = copy(δ)

for k=1:12, i=1:N
    δ[k, i]    = Δ(k, i  )
    δ_p1[k, i] = Δ(k, i+1)
end
γ = (δ_p1 - δ) / dt

α = 1e6
β = 1.0

Λ_func   = x -> 1.0 + β / 2.0 * (tanh(α * x) - 1.0)
∂Λ_func  = x -> α * β * (sech(α*x)^2.0)  / 2.0
∂∂Λ_func = x -> - α^2.0 * β * (sech(α*x)^2.0) * tanh(α*x)




function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)
@inline mod12(n) = mod(n-1, 12) + 1



∇ϵ = zeros(dtype, N, 24)
for i=1:N, j=1:12
    ∇ϵ[i, j+12] = - δ[j, i]
end

#∇∇ϵ = zeros(dtype, N, 24, 24)
#for i=1:N, j=1:24, k=1:24
#    if 
#    ∇∇ϵ[i, j, k] = - δ[j, i]
#end


function g_and_∇g(;h, Q_ph, θ, θ_p1, S_ph, B_ph)

    repeat_fill!(h, h[1:12])
    repeat_fill!(Q_ph, Q_ph[1:12])

    h_p1 = circshift(h, -1)

    ∂h∂t = (h_p1 - h) / dt

    Λ   =   Λ_func.(∂h∂t)
    ∂Λ  =  ∂Λ_func.(∂h∂t)
    ∂∂Λ = ∂∂Λ_func.(∂h∂t)
    
    #println(Λ)   
    #println(∂Λ)   
    #println(∂∂Λ)   
 
    ϵ =  (
        h    / dt2 .* ( θ_p1 .* (1.0 .- Λ) - θ .* (1.0 .+ Λ) )
        + h_p1 / dt2 .* ( θ_p1 .* (1.0 .+ Λ) - θ .* (1.0 .- Λ) )
        - S_ph - B_ph - Q_ph
    )

    for i=1:N, j=1:12
        ∇ϵ[i, j] =  (
            δ[j, i]    / dt2 * ( θ_p1[i] * (1.0 - Λ[i]) - θ[i] * (1.0 + Λ[i]) )
            + δ_p1[j, i] / dt2 * ( θ_p1[i] * (1.0 + Λ[i]) - θ[i] * (1.0 - Λ[i]) )
            + ∂h∂t[i] * (θ[i] + θ_p1[i]) / 2.0 * γ[j, i] * ∂Λ[i]
        )
    end

    g  = ∇ϵ' * ϵ
    ∇g = ∇ϵ' * ∇ϵ

    # Add ϵ ∂^2ϵ/∂h^2 to ∇g
    for i=1:12, j=1:12
        ∇g[i, j] += sum(
            ϵ .* γ[j, :] .* γ[i, :] .* θ_p1 .* (2.0 * ∂Λ - ∂∂Λ .* ∂h∂t)
        )
    end

    return g, ∇g 
end


# Assign h and Q initial condition
fn_HQ_REAL = joinpath(data_path, "LR_halfgrid_hQ_fit.jl.nc")
output_h = ncread(fn_HQ_REAL, "h") 
output_Q = ncread(fn_HQ_REAL, "Q_Delta")

println("## Linear Regression h")
println(output_h[iii, jjj, :])
println("## Linear Regression Q")
println(output_Q[iii, jjj, :] )

#output_h .*= rand(size(output_h)...) * 0.0 
#output_Q .*= rand(size(output_h)...) * 0.0


rng1 = collect(beg_t : beg_t + (N-1))
rng2 = rng1 .+ 1

h_long_vec = zeros(dtype, N)
Q_long_vec = zeros(dtype, N)
x0 = zeros(dtype, 24)
for i = 1:length(rlons), j = 1:length(rlats)

#    if j == 1
#        @printf("Doing lon[%d]: %.1f\n", i, rlons[i])
#    end

    if i != iii || j != jjj
        continue
    end


    if missing_places[i, j, 1]
        continue
    end

    _S_ph = (S[i, j, rng1] + S[i, j, rng2]) / 2.0
    _B_ph = (B[i, j, rng1] + B[i, j, rng2]) / 2.0
    _θ    = θ[i, j, rng1]
    _θ_p1 = θ[i, j, rng2]

    _g_and_∇g = function(x)
        h_long_vec[1:12] = x[ 1:12]
        Q_long_vec[1:12] = x[13:24]
        return g_and_∇g(
            h = h_long_vec,
            Q_ph = Q_long_vec,
            θ = _θ,
            θ_p1 = _θ_p1,
            S_ph = _S_ph,
            B_ph = _B_ph
        )
    end

    global x0
    x0[ 1:12] = output_h[i, j, :] 
    x0[13:24] = output_Q[i, j, :] 

    x_opt = Newton(
        g_and_∇g = _g_and_∇g,
        η        = η,
        x0       = x0,
        max      = 1000
    )

    output_h[i, j, :] = x_opt[ 1:12]
    output_Q[i, j, :] = x_opt[13:24]
end


#println(rlons[125], "  ", rlats[137])
println("## h")
println(output_h[iii, jjj, :])
println("## Q")
println(output_Q[iii, jjj, :] )



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
