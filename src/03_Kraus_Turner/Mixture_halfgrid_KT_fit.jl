include("simulate_given_dhdt.jl")
include("../lib/Newton.jl")

using LinearAlgebra
using Statistics



N       = Int((cycles - 2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

Δt2 = 2.0 * Δt
println("Δt = ", Δt)


output_h = zeros(dtype, 12)
output_Q = copy(output_h)

η = 1e-5

# Construct δ* function
Δ = (i, j) -> (i == (mod(j-1, 12) + 1)) ? 1 : 0
δ    = zeros(dtype, 12, N)
δ_p1 = copy(δ)

for k=1:12, i=1:N
    δ[k, i]    = Δ(k, i  )
    δ_p1[k, i] = Δ(k, i+1)
end
γ = (δ_p1 - δ) / Δt



Λ_func   = (x, a, b) -> 1.0 + b / 2.0 * (tanh(x/a) - 1.0)
∂Λ_func  = (x, a, b) -> b / a * (sech(x/a)^2.0)  / 2.0
∂∂Λ_func = (x, a, b) -> - (b / a^2.0)  * (sech(x/a)^2.0) * tanh(x/a)

ab_pairs = [
    [1e-5, 1.00],
    [1e-11, 1.00],
]

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
    ∇ϵ[i, j+12] = - δ[j, i] * 1
end

function g_and_∇g(;h, Q_ph, θ, θ_p1, S_ph, B_ph, a, b)

    repeat_fill!(h, h[1:12])
    repeat_fill!(Q_ph, Q_ph[1:12])

    h_p1 = circshift(h, -1)

    ∂h∂t = (h_p1 - h) / Δt

    Λ   =   Λ_func.(∂h∂t, a, b)
    ∂Λ  =  ∂Λ_func.(∂h∂t, a, b)
    ∂∂Λ = ∂∂Λ_func.(∂h∂t, a, b)
    
    ϵ =  (
        h    / Δt2 .* ( θ_p1 .* (1.0 .- Λ) - θ .* (1.0 .+ Λ) )
        + h_p1 / Δt2 .* ( θ_p1 .* (1.0 .+ Λ) - θ .* (1.0 .- Λ) )
        - S_ph - B_ph - Q_ph
    )

    for i=1:N, j=1:12
        ∇ϵ[i, j] =  (
            δ[j, i]    / Δt2 * ( θ_p1[i] * (1.0 - Λ[i]) - θ[i] * (1.0 + Λ[i]) )
            + δ_p1[j, i] / Δt2 * ( θ_p1[i] * (1.0 + Λ[i]) - θ[i] * (1.0 - Λ[i]) )
            + ∂h∂t[i] * (θ[i] + θ_p1[i]) / 2.0 * γ[j, i] * ∂Λ[i]
        )
    end

    g  = ∇ϵ' * ϵ
    ∇g = ∇ϵ' * ∇ϵ

    # Add ϵ ∂^2ϵ/∂h^2 to ∇g
    #=
    for i=1:12, j=1:12
        ∇g[i, j] += sum(
            ϵ .* γ[j, :] .* γ[i, :] .* θ_p1 .* (2.0 * ∂Λ - ∂∂Λ .* ∂h∂t)
        )
    end
    =#
    return g, ∇g 
end


rng1 = collect(beg_t : beg_t + (N-1))
rng2 = rng1 .+ 1

h_long_vec = zeros(dtype, N)
Q_long_vec = zeros(dtype, N)
x_mem = zeros(dtype, 24)

let
    _S_ph = (S[rng1] + S[rng2]) / 2.0
    _B_ph = (B[rng1] + B[rng2]) / 2.0

    _θd   = _S_ph * 0.0
    _θd[1:12] .= 273.15
    repeat_fill!(_θd, _θd[1:12])

    _θs    = θs[rng1] - _θd
    _θs_p1 = θs[rng2] - _θd

    global x_mem
    x_mem[ 1:12] = output_h[:] 
    x_mem[13:24] = output_Q[:] 


    # Gradually change entrainment condition
    for i in 1:size(ab_pairs)[1]
        println("Now we are doing ab_pair: ", ab_pairs[i])
        println(x_mem)
        _g_and_∇g = function(x)
            h_long_vec[1:12] = x[ 1:12]
            Q_long_vec[1:12] = x[13:24]
            return g_and_∇g(
                h = h_long_vec,
                Q_ph = Q_long_vec,
                θ = _θs,
                θ_p1 = _θs_p1,
                S_ph = _S_ph,
                B_ph = _B_ph,
                a    = ab_pairs[i][1],
                b    = ab_pairs[i][2]
            )
        end


        x_mem[:] = Newton(
            g_and_∇g = _g_and_∇g,
            η        = η,
            x0       = x_mem,
            max      = 1000
        )

        #println((x_mem[2:12] - x_mem[1:11]) / Δt)
        println("Q std: ", std(x_mem[13:24]))
    end

    output_h[:] = x_mem[ 1:12]
    output_Q[:] = x_mem[13:24]

end


#println(rlons[125], "  ", rlats[137])
println("## h")
println(output_h)
println("## Q")
println(output_Q)
