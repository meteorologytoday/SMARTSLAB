include("../../lib/Newton.jl")

using NetCDF
using LinearAlgebra
using Statistics

Module NewtonApproach  

Λ_func   = (x, a, b) -> 1.0 + b / 2.0 * (tanh(x/a) - 1.0)
∂Λ_func  = (x, a, b) -> b / a * (sech(x/a)^2.0)  / 2.0
∂∂Λ_func = (x, a, b) -> - (b / a^2.0)  * (sech(x/a)^2.0) * tanh(x/a)

Δ = (i, j, p) -> (i == (mod(j-1, p) + 1)) ? 1 : 0

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)
@inline mod12(n) = mod(n-1, 12) + 1

mutable struct bundle
    δ      :: Array{T, 2}
    δ_p1   :: Array{T, 2}
    γ      :: Array{T, 2}
    ∇ϵ     :: Array{T, 2}
    N      :: Integer
    period :: Integer

    tool_bundle = function(;N::Integer, period::Integer) 

        δ    = zeros(dtype, period, N)
        δ_p1 = copy(δ)

        for k=1:period, i=1:N
            δ[k, i]    = Δ(k, i  , period)
            δ_p1[k, i] = Δ(k, i+1, period)
        end
        γ = (δ_p1 - δ) / dt

        ∇ϵ = zeros(dtype, N, period*2)
        for i=1:N, j=1:period
            ∇ϵ[i, j+period] = - δ[j, i]
        end

        #∇∇ϵ = zeros(dtype, N, 24, 24)
        #for i=1:N, j=1:24, k=1:24
        #    if 
        #    ∇∇ϵ[i, j, k] = - δ[j, i]
        #end


        new(δ, δ_p1, γ, ∇ϵ, N, period)
    end
end

function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end

function g_and_∇g(bundle; h, Q_ph, θ, θ_p1, S_ph, B_ph, a, b, Δt)


    repeat_fill!(h, h[1:12])
    repeat_fill!(Q_ph, Q_ph[1:12])

    h_p1 = circshift(h, -1)

    ∂h∂t = (h_p1 - h) / Δt

    Λ   =   Λ_func.(∂h∂t, a, b)
    ∂Λ  =  ∂Λ_func.(∂h∂t, a, b)
    ∂∂Λ = ∂∂Λ_func.(∂h∂t, a, b)
    
    ϵ =  (
        h    / 2Δt .* ( θ_p1 .* (1.0 .- Λ) - θ .* (1.0 .+ Λ) )
        + h_p1 / 2Δt .* ( θ_p1 .* (1.0 .+ Λ) - θ .* (1.0 .- Λ) )
        - S_ph - B_ph - Q_ph
    )

    for i=1:N, j=1:bundle.period
        bundle.∇ϵ[i, j] =  (
            δ[j, i]    / 2Δt * ( θ_p1[i] * (1.0 - Λ[i]) - θ[i] * (1.0 + Λ[i]) )
            + δ_p1[j, i] / 2Δt * ( θ_p1[i] * (1.0 + Λ[i]) - θ[i] * (1.0 - Λ[i]) )
            + ∂h∂t[i] * (θ[i] + θ_p1[i]) / 2.0 * γ[j, i] * ∂Λ[i]
        )
    end

    g  = bundle.∇ϵ' * ϵ
    ∇g = bundle.∇ϵ' * bundle.∇ϵ

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

function fit(
    pts_per_year :: Integer,
    F            :: Array{T},
    θ            :: Array{T},
    Δt           :: T,
    S,
    B,
    ϕ = 0
) where T <: AbstractFloat
 
    if mod(length(F), pts_per_year) != 0
        throw(ArgumentError("Data length should be multiple of [pts_per_year]"))
    end

    years = Int(length(F) / pts_per_year)

    N       = (years-2) * pts_per_year      # Discard the first and last year
    beg_t   = pts_per_year + 1              # Jan of second year

    rng1 = collect(beg_t:beg_t+N-1)
    rng2 = rng1+1
 
    _S_ph = (S[rng1] + S[rng2]) / 2.0
    _B_ph = (B[rng1] + B[rng2]) / 2.0

    _θd   = _S_ph * 0.0
    _θd[1:12] .= 273.15
    repeat_fill!(_θd, _θd[1:12])

    _θ    = θ[i, j, rng1] - _θd
    _θ_p1 = θ[i, j, rng2] - _θd

    global x_mem
    x_mem[ 1:12] = output_h[i, j, :] 
    x_mem[13:24] = output_Q[i, j, :] 


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
                θ = _θ,
                θ_p1 = _θ_p1,
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

        println((x_mem[2:12] - x_mem[1:11]) / dt)
 
        println("Q std: ", std(x_mem[13:24]))
    end

    output_h[i, j, :] = x_mem[ 1:12]
    output_Q[i, j, :] = x_mem[13:24]

end

