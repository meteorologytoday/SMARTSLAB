#using LinearAlgebra
#using Statistics

include("../../lib/Newton.jl")


module NewtonApproach 
using ForwardDiff
using ..NewtonMethod

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)

function repeat_fill!(to::AbstractArray, fr::AbstractArray)

    len_fr = length(fr)
    println("Repeat:", len_fr, ", ", length(to))
    println(typeof(fr), ", ", typeof(to))
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
    println("Repeat done")
end


function fit(;
    N        :: Integer,
    period   :: Integer,
    beg_t    :: Integer,
    Δt       :: T,
    init_h   :: Array{T},
    init_Q   :: Array{T},
    init_θd  :: T,
    θ        :: Array{T},
    S        :: Array{T},
    B        :: Array{T},
    a        :: T,
    max      :: Integer,
    η        :: T,
    σ_ϵ      :: T,
    σ_Q      :: T,
    σ_h      :: T,
    σ_θd     :: T,
    h_rng    :: Array{T},
    θd_rng   :: Array{T},
    verbose  :: Bool = false
) where T <: AbstractFloat


    if mod(length(θ), period) != 0 || mod(N, period) !=0
        throw(ArgumentError("Data length should be multiple of [period]"))
    end

    reduced_years = Int(N / period)

    rng1 = collect(beg_t:beg_t+N-1)
    rng2 = rng1 .+ 1

    σ²_ϵ  = σ_ϵ  ^ 2
    σ²_Q  = σ_Q  ^ 2
    σ²_h  = σ_h  ^ 2
    σ²_θd = σ_θd ^ 2

    # Extract fixed data
    _S_ph = (S[rng1] + S[rng2]) / 2.0
    _B_ph = (B[rng1] + B[rng2]) / 2.0

    _θ    = θ[rng1]
    _θ_p1 = θ[rng2]

    _∂θ∂t_ph = (_θ_p1 - _θ) / Δt
    _θ_ph    = (_θ_p1 + _θ) / 2.0

    x_mem = zeros(T, 2*period + 1)

    x_mem[ 1        :   period] = init_h 
    x_mem[period+1  : 2*period] = init_Q
    x_mem[2*period+1]           = init_θd

    calLogPost = function(x)

        local L = 0.0

        h    = repeat(x[1:period],          outer=(reduced_years,))
        Q_ph = repeat(x[period+1:2*period], outer=(reduced_years,))
        θd   = x[end]

        h_p1 = circshift(h, -1)

        h_p1 = circshift(h, -1)
        h_ph = (h_p1 + h) / 2.0

        ∂h∂t = (h_p1 - h) / Δt

        we = ((∂h∂t .≥ 0) + (∂h∂t .< 0) .* a) .* ∂h∂t
        
        ϵ =  (
            h_ph .* _∂θ∂t_ph
            + (_θ_ph .- θd) .* we
            -  _S_ph  - _B_ph - Q_ph
        )

        L += - ϵ' * ϵ / σ²_ϵ

        # Add prior of h
        _h = h[1:period]
        L += - sum( ((_h .< h_rng[1]) .* (_h .- h_rng[1])).^2 ) / σ²_h
        L += - sum( ((_h .> h_rng[2]) .* (_h .- h_rng[2])).^2 ) / σ²_h

        # Add prior of Q
        L += - sum( (Q_ph[1:period]).^2 ) / σ²_Q

        # Add prior of θd
        L += (θd < θd_rng[1]) ? - (θd - θd_rng[1])^2.0 / σ²_θd : 0.0
        L += (θd > θd_rng[2]) ? - (θd - θd_rng[2])^2.0 / σ²_θd : 0.0

        #=
        # Combined version (before applied product rule)
        ϵ = (
            (h_p1 .* _θ_p1 - h .* _θ) / Δt - _S_ph - _B_ph - Q_ph
        )
        =#

        return L
    end

    f_and_∇f = function(x)
        f  = ForwardDiff.gradient(calLogPost, x)
        ∇f = ForwardDiff.hessian( calLogPost, x)

        return f, ∇f
    end

    x_mem[:] = NewtonMethod.fit(;
        f_and_∇f = f_and_∇f,
        η        = η,
        x0       = x_mem,
        max      = max,
        verbose  = verbose
    )


    return x_mem

end

end
