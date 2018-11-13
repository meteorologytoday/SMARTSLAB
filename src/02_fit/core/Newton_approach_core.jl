using LinearAlgebra
using Statistics

include("../../lib/Newton.jl")


module NewtonApproach 

using ..NewtonMethod

Λ_func   = (x, a, b) -> 1.0 + b / 2.0 * (tanh(x/a) - 1.0)

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)

we_func = (x, a) -> (x >= 0) ? x : a*x


function repeat_fill!(to::AbstractArray, fr::AbstractArray)
    len_fr = length(fr)
    for i = 1 : length(to)
        to[i] = fr[mod(i-1, len_fr)+1]
    end 
end


function fit(;
    N        :: Integer,
    period   :: Integer,
    beg_t    :: Integer,
    Δt       :: T,
    init_h   :: Array{T},
    init_Q   :: Array{T},
    θ        :: Array{T},
    S        :: Array{T},
    B        :: Array{T},
    θd       :: T,
    a        :: T,
    max      :: Integer,
    η        :: T,
    verbose  :: Bool = false
) where T <: AbstractFloat

    if mod(length(θ), period) != 0
        throw(ArgumentError("Data length should be multiple of [pts_per_year]"))
    end

    years = Int(length(S) / period)

    rng1 = collect(beg_t:beg_t+N-1)
    rng2 = rng1 .+ 1


    # Extract fixed data

    _S_ph = (S[rng1] + S[rng2]) / 2.0
    _B_ph = (B[rng1] + B[rng2]) / 2.0

    _θd   = zeros(T, N)
    _θd[1:period] .= θd
    repeat_fill!(_θd, _θd[1:period])

    _θ    = θ[rng1]
    _θ_p1 = θ[rng2]

    _∂θ∂t_ph = (_θ_p1 - _θ) / Δt

    x_mem = zeros(T, 2*period) 
    x_mem[ 1       : period] = init_h 
    x_mem[period+1 : end   ] = init_Q

    h    = zeros(T, N)
    Q_ph = zeros(T, N)

    calϵ2 = function(x)
        repeat_fill!(h,    x[ 1:12])
        repeat_fill!(Q_ph, x[13:24])

        h_p1 = circshift(h, -1)
        h_ph = (h_p1 + h) / 2.0

        ∂h∂t = (h_p1 - h) / Δt

        we =   we_func.(∂h∂t, a)
        
        ϵ =  (
            h_ph .* ∂θ∂t_ph
            + (θ_ph .- θd) .* we
            - S_ph - B_ph - Q_ph
        )

        return ϵ' * ϵ
    end

    f_and_∇f = function(x)
        f  = ForwardDiff.gradient(calϵ2, x)
        ∇f = ForwardDiff.hessian( calϵ2, x)

        return f, ∇f
    end


    x_mem[:] = Newton.fit(;
        f_and_∇f = f_and_∇f,
        η        = η,
        x0       = x_mem,
        max      = max,
        verbose  = verbose
    )

    return x_mem[ 1:period], x_mem[period+1:end]

end

end
