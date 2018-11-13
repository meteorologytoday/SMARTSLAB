#using LinearAlgebra
#using Statistics

include("../../lib/Newton.jl")


module NewtonApproach 
using ForwardDiff
using ..NewtonMethod

eucLen    = x -> (sum(x.^2.0))^(0.5)
normalize = x -> x / eucLen(x)

we_func = (x, a) -> (x >= 0) ? x : a*x


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
        println("Enter,", typeof(h))
        println("Enter,", typeof(Q_ph))
        println("Enter,", x)
        repeat_fill!(h,    x[ 1:12])
        repeat_fill!(Q_ph, x[13:24])

        println("Enter1")
        h_p1 = circshift(h, -1)
        h_ph = (h_p1 + h) / 2.0

        ∂h∂t = (h_p1 - h) / Δt

        println("Enter2")
        we =   we_func.(∂h∂t, a)
        println("Enter3")
        
        ϵ =  (
            h_ph .* ∂θ∂t_ph
            + (θ_ph .- θd) .* we
            - S_ph - B_ph - Q_ph
        )

        return ϵ' * ϵ
    end

    f_and_∇f = function(x)
        println("Now we are here")
        f  = ForwardDiff.gradient(calϵ2, x)
        println("Now we are here")
        ∇f = ForwardDiff.hessian( calϵ2, x)

        return f, ∇f
    end

    println(x_mem)
    println(f_and_∇f(x_mem))
    println("Gonna call newton.fit")
    x_mem[:] = NewtonMethod.fit(;
        f_and_∇f = f_and_∇f,
        η        = η,
        x0       = x_mem,
        max      = max,
        verbose  = verbose
    )
    println("done")


    return x_mem[ 1:period], x_mem[period+1:end]

end

end
