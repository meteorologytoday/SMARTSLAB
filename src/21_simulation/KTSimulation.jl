
module KTSimulation

function f(;
    θ    :: G,
    θd   :: G,
    F    :: G,
    Q    :: G,  
    h    :: G,
    ∂h∂t :: G
) where G <: AbstractFloat

    return ((F + Q) - (θ - θd) * (∂h∂t > 0.0) * ∂h∂t) / h

end

function repeat_fill!(
    to   :: Array{T, 1},
    from :: Array{T, 1}
) where T <: Any

    period = length(from)
    for i = 1:length(to)
        to[i] = from[ mod(i-1, period) + 1]
    end
end

function repeat_fill(
    arr     :: Array{T, 1},
    ret_len :: Int
) where T <: Any

    arr2 = zeros(T, ret_len)
    repeat_fill!(arr2, arr)

    return arr2
end


function linear_interpolate(a::Array{G, 1}, i::G) where G <: AbstractFloat
    return (i - floor(i)) * a[ceil(Int, i)] + (ceil(i) - i) * a[floor(Int, i)]
end


function run(;
    θ_init  :: G,
    θd      :: G,
    Δt      :: G,
    F       :: Array{G, 1},
    Q       :: Array{G, 1},
    h       :: Array{G, 1},
    period  :: Int,
    ret_len :: Int
) where G <: AbstractFloat

    ∂h∂t = (circshift(h, -1) - circshift(h, 1)) / (2.0 * Δt)
 
    θ = zeros(G, ret_len)
    Q = repeat_fill(Q, ret_len)
    F = repeat_fill(F, ret_len)
    h = repeat_fill(h, ret_len)
    ∂h∂t = repeat_fill(∂h∂t, ret_len)

    θ[1] = θ_init 
    for i = 1 : ret_len-1

        # Euler method
        Δθ = Δt * f(;
            θ    = θ[i],
            θd   = θd,
            F    = F[i],
            Q    = Q[i],  
            h    = h[i],
            ∂h∂t = ∂h∂t[i]
        )

        θ[i+1] = θ[i] + Δθ

    end 

    return Dict(
        "θ" => θ,
        "F" => F,
        "Q" => Q,
        "h" => h,
        "∂h∂t" => ∂h∂t,
    )
end

end
