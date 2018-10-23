function extend(a::AbstractArray, len::Int)
    if length(a) < len
        a = repeat(a, outer=ceil(Int, len / length(a)))
    end

    return a[1:len]
end



function f(;
    h       :: TYPE,
    θs      :: TYPE,
    θd      :: TYPE,
    β       :: TYPE,
    S       :: TYPE,
    B       :: TYPE,
    ∂h∂t    :: TYPE,
) where TYPE <: AbstractFloat

    #Λ = (∂h∂t > 0.0) ? 1.0 : 0
    Λ = 0.0
    return (S * (1.0 - exp(-β * h)) + B - (θs - θd) * ∂h∂t * Λ) / h

end


function simulation(;
    θs_init :: TYPE,
    θd      :: TYPE,
    Δt      :: TYPE,
    β       :: TYPE,
    S       :: Array{TYPE, 1},
    B       :: Array{TYPE, 1},
    h       :: Array{TYPE, 1},
    ret_len :: Integer
) where TYPE <: AbstractFloat

    ∂h∂t = (circshift(h, -1) - circshift(h, 1)) / (2.0 * Δt)
    ∂h∂t[1] = (h[2] - h[1]) / Δt
    ∂h∂t[end] = (h[end] - h[end - 1]) / Δt

    θs = zeros(TYPE, ret_len)
    θs[1] = θs_init

    for i = 1 : ret_len-1
        
        ∂θs∂t = f(
            h=h[i],
            θs=θs[i],
            θd=θd,
            β=β,
            S=S[i],
            B=B[i],
            ∂h∂t=∂h∂t[i],
        )

        θs[i+1] = θs[i] + ∂θs∂t * Δt

    end

    return θs
end




