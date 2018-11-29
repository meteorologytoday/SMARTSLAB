function LR_shallow_water!(;
    period :: Integer,
    F            :: Array{T},
    θ            :: Array{T},
    Δt           :: T,
    ϕ = 0
) where T <: AbstractFloat


if mod(length(F), period) != 0
    throw(ArgumentError("Data length should be multiple of [period]"))
end

years = Int(length(F) / period)

N       = (years-2) * period      # Discard the first and last year
beg_t   = period + 1              # Jan of second year

params  = period * 2

if ϕ == 0
    ϕ = zeros(dtype, N, params)             # for h and Q
elseif  (! isa(ϕ, Array{T,2})) || ( size(ϕ) != (N, params) )
    throw(ArgumentError("ϕ is not the correct type or shape"))
end

mod_cycle(n) = mod(n-1, period) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

ϕ .= 0.0

#println(length(θ))
for t = 1:N

    # ϕ_h
    ϕ[t, mod_cycle(t  )] = - θ[(beg_t + t - 1)    ] / Δt
    ϕ[t, mod_cycle(t+1)] =   θ[(beg_t + t - 1) + 1] / Δt

    # ϕ_Q
    ϕ[t, 12 + mod_cycle(t)] =  - 1.0

end


#prtArr(ϕ)


# ϕ β = F => β = ϕ \ F

_F =  (F[rng1] + F[rng2]) / 2.0
β = ϕ \ _F

#_F =  F[rng1]
#println(_F)
# Solve normal equation




#println(ϕ * β - _F)


#∂θs∂t = ((circshift(θs, -1) - θs)[beg_t:beg_t+N-1]) / Δt
#println(sum(_F .* ∂θs∂t) / sum(∂θs∂t.^2.0))
#println(_F ./ ∂θs∂t)

return β[1:period], β[period+1:end]

end







