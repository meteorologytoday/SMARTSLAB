

function LR_shallow_water(;
    years        :: Integer,
    pts_per_year :: Integer,
    F            :: Array{T},
    θ            :: Array{T},
    Δt           :: T,
) where T <: AbstractFloat


N       = (years-2) * pts_per_year      # Discard the first and last year
beg_t   = pts_per_year + 1              # Jan of second year

params  = pts_per_year * 2

ϕ = zeros(dtype, N, params)             # for h and Q

mod_cycle(n) = mod(n-1, pts_per_year) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

ϕ .= 0.0

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

return β[1:pts_per_year], β[pts_per_year+1:end]

end







