function SLAB!(;
    period        :: Integer,
    F             :: Array{T},
    θ             :: Array{T},
    Δt            :: T,
    reinterpolate :: Bool = false,
    ϕ = 0
) where T <: AbstractFloat


if mod(length(F), period) != 0
    throw(ArgumentError("Data length should be multiple of [period]."))
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

for t = 1:N

    # ϕ_hΔ
    ϕ[t, mod_cycle(t)] = (θ[(beg_t + t - 1) + 1] - θ[(beg_t + t - 1)]) / Δt

    # ϕ_QΔ
    ϕ[t, 12 + mod_cycle(t)] =  - 1.0

end



#prtArr(ϕ)


# ϕ β = F => β = ϕ \ F

_F =  (F[rng1] + F[rng2]) / 2.0
β = ϕ \ _F

if reinterpolate
    if mod(period, 2) == 0
        h = β[1:period]
        Q = β[period+1:2*period]
        β[       1:  period] = ( circshift(h, -1) + h ) / 2.0
        β[period+1:2*period] = ( circshift(Q, -1) + Q ) / 2.0
    else
        throw(ErrorException("Have not implement interpolation for odd period yet."))
    end
end



return β

end







