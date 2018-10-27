

module LR_breakdown

@inline modc(n, c) = mod(n-1, c) + 1
function int2bin_arr(
    x      :: Integer,
    digits :: Integer,
)


arr = zeros(Integer, digits)


for i in 1:digits
    arr[end-i+1] = x >> (i-1) & 0b1
end

return arr
end



function gen_bin_arr(
    digits :: Integer
)

    max = 1 << digits
    
    arr = zeros(Integer, max, digits)

    for i = 1 : max
        arr[i, :] = int2bin_arr(i-1, digits)
    end

    return arr
end


function set_ϕ!(;
    ϕ  :: Array{T, 2},
    Λ  :: Array{Integer, 1},
    θ  :: Array{T, 1},
    Δt :: T
) where T <: AbstractFloat

N, M = size(ϕ)

for i = 1 : N

    ic  = modc(i, M)
    iic = modc(i+1, M)
    
    λ  = Λ[ic]

    ϕ[i,  ic] = (θ[i+1] * (1.0 - λ) - θ[i] * (1.0 + λ)) / (2.0 * Δt)
    ϕ[i, iic] = (θ[i+1] * (1.0 + λ) - θ[i] * (1.0 - λ)) / (2.0 * Δt)

end

end

function solve(;
    years        :: Integer,
    pts_per_year :: Integer,
    F            :: Array{T},
    θ            :: Array{T},
    Δt           :: T,
) where T <: AbstractFloat

params  = pts_per_year * 2

Λ         = gen_bin_arr(pts_per_year)
ϵ2_record = zeros(T, size(Λ)[1])
matches   = zeros(Integer, size(Λ)[1])
β         = zeros(T, size(Λ)[1], params)

ϵ2_record .= Inf


N       = (years-2) * pts_per_year      # Discard the first and last year
beg_t   = pts_per_year + 1              # Jan of second year



ϕ        = zeros(T, N, params)                   # h and Q
ϕ_subset = zeros(T, N, pts_per_year)             # h

mod_cycle(n) = mod(n-1, pts_per_year) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

for t = 1:N
    # ϕ_Q
    ϕ[t, 12 + modc(t, pts_per_year)] =  - 1.0
end

_F =  (F[rng1] + F[rng2]) / 2.0

for k = 1:size(Λ)[1]
    set_ϕ!(
        ϕ  = ϕ_subset,
        Λ  = Λ[k, :],
        θ  = θ,
        Δt = Δt
    )

    ϕ[:, 1:pts_per_year] = ϕ_subset

    # ϕ β = F => β = ϕ \ F
    β[k, :] = ϕ \ _F

    # Test if this solution matches Λ constraint
    h = β[k, 1:pts_per_year]
    Λ_LR = convert(Array{Integer}, (circshift(h, -1) - h) .> 0)


    matches[k] = sum(Λ_LR .== Λ[k, :]) == pts_per_year
    if !matches[k]
        #continue
    end
    # Record residue
    ϵ = ϕ * β[k, :] - _F
    ϵ2_record[k] = ϵ' * ϵ

    println("h:", h, "; Λ:", Λ[k,:], "; ϵ2:", ϵ2_record[k],)    

end




end


end
