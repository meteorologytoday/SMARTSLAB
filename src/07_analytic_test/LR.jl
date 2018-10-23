#include("analytic_given_dhdt.jl")
include("simulate_given_dhdt.jl")


N       = Int((years-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

ϕ = zeros(dtype, N, 1)
β = zeros(dtype, 1)
@inline mod12(n) = mod(n-1, 12) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

F = S + B
θs .-= θd

#println(θs / ρ / c_p)
#println(θs / Δt)
let
    ϕ .= 0.0
   
    println("∂θs∂t:") 
    ∂θs∂t = ((circshift(θs, -1) - θs)[beg_t:beg_t+N-1]) / Δt

    for t = 1:N

        # ϕ_h
        ϕ[t, 1] = ∂θs∂t[t]

    end


    prtArr(ϕ)
    

    _F =  (F[rng1] + F[rng2]) / 2.0
    _F =  F[rng1]
    #println(_F)
    # Solve normal equation
    # ϕ β = F => β = ϕ \ F
    β[:] = ϕ \ _F

    #println(ϕ * β - _F)

    

    #println(sum(_F .* ∂θs∂t) / sum(∂θs∂t.^2.0))
    #println(_F ./ ∂θs∂t)
    
end

println("## h")
println(β[1])

#println("## Q")
#println(β[13:24])


