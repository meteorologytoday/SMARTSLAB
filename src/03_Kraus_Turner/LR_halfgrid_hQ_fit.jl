include("simulate_given_dhdt.jl")


N       = Int((years-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year

ϕ = zeros(dtype, N, 12 * 2)   # for h and Q
β = zeros(dtype, 24)
@inline mod12(n) = mod(n-1, 12) + 1

rng1 = collect(beg_t:(beg_t + N -1))
rng2 = rng1 .+ 1

F = S + B
θs .-= θd

#println(θs / ρ / c_p)
#println(θs / Δt)
let
    ϕ .= 0.0

    for t = 1:N

        # ϕ_h
        ϕ[t, mod12(t  )] = - θs[(beg_t + t - 1)    ] / Δt
        ϕ[t, mod12(t+1)] =   θs[(beg_t + t - 1) + 1] / Δt

        # ϕ_Q
        ϕ[t, 12 + mod12(t  )] =  - 1.0

    end


#    prtArr(ϕ)
    

    _F =  (F[rng1] + F[rng2]) / 2.0
    _F =  F[rng1]
    println(_F)
    # Solve normal equation
    # ϕ β = F => β = ϕ \ F
    β[:] = ϕ \ _F

end

println("## h")
println(β[1:12])

println("## Q")
println(β[13:24])


