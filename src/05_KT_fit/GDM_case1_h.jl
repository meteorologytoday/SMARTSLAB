include("config.jl")
include("NetCDFHelper.jl")
include("BacktrackingLineSearchStruct.jl")
using NetCDF
using LinearAlgebra

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year


dt2 = 2.0 * dt

# Gradient Descent method parameters
iter_max = 100

BLSS = BacktrackingLineSearchStruct(iter_max, 0.1, 0.8, 10.0, 1e-3)


# Assign h and Q initial condition
h = 100.0
∂θ∂t = (T_star[:, :, beg_t+1:beg_t+N] - T_star[:, :, beg_t-1:beg_t+N-2]) / dt2
S = S[:, :, beg_t:beg_t+N-1]
B = B[:, :, beg_t:beg_t+N-1]

idx = isfinite.(T_star[:, :, beg_t:beg_t+N-1])

∂θ∂t = ∂θ∂t[idx]
S = S[idx]
B = B[idx]

println(size(∂θ∂t))
println(size(S))

println(sum(isnan.(S)))


function getϵandϵ2(h, ∂θ∂t, S, B)
    ϵ  = h * ∂θ∂t - S - B
    ϵ2 = sum(ϵ.^2.0)/length(∂θ∂t)
    return ϵ, ϵ2
end

# For each iteration
for k = 1 : iter_max
    global h, ∂θ∂t

    if_update = false

    ϵ, ϵ2 = getϵandϵ2(h, ∂θ∂t, S, B)

    ∂ϵ∂h       = ∂θ∂t
    ∂ϵ2∂h      = sum(ϵ .* ∂ϵ∂h) / length(ϵ)
    ∂ϵ2∂h_unit = ∂ϵ2∂h / sum(∂ϵ2∂h.^2.0)^0.5
   
    abs_∇ = abs(∂ϵ2∂h)


 
    # 1: Test if this iteration should stop
    if abs_∇ < BLSS.η
        @printf("The stop condition is met: abs(∂ϵ2∂h) = %f\n", abs_∇)
        break
    end

    # 2: Compare values of ϵ2 of changing h and Taylor extrapolation
    new_h  = h - BLSS.t * ∂ϵ2∂h_unit

    _, ϵ2_from_new_h  = getϵandϵ2(new_h, ∂θ∂t, S, B)
    ϵ2_from_taylor = ϵ2 - BLSS.α * BLSS.t * ∂ϵ2∂h_unit

    if ϵ2_from_new_h > ϵ2_from_taylor
        BLSS.t *= BLSS.β
    else
        h = new_h
        if_update = true
    end


    @printf("Iter: %d. Update? %s. h: %.2f, ϵ2_from_new_h: %f, ϵ2_from_taylor: %f, BLSS.t: %f, abs_∇: %f\n",
        k,
        ((if_update == true) ? "YES" : "NO"),
        h,
        ϵ2_from_new_h,
        ϵ2_from_taylor,
        BLSS.t,
        abs_∇
    )
end

@printf("Result h: %.2f m.\n", h)

   

