include("config.jl")
include("NetCDFHelper.jl")

using NetCDF
using LinearAlgebra

N       = Int((nyrs-2) * 12) # Discard the first and last year
beg_t   = Int(13)            # Jan of second year



output_converge = 0

dt2 = 2.0 * dt

# Gradient Descent method parameters
iter_max = 1000
#update_bounds = [0.16, 0.16]
ϵ_converge_ratio_threshold = 1e-11
converge_count_threshold = 10
ϵ_mem = Inf
η = 1e-8

println(beg_t+1)
println(beg_t+N - 1)

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

converge_count = 0
# For each iteration
for k = 1 : iter_max
    global h, ∂θ∂t, ϵ_mem, converge_count
    # Calculate ϵ and ϵ^2
    ϵ =  h * ∂θ∂t - S - B
    ϵ2_sum = sum(ϵ.^2.0)
    ϵ_now = (ϵ2_sum / N)^0.5

    #@printf("Determine convergence condition.\n")
    # Determine if iteration should stop or not
    Δϵ_ratio = (ϵ_now - ϵ_mem) / ϵ_mem

    @printf("Iter: %d. h: %.2f, ϵ_now: %f; Δϵ_ratio * 100: %f\n", k, h, ϵ_now, Δϵ_ratio * 100.0)

    #@printf("Converging?\n")
    if Δϵ_ratio < 0.0 && abs(Δϵ_ratio) < ϵ_converge_ratio_threshold  # if it is converging
        converge_count +=1
    else # if it is diverging
        converge_count = 0
    end

    #@printf("Converge count threshold reached?\n")
    if converge_count >= converge_count_threshold
        break
    end
   
    #@printf("Gradient descent.\n")
    # Calculate ∂Post∂h, ∂Post∂Q
    # Assume flat prior for now
    
    ∂ϵ∂h = ∂θ∂t
    ∂Post∂h = - sum(ϵ .* ∂ϵ∂h)

    #@printf("Update h and Q.\n")
    # Update h
    h += ∂Post∂h * η

    ϵ_mem = ϵ_now 
end

@printf("Result h: %.2f m.\n", h)

   
