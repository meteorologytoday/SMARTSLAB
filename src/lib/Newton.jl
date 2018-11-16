module NewtonMethod

export NotConvergeException, Newton

struct NotConvergeException <: Exception
end

function fit(;
    f_and_∇f  :: Function,    # Target function g and its Jacobian
    η         :: T,           # Threshold
    x0        :: Array{T, 1}, # Initial guess
    max       :: Integer,     # Maximum iteration
    verbose   :: Bool=false
) where T <: AbstractFloat
    #println("Newton method!")
    #println(typeof(f_and_∇f))
    local x = x0 * 1.0
    local if_converge = false

    for i = 1:max
        if verbose
            println("Newton method iteration: ", i)
            println("x: ", x)
        end
        f, ∇f = f_and_∇f(x)
        Δx = - ∇f \ f

        eulen_Δx = (Δx' * Δx)^0.5
        verbose && println("|Δx| = ", eulen_Δx)
        if eulen_Δx >= η
            x += Δx
        else
            if_converge = true
            break
        end
    end

    if if_converge == false
       throw(NotConvergeException()) 
    end

    return x
end


end
