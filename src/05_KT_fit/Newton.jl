#using LinearAlgebra

function Newton(;
    #g  :: Function,             # Target function g which receives a vector
    #∇g :: Function,             # Jacobian of g which receives a vector
    g_and_∇g  :: Function,       # Target function g and its Jacobian
    η  :: T,        # Threshold
    x0 :: Array{T, 1}, # Initial guess
    max:: Integer               # Maximum iteration
) where T <: AbstractFloat

    local x = x0 * 1.0
    local if_converge = false

    for i = 1:max
        g, ∇g = g_and_∇g(x)
        Δx = - ∇g \ g
        println(x)
        if (Δx' * Δx)^0.5 >= η
            x += Δx
        else
            if_converge = true
            break
        end
    end

    if if_converge == false
       throw(Excpetion("Not converge")) 
    end

    return x
end
