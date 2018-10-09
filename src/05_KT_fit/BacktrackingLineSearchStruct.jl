
# This ηCtl implements backtrack line search.
# The algorithm is referred to "Convex Optimization p.464" by Boyd and Vandengerghe.


mutable struct BacktrackingLineSearchStruct

    N     :: Integer
    α     :: AbstractFloat           # the tolerance factor must in open interval (0, 0.5)
    β     :: AbstractFloat           # the shrinking factor must in open interval (0, 1.0)
    t     :: AbstractFloat           # step size
    η     :: AbstractFloat           # stopping criterion
    R     :: Array{AbstractFloat,1}  # array of residue
    ptr   :: Integer                 # pointer to the latest R


    function BacktrackingLineSearchStruct(N, α, β, t, η)
        return new(N, α, β, t, η, zeros(Float32, N), 1)
    end
end


function rewind(ctl :: BacktrackingLineSearchStruct)
    if N > 1
        ctl.R[ptr] = NaN
        ctl.N -= 1
    end
end

function shrink(ctl :: BacktrackingLineSearchStruct)
    ctl.t *= ctl.β
end
