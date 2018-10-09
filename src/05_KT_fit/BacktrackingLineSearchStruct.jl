
# This struct implements backtrack line search.
# The algorithm is referred to "Convex Optimization p.464" by Boyd and Vandengerghe.

mutable struct BacktrackingLineSearchStruct
    N     :: Integer
    α     :: AbstractFloat           # the tolerance factor must in open interval (0, 0.5)
    β     :: AbstractFloat           # the shrinking factor must in open interval (0, 1.0)
    t     :: AbstractFloat           # step size
    η     :: AbstractFloat           # stopping criterion
    BacktrackingLineSearchStruct(; N, α, β, t, η) = new(N, α, β, t, η)
end

