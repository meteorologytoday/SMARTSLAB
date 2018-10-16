include("Newton.jl")

using LinearAlgebra

g = x -> [2.0 * x[1] + 1, x[2]]
Jg = x -> (I + zeros(2,2))  * 2.0


ans = Newton(
    g  = g,
    Jg = Jg,
    η  = 1e-45,
    x0 = [2.0, 100.0],
    max = 1000
)

println(ans)
