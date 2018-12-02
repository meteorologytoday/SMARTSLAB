include("Newton.jl")

using .NewtonMethod
using LinearAlgebra
using ForwardDiff

g = x -> sum(x.^2.0)

f_and_∇f = x -> (ForwardDiff.gradient(g, x), ForwardDiff.hessian(g, x))



ans = NewtonMethod.fit(;
    f_and_∇f = f_and_∇f,
    η  = 1e-45,
    x0 = [2.0, 100.0],
    max = 1000,
    verbose = true
)

println(ans)
