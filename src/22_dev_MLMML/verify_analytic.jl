include("MLMML.jl")
using .MLMML
using PyPlot


PERIOD = 86400.0 * 360.0
J0 = - 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p
ω  = 2π / PERIOD
h0 = 30.0
m  = 30.0 / 5000.0
k  = MLMML.getTKE(fric_u=MLMML.getFricU(ua=0.0))

sol = t -> (m * h0)^(-1.0) * (√(k^2.0 + (2.0*m*h0^2.0*(-J0))/ω * (1.0 - cos(ω*t))) - k)
ts = collect(Float64, range(0.0, stop=180.0*86400.0, length=1000))
analytical_sol = sol.(ts)

using DifferentialEquations
f(h, p, t) = - h * J0 * sin(ω*t) / (m * (h^2.0 - h0^2.0) / 2.0 + k)
prob = ODEProblem(f, h0, (ts[1], ts[end]))
numerical_sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)



plt[:figure]()
plt[:plot](ts/86400.0, - analytical_sol, label="analytical")
plt[:plot](numerical_sol.t/86400.0, - (numerical_sol.u .- h0), label="numerical")
plt[:legend]()
plt[:show]()
