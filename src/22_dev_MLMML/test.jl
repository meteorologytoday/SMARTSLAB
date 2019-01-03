include("MLMML.jl")

using .MLMML

D  = 1000.0
N  = 1001
zs = collect(Float64, range(0.0, stop=D, length=N))


obj = MLMML.OceanColumn(zs);
