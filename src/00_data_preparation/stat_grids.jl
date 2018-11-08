include("../01_config/general_config.jl")

tos = readModelVar("tos", (:, :, 1))

grid_points = length(tos)
oc_grid_points = sum(isfinite.(tos))

println("tos has ", grid_points, " points.")
println("tos has ", oc_grid_points, " ocean points.")





