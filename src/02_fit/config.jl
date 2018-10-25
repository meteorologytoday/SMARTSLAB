rng = 1:120
rng = Colon()

F = readModelVar("hfds", (:, :, rng))
θ = readModelVar("tos", (:, :, rng)) * ρ * c_p

println()

