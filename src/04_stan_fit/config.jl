rng = 1:1200
#rng = Colon()
years = 100
F = readModelVar("hfds", (:, :, rng))
θ = readModelVar("tos", (:, :, rng)) * ρ * c_p

println()

unconverge_value = -999999.0

