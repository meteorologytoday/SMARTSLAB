

function LinearRegression(
    x :: Array{T, 1},
    y :: Array{T, 1}
) where T <: AbstractFloat 

N = length(x)

ϕ = zeros(N, 2)
ϕ[:, 1] .= 1.0
ϕ[:, 2] = x

# ϕ β = y => β = ϕ \ y

return ϕ \ y

end
