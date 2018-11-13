using PyPlot




Λ = (x, a, b) -> 1.0 + b / 2.0 * (tanh(x/a) - 1.0)


ab = [
    1.0 0.0;
    1.0 0.5;
    1.0 1.0;
    0.5 1.0;
    1e-2 1.0;
]

x = range(-1, stop=1, length=1000) |> collect
for i = 1:size(ab)[1]
    plt[:plot](x, Λ.(x, ab[i,1], ab[i,2]))
end
plt[:legend]()
plt[:show]()
