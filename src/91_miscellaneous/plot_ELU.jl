using PyPlot


ELU(x, α) = 

Λ = (x, a, b) -> 1.0 + b / 2.0 * (tanh(x/a) - 1.0)


ab = [
    1.0 0.0;
    1.0 0.5;
    1.0 1.0;
    0.5 1.0;
    1e-2 1.0;
]

fig, ax = plt[:subplots](1, 2, sharex=true, figsize=(16,10))


x = range(-1, stop=1, length=1000) |> collect
for i = 1:size(ab)[1]
    y = Λ.(x, ab[i,1], ab[i,2])
    ax[1][:plot](x, y)
    ax[2][:plot](x, y .* x)
end
ax[1][:legend]()
ax[2][:legend]()
plt[:show]()
