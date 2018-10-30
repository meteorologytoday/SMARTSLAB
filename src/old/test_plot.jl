using PyPlot
x = collect(1:10)
y = x.^2
plt[:plot](x, y)
plt[:show]()
