include("MLMML.jl")

using .MLMML

D  = 1000.0
N  = 11
zs = collect(Float64, range(0.0, stop=-D, length=N))


oc = MLMML.OceanColumn(zs);

# setup b anomaly

#MLMML.stepOceanColumn


using PyPlot

fig, ax = plt[:subplots](1, 1, figsize=(8,6))

for i=1:length(oc.bs)
    ax[:plot]([oc.bs, oc.bs], [oc.zs[i], oc.zs[i+1]], color="k",)
end

plt[:show]()
