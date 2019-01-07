include("MLMML.jl")

using Printf
using .MLMML

D  = 2000.0
N  = 21
zs = collect(Float64, range(0.0, stop=-D, length=N))

t = collect(Float64, range(0.0, step=2.0, stop=360.0)) * 86400.0
Δt = t[2] - t[1]

J0 = - 1376.0/4.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p
J = J0 * cos.(2π * (t .- t[1]) / (86400.0*360))

oc = MLMML.OceanColumn(zs);
h_rec = []
b_rec = []

# setup b anomaly


@printf("Importing PyPlot... ")
using PyPlot
@printf("done.\n")
fig, ax = plt[:subplots](1, 1, figsize=(8,6))

for k = 1:length(t)
    println("k = ", k)
    MLMML.stepOceanColumn!(
        oc=oc,
        ua=1.0,
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h)
    push!(b_rec, oc.b_ML)
    ax[:plot]([oc.b_ML, oc.b_ML], [oc.zs[1], oc.zs[1] - oc.h], color="r", markersize=2, marker="+", linewidth=0.5)
    for i=1:length(oc.bs)
#        ax[:plot]([oc.bs, oc.bs], [oc.zs[i], oc.zs[i+1]], color="k", markersize=2, marker="o", linewidth=0.5)
    end
end

ax[:set_ylim]([-20, 1])
plt[:show]()


fig, ax = plt[:subplots](3, 1, figsize=(8,6), sharex=true)

ax[1][:plot](t/86400.0, J)
ax[2][:plot](t/86400.0, b_rec)
ax[3][:plot](t/86400.0, h_rec)

plt[:show]()
