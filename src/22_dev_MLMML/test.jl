include("MLMML.jl")

using Printf
using .MLMML

D  = 1000.0
N  = 5001
zs = collect(Float64, range(0.0, stop=-D, length=N))
bg_Ts = 30.0 * (1.0 .+ (zs[1:end-1] + zs[2:end]) / 2D)

t = collect(Float64, range(0.0, step=1.0, stop=360.0*10)) * 86400.0
Δt = t[2] - t[1]

J0 = - 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p
J = J0 * cos.(2π * (t .- t[1]) / (86400.0*360))

U10 = zeros(Float64, length(t))
U10 .= 1.0 .+ 2.0 * rand(length(U10)) 

oc = MLMML.OceanColumn(zs);
h_rec = []
b_rec = []

bs_rec = zeros(Float64, length(oc.bs), length(t))

bg_bs = Ts * MLMML.α * MLMML.g

oc.bs += bg_bs
oc.b_ML = oc.bs[1]

MLMML.correctOceanColumn!(oc)
# setup b anomaly


@printf("Importing PyPlot... ")
using PyPlot
@printf("done.\n")
#fig, ax = plt[:subplots](1, 1, figsize=(8,6))

bs_rec[:, 1] = oc.bs
for k = 1:length(t)-1
    println("iteration = ", k)
    MLMML.stepOceanColumn!(
        oc=oc,
        ua=U10[k],
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h)
    push!(b_rec, oc.b_ML)
    bs_rec[:, k+1] = oc.bs - bg_bs
    #ax[:plot]([oc.b_ML, oc.b_ML], [oc.zs[1], oc.zs[1] - oc.h], color="r", markersize=2, marker="+", linewidth=0.5)
    for i=1:length(oc.bs)
#        ax[:plot]([oc.bs, oc.bs], [oc.zs[i], oc.zs[i+1]], color="k", markersize=2, marker="o", linewidth=0.5)
    end
end

#ax[:set_ylim]([-20, 1])
#plt[:show]()


fig, ax = plt[:subplots](1, 1, figsize=(12, 6))
t_day = t / 86400.0
cmap = plt[:get_cmap]("jet")
clevs = range(0.0, stop=0.01, length=20) |> collect
cb = ax[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec, clevs, cmap=cmap, extend="both", zorder=1)
plt[:colorbar](cb)

ax[:plot](t_day[2:end], - h_rec , "r--", linewidth=2, zorder=10)
ax[:set_ylim]([-1000, 0])
plt[:show]()

#=
fig, ax = plt[:subplots](3, 1, figsize=(8,6), sharex=true)

ax[1][:plot](t/86400.0, J, label="- J")
ax[1][:plot]([t[1], t[end]]/86400.0, [0, 0], "k--")
ax[2][:plot](t/86400.0, b_rec, label="b")
ax[3][:plot](t/86400.0, h_rec, label="h")

ax[1][:legend]()
ax[2][:legend]()
ax[3][:legend]()

plt[:show]()
=#
