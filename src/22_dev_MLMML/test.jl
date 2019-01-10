include("MLMML.jl")

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

using Printf
using Statistics: mean
using .MLMML

D  = 5000.0
N  = 5001
zs = collect(Float64, range(0.0, stop=-D, length=N))


PERIOD = 360.0 * 86400.0
TOTAL_TIME = 0.25 * PERIOD

t = collect(Float64, range(0.0, step=1 * 86400.0, stop=TOTAL_TIME))
Δt = t[2] - t[1]
t_day = t / 86400.0

J0 = - 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p
J = J0 * (-sin.(2π * (t .- t[1]) / (86400.0*360)))

U10 = zeros(Float64, length(t))
U10 .= 0.0 #.+ 2.0 * rand(length(U10)) 

oc = MLMML.OceanColumn(zs);
h_rec = []
b_rec = []

bs_rec = zeros(Float64, length(oc.bs), length(t))

bs = copy(oc.bs)
b_ML_init = 1.0
h_init = MLMML.h_min
b_slope = 30.0 / D * MLMML.g * MLMML.α
for i = 1:length(bs)
    z = (oc.zs[i] + oc.zs[i+1]) / 2.0
    if z > -h_init
        bs[i] = b_ML_init
    else
        bs[i] = b_ML_init - b_slope * (-z - h_init)
    end
end

MLMML.setBuoyancy!(oc, bs=bs, b_ML=b_ML_init, h=h_init)

plt[:figure]()
plt[:plot](oc.bs, (oc.zs[1:end-1] + oc.zs[2:end]) / 2.0, "k-")

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
    bs_rec[:, k+1] = oc.bs
end


bs_rec_mean = mean(bs_rec, dims=2)
#bs_rec -= repeat(bs_rec_mean, outer=(1, size(bs_rec)[2]))

#ax[:set_ylim]([-20, 1])
#plt[:show]()



gs0 = GS.GridSpec(1, 2, width_ratios=[100,5])
gs_l = GS.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[1], height_ratios=[1, 4])

fig = plt[:figure](figsize=(12, 6))
ax1 = plt[:subplot](gs_l[1])
ax2 = plt[:subplot](gs_l[2])
cax = plt[:subplot](gs0[2])




ax1[:plot](t_day, J, label="J")
ax1[:plot]([t[1], t[end]]/86400.0, [0, 0], "k--")

cmap = plt[:get_cmap]("jet")
clevs = range(0.07, stop=0.09, length=20) |> collect
cb = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
plt[:colorbar](cb, cax=cax)

ax2[:plot](t_day[2:end], - h_rec , "r--", linewidth=2, zorder=10)
ax2[:set_ylim]([-200, 0])

tlim = [t_day[1], t_day[end]]
ax1[:set_xlim](tlim)
ax2[:set_xlim](tlim)


ax1[:set_ylabel]("Insolation Flux")
ax2[:set_ylabel]("Z [m]")

ax1[:set_xticks](collect(range(0.0, step=PERIOD/4.0, stop=TOTAL_TIME))/86400.0)  
ax2[:set_xticks](collect(range(0.0, step=PERIOD/4.0, stop=TOTAL_TIME))/86400.0)  

using Formatting
fig[:suptitle](format("N = {:d}", N))
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
