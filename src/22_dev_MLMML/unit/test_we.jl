include("../MLMML.jl")

using Printf
using Statistics: mean
using .MLMML

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

D  = 1000.0
N  = 1001
zs = collect(Float64, range(0.0, stop=-D, length=N))

Δb_init = 0.5 * 10.0 * MLMML.α * MLMML.g  
b_ML_init = 1.0
h_init = 50.0

b_slope = 10.0 / D * MLMML.g * MLMML.α

PERIOD = 360.0 * 86400.0
TOTAL_TIME = .25 * PERIOD

ω = 2π/360.0/86400.0
t = collect(Float64, range(0.0, step=1, stop=90)) * 86400.0
Δt = t[2] - t[1]
t_day = t / 86400.0

J0 = 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p

J = J0 * (t * 0 .+ 1.0)
E = - J0 * t

J = J0 / 90.0 / 86400.0 * t
E = - J .* t / 2

#J = J0 * sin.(ω*t)
#E = J0/ω * (cos.(ω*t) .- 1.0)

U10 = zeros(Float64, length(t))
U10 .= 0.0 #.+ 2.0 * rand(length(U10)) 

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = b_slope,
    b_ML    = b_ML_init,
    h       = h_init,
    Δb      = Δb_init,
    K       = 0.0
)

#if ! MLMML.checkDiffusionStability(oc, Δt=Δt)
#    throw(ErrorException("Stability criteria does not fulfill."))
#end



h_rec = [h_init]
b_rec = [b_ML_init]
hb_rec = [MLMML.getIntegratedBuoyancy(oc)]
we_rec = []
Δb_rec = [oc.b_ML - oc.bs[oc.FLDO]]
bs_rec = zeros(Float64, length(oc.bs), length(t))


bs_rec[:, 1] = oc.bs
for k = 1:length(t)-1
    println("iteration = ", k)
    println("oc.h = ", oc.h)
    info = MLMML.stepOceanColumn!(
        oc=oc,
        ua=U10[k],
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h)
    push!(b_rec, oc.b_ML)
    push!(Δb_rec, info[:Δb])
    push!(hb_rec, MLMML.getIntegratedBuoyancy(oc))#, target_z=-oc.h))
    push!(we_rec, ( info[:flag] == :we ) ? info[:val] : NaN)
    bs_rec[:, k+1] = oc.bs
end


bs_rec_mean = mean(bs_rec, dims=2)
#bs_rec -= repeat(bs_rec_mean, outer=(1, size(bs_rec)[2]))

#ax[:set_ylim]([-20, 1])
#plt[:show]()

plt[:figure]()
for i = 1:5:size(bs_rec)[2]
    plt[:plot](bs_rec[:, i], (oc.zs[1:end-1] + oc.zs[2:end]) / 2.0, "k-")
end


fig, ax = plt[:subplots](5, 1, sharex=true)

fig[:suptitle]("Diagnose of an idealize case")

ax[1][:plot](t_day, - h_rec, label="h")
ax[2][:plot](t_day, hb_rec, "k-", label="hb")
ax[2][:plot](t_day, hb_rec[1] .+ E, "r--", label="hb theory")
ax[3][:plot](t_day[1:end-1], we_rec, label="we")

ax[4][:plot](t_day, Δb_rec, label="Δb")
ax[4][:plot](t_day, ( b_slope/2.0 * (h_rec.^2 .- h_rec[1]^2.0) .+ Δb_init * h_init + E) ./ h_rec, label="Δb theory if trust num model")
#ax[4][:plot](t_day, , label="Δb theory")
ax[5][:plot](t_day, J, label="J")

for a in ax
    a[:legend]()
end

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
#cb = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec, clevs, cmap=cmap, extend="both", zorder=1, antialiased=false)
cb = ax2[:contourf](t_day, (zs[1:end-1] + zs[2:end]) / 2.0, bs_rec, cmap=cmap, extend="both", zorder=1, antialiased=false)
plt[:colorbar](cb, cax=cax)

ax2[:plot](t_day, - h_rec , "r--", linewidth=2, zorder=10)
ax2[:set_ylim]([-500, 0])

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
