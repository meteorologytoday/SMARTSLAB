include("../MLMML.jl")
include("genLineCoord.jl")

using Printf
using Statistics: mean
using .MLMML
using Formatting

@printf("Importing PyPlot... ")
using PyCall
using PyPlot
@pyimport matplotlib.gridspec as GS
@printf("done.\n")

D  = 1000.0
N  = 1001
zs = collect(Float64, range(0.0, stop=-D, length=N))

Δb_0 = 0.5 * 10.0 * MLMML.α * MLMML.g  
b_ML_0 = 1.0
h_ML_0 = 50.0

b_slope = 2.0 / 4000.0 * MLMML.g * MLMML.α

PERIOD = 360.0 * 86400.0
TOTAL_TIME = .25 * PERIOD

ω = 2π/360.0/86400.0
t = collect(Float64, range(0.0, step=1, stop=90)) * 86400.0
Δt = t[2] - t[1]
t_day = t / 86400.0

U10 = zeros(Float64, length(t))
U10 .= 0.0 #.+ 2.0 * rand(length(U10)) 

J0 = 100.0 * MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p

J = J0 * (t * 0 .+ 1.0)
E = - J0 * t


J_func(t) = J0 / 90.0 / 86400.0 * t 
E_func(t) = - J_func(t) * t/2.0

J = J_func.(t)
E = E_func.(t)

n = 0.2
# Numerical solution (should be almost analytic)
using DifferentialEquations
f(h_ML, p, t) = h_ML * n * J_func(t) / ( b_slope * (h_ML^2.0 - h_ML_0^2.0) / 2.0 + MLMML.getTKE(fric_u=MLMML.getFricU(ua=0.0)) + Δb_0 * h_ML_0 + E_func(t))
prob = ODEProblem(f, h_ML_0, (t[1], t[end]))
#num_h_ML     = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
num_h_ML     = solve(prob, reltol=1e-5)
num_dh_ML_dt = [f(num_h_ML.u[i], 0 , num_h_ML.t[i]) for i = 1:length(num_h_ML.u)]


#J = J0 * sin.(ω*t)
#E = J0/ω * (cos.(ω*t) .- 1.0)

oc = MLMML.makeSimpleOceanColumn(
    zs      = zs,
    b_slope = b_slope,
    b_ML    = b_ML_0,
    h_ML    = h_ML_0,
    Δb      = Δb_0,
    K       = 0.0
)

#if ! MLMML.checkDiffusionStability(oc, Δt=Δt)
#    throw(ErrorException("Stability criteria does not fulfill."))
#end



h_rec = [h_ML_0]
b_rec = [b_ML_0]
hb_rec = [MLMML.getIntegratedBuoyancy(oc)]
we_rec = []
Δb_rec = [oc.b_ML - oc.bs[oc.FLDO]]
bs_rec = zeros(Float64, length(oc.bs), length(t))
FLDO_rec = [oc.FLDO]

bs_rec[:, 1] = oc.bs
for k = 1:length(t)-1
    println("iteration = ", k)
    println("oc.h_ML = ", oc.h_ML)
    info = MLMML.stepOceanColumn!(
        oc=oc,
        ua=U10[k],
        B0=0.0,
        J0=J[k],
        Δt=Δt,
    )
    push!(h_rec, oc.h_ML)
    push!(b_rec, oc.b_ML)
    push!(Δb_rec, info[:Δb])
    push!(hb_rec, MLMML.getIntegratedBuoyancy(oc))
    push!(we_rec, ( info[:flag] == :we ) ? info[:val] : NaN)
    push!(FLDO_rec, oc.FLDO)
    bs_rec[:, k+1] = oc.bs
end


plt[:figure]()
for i = 1:5:size(bs_rec)[2]
    plt[:plot](bs_rec[:, i], (oc.zs[1:end-1] + oc.zs[2:end]) / 2.0, "k-")
end


gs0 = GS.GridSpec(1, 2, width_ratios=[1,1])
gs_l = GS.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[1])

fig = plt[:figure](figsize=(16, 12))

ax = []
for i = 1:5
    push!(ax, plt[:subplot](gs_l[i]))
end

bs_ax = plt[:subplot](gs0[2])


fig[:suptitle](
    "Diagnose of an idealize case (Linear Forcing)\n"
    * "Δb = " * format("{:.2e}", Δb_0) * " \$\\mathrm{m}\\,\\mathrm{s}^{-2}\$, "
    * "b_slope = " * format("{:.2e}", b_slope) * " \$\\mathrm{s}^{-2}\$"
)

ax[1][:plot](t_day, J * 1e8, "k-", label="\$J_0\$ surface forcing")
ax[2][:plot](t_day, - h_rec, "k-", label="h_ML simulated")
ax[2][:plot](num_h_ML.t/86400.0, - num_h_ML.u, "r--", label="h_ML in theory (numerically solve diff eqn)")

ax[3][:plot](t_day, hb_rec, "k-", label="\$\\int b\\,\\mathrm{d} z\$ simulated")
ax[3][:plot](t_day, hb_rec[1] .+ E, "r--", label="\$\\int b\\,\\mathrm{d} z\$ in theory (analytic solution)")
ax[4][:plot](t_day[1:end-1], we_rec * 1e5, "k-", label="\$w_e\$ simulated")
ax[4][:plot](num_h_ML.t/86400, num_dh_ML_dt * 1e5, "r--", label="\$w_e\$ (evaluate \"h_ML in theory\" into analytic eqn of \$w_e\$)")

ax[5][:plot](t_day, Δb_rec, "k-", label="Δb simulated")
ax[5][:plot](t_day, ( b_slope/2.0 * (h_rec.^2 .- h_rec[1]^2.0) .+ Δb_0 * h_ML_0 + E) ./ h_rec, "r--", label="\$\\Delta b\$ in theory (evaluate \"h_ML in theory\" into analytic eqn of \$\\Delta b\$)")
#ax[4][:plot](t_day, , label="Δb theory")


ax[1][:set_ylabel]("Buoyancy Flux\n[\$\\times\\,10^{-8}\\,\\mathrm{m}\\,\\mathrm{s}^{-3}\$]")
ax[2][:set_ylabel]("Z\n[m]")
ax[3][:set_ylabel]("Integrated\nBuoyancy\n[\$\\mathrm{m}^2 \\, \\mathrm{s}^{-2}\$]")
ax[4][:set_ylabel]("\$w_e\$\n[\$\\times\\,10^{-5}\\,\\mathrm{m}\\,\\mathrm{s}^{-1}\$]")
ax[5][:set_ylabel]("Δb\n[\$\\mathrm{m}\\,\\mathrm{s}^{-2}\$]")


for a in ax
    a[:legend]()
end

for i = 1:10:size(bs_rec)[2]
    println(i)
    x, z = genLineCoord(b_ML=b_rec[i], h_ML=h_rec[i], zs=oc.zs, bs=bs_rec[:,i], FLDO=FLDO_rec[i])
    if i == 1
        bs_ax[:plot](x, z, "k-", label="Init")
    else
        bs_ax[:plot](x, z, "--")
    end
    bs_ax[:text](b_rec[i], 5, "$(i-1)", va="bottom", ha="center")
end
bs_ax[:legend]()

bs_ax[:set_ylabel]("Z [m]")
bs_ax[:set_xlabel]("Buoyancy [\$\\mathrm{m}\\,\\mathrm{s}^{-2}\$]")
bs_ax[:set_ylim]([-100,10])
bs_ax[:set_xlim]([.98,1.005])



# Hovmoeller Diagram
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
