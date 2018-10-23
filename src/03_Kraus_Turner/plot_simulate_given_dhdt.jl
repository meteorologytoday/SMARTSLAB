include("simulate_given_dhdt.jl")
using PyPlot

fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(12, 8))

ax1_twin = ax[1][:twinx]()
ax2_twin = ax[2][:twinx]()

#ax[1][:set_title]("S")
ax[1][:plot](t_year, S, color="k", label="S")
ax[1][:legend]()


ax[2][:plot](t_year, Ts, color="k", label="Ts")
ax[2][:plot](t_year, Ts_analytic, color="r", label="Ts analytic", dashes=(5,3))
ax[2][:legend]()

ax[3][:plot](t_year, h, color="k", label="h")
ax[3][:legend]()

plt[:show]()

