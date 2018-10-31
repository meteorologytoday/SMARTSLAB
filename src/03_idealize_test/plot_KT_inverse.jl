#include("KT_inverse.jl")
include("KT_simulate.jl")



using PyPlot

ms = 5
t = collect(1:12)

fig, ax = plt[:subplots](2, 2, sharex=true, figsize=(12, 8))

ax11_twin = ax[1, 1][:twinx]()
ax21_twin = ax[2, 1][:twinx]()





ax[1,1][:set_xlim](0, length(t_year))
ax[1,2][:set_xlim](0, 13)

#ax[1,1][:set_xlim](0, length(t_year))

ax[1,1][:set_title]("Solar Radiation input.")
ax[1,1][:plot](t_year, S, color="k", label="S")
ax[1,1][:legend]()


ax[2,1][:set_title]("Simulated SST and analytic solution.")
ax[2,1][:plot](t_year, Ts, color="k", label="Ts")
ax[2,1][:plot](t_year, Ts_analytic, color="r", label="Ts analytic", dashes=(5,3))
ax[2,1][:legend]()
#=
ax[1,2][:set_title]("Inversed h")

ax[1,2][:plot](t, omlmax, color="k", marker="o", ms=ms, dashes=(5,2), label="true h")
ax[1,2][:plot](t, data["degenerate"]["h"], color="k", marker="o", ms=ms, dashes=(5,3), label="degenerate")
ax[1,2][:plot](t, data["init_zero"]["h"], color="b", marker="s", ms=ms, label="init_zero")
ax[1,2][:plot](t, data["init_omlmax"]["h"], color="r", marker="d", ms=ms, label="init_omlmax")
ax[1,2][:legend]()
=#

#=
ax[3][:set_title]("Given h and inverted h using linear regression.")
ax[3][:plot](t_year, , color="k", label="h")
ax[3][:plot](t_year, LR_h, color="r", dashes=(5, 3), label="LR_h")
ax[3][:legend]()
=#
plt[:show]()

