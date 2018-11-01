include("KT_simulate_h_constant.jl")
include("KT_inverse.jl")

using PyPlot

ms = 5
t = collect(1:12)

fig, ax = plt[:subplots](2, 3, figsize=(15, 12))

ax[1,1][:set_xlim](0, 13)
ax[2,1][:set_xlim](0, 13)

ax[1,2][:set_xlim](0, 13)
ax[2,2][:set_xlim](0, 13)
ax[1,3][:set_xlim](0, 13)
ax[2,3][:set_xlim](0, 13)

ax[1,2][:set_xticks](collect(1:12))
ax[2,2][:set_xticks](collect(1:12))
ax[1,3][:set_xticks](collect(1:12))
ax[2,3][:set_xticks](collect(1:12))


#ax[1,1][:set_xlim](0, length(t_year))

ax[1,1][:set_title]("Solar Radiation input [\$\\mathrm{W} \\, \\mathrm{m}^{-2}\$]")
ax[1,1][:set_xlabel]("year")
ax[1,1][:plot](t, avg_S, color="k", label="S")
ax[1,1][:legend]()


ax[2,1][:set_title]("Simulated SST and analytic solution [\$\\mathrm{K}\$]")
ax[2,1][:set_xlabel]("year")
ax[2,1][:plot](t, avg_Ts, color="k", label="Ts")
ax[2,1][:plot](t, avg_Ts_analytic, color="r", label="Ts analytic", dashes=(5,3))
ax[2,1][:legend]()

ax[1,2][:set_title]("MLD [\$\\mathrm{m}\$]")
ax[2,1][:set_xlabel]("month")
ax[1,2][:plot](t, data["degenerate"]["h"], color="k", marker="o", ms=ms, label="degenerate")
ax[1,2][:plot](t, data["init_zero"]["h"], color="b", marker="s", ms=ms, label="init_zero")
ax[1,2][:plot](t, data["init_omlmax"]["h"], color="r", marker="d", ms=ms, label="init_omlmax")
ax[1,2][:plot](t, omlmax, color="k", marker="o", ms=ms, dashes=(5,2), label="true h")
ax[1,2][:legend]()


ax[2,3][:set_title]("Inversed Q [\$\\mathrm{W} \\, \\mathrm{m}^{-2}\$]")
ax[2,3][:set_xlabel]("month")
ax[2,3][:plot](t, data["degenerate"]["Q"], color="k", marker="o", ms=ms, label="degenerate")
ax[2,3][:plot](t, data["init_zero"]["Q"], color="b", marker="s", ms=ms, label="init_zero")
ax[2,3][:plot](t, data["init_omlmax"]["Q"], color="r", marker="d", ms=ms, label="init_omlmax")
ax[2,3][:legend]()

ax[1,3][:set_title]("MLD [\$\\mathrm{m}\$]\n(magnified for degenerate and init_zero scenarios)")
ax[1,3][:set_xlabel]("month")
ax[1,3][:plot](t, data["degenerate"]["h"], color="k", marker="o", ms=ms, dashes=(5,3), label="degenerate")
ax[1,3][:plot](t, data["init_zero"]["h"], color="b", marker="s", ms=ms, label="init_zero")
ax[1,3][:legend]()

ax[2,2][:set_visible](false)

plt[:show]()
fig[:savefig](joinpath("img", "idealized_KT_inverse_h_constant.png"), dpi=100)


