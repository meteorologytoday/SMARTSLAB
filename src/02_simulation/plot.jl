include("config.jl")

@printf("Running %s\n", basename(@__FILE__))

using PyPlot
using DataFrames
using CSV

df = CSV.read(joinpath(data_path, "simulate.jl_SST.csv"))

fig, ax = plt[:subplots](3, 1, sharex=true, figsize=(20, 10))

ax[1][:set_xlim](1,240)
ax[1][:set_ylim](285,310)

t = collect(1:length(df[:SST_REAL]))

ax[1][:set_title]("SST budget")
ax[1][:plot](t, df[:SST_REAL], color="r", linewidth=2.0, linestyle="-", label="Real")
ax[1][:plot](t, df[:SST_RK4], color="g", linewidth=2.0, dashes=(7,3), label="RK4")
ax[1][:plot](t, df[:SST_EULER], color="b", linewidth=2.0, dashes=(7,3,3,3), label="EULER")

ax[1][:legend]()

ax[2][:set_title]("F, Q, h")
ax[2][:plot](t, df[:F], color="k", linewidth=2.0, linestyle="-", label="F")
ax[2][:plot](t, df[:Q] ./ 100.0, color="k", linewidth=2.0, linestyle="--", label="Q/100")

ax2_twin = ax[2][:twinx]()
ax2_twin[:plot](t, df[:h], color="k", linewidth=2.0, linestyle="-.", label="h")

ax[2][:legend]()


ax[3][:set_title]("F/h, Q/h, - T/h * dh/dt")
ax[3][:plot](t, df[:F] ./ df[:h], color="k", linewidth=2.0, linestyle="-", label="F/h")
ax[3][:plot](t, df[:Q] ./ df[:h], color="k", linewidth=2.0, linestyle="--", label="Q/h")
ax[3][:plot](t, - df[:SST_REAL] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="r", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (Real)")
ax[3][:plot](t, - df[:SST_RK4] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="g", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (RK4)")
ax[3][:plot](t, - df[:SST_EULER] ./ df[:h] .* df[:dhdt] .* (ρ*c_p), color="b", linewidth=2.0, linestyle="-.", label="T/h * dh/dt (EULER)")
ax[3][:legend]()


#plot(df[:SST_simulated], color="r", linewidth=2.0, linestyle="-")


plt[:show]()

fig[:savefig]("")
