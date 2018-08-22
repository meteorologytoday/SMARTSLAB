include("config.jl")

using PyPlot
using DataFrames
using CSV

df = CSV.read(joinpath(data_path, "simulate.jl_SST.csv"))

plot(df[:SST_real], color="k", linewidth=2.0, linestyle="-")
plot(df[:SST_simulated], color="r", linewidth=2.0, linestyle="-")
title("Simulated SST at")

show()
