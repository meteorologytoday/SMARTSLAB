include("config.jl")

using PyPlot
using DataFrames
using CSV

df = CSV.read(joinpath(data_path, "simulate.jl_SST.csv"))

plot(df[:SST], color="k", linewidth=2.0, linestyle="--")
title(L"Simulated SST at")

show()
