using PyPlot
using Printf
using NetCDF

@printf("Running %s\n", basename(@__FILE__))

data_path = joinpath(dirname(@__FILE__), "..", "..", "data")

fn_KT = joinpath(data_path, "GDM_KT_fit.jl.nc")
fn_HQ = joinpath(data_path, "GDM_HQ_fit.jl.nc")
fn_HQ_REAL = joinpath(data_path, "case3_hQ.jl.nc")


rlons = ncread(fn_KT, "rlon")
rlats = ncread(fn_KT, "rlat")
ts    = collect(1:12)

rlon_i = 68
rlat_i = 124


h_KT = ncread(fn_KT, "h")[rlon_i, rlat_i, :]
Q_KT = ncread(fn_KT, "Q")[rlon_i, rlat_i, :]

h_HQ = ncread(fn_HQ, "h")[rlon_i, rlat_i, :]
Q_HQ = ncread(fn_HQ, "Q")[rlon_i, rlat_i, :]

h_HQ_REAL = ncread(fn_HQ_REAL, "h")[rlon_i, rlat_i, :]
Q_HQ_REAL = ncread(fn_HQ_REAL, "Q")[rlon_i, rlat_i, :]


fig, ax = plt[:subplots](2, 1, sharex=true, figsize=(12, 8))

ax[1][:set_ylim](-100, 100)


ax[1][:plot](ts, h_KT, label="KT", color="r")
ax[1][:plot](ts, h_HQ, label="HQ", color="b")
ax[1][:plot](ts, h_HQ_REAL, label="HQ_REAL", color="k", dashes=(7,3))
ax[1][:legend]()

ax[2][:plot](ts, Q_KT, label="KT", color="r")
ax[2][:plot](ts, Q_HQ, label="HQ", color="b")
ax[2][:plot](ts, Q_HQ_REAL, label="HQ_REAL", color="k", dashes=(7,3))
ax[2][:legend]()

plt[:show]()
