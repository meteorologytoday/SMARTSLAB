include("single_point.jl")
#include("config.jl")

#=
#lat_i = 137
lon_i = 125
#idx = (lon_i, lat_i, :)

S     = ncread(fn, "S")[idx] / ρ / c_p
B     = ncread(fn, "B")[idx] / ρ / c_p
SST   = ncread(fn, "tos")[idx]
=#

ret_len = 12 * 10


h_init  = 30.0
Ts_init = SST[1]
Td      = 273.15
β       = Inf
S       = 

h, Ts = KrausTurnerMLD_simplified(
    h_init,
    Ts_init,
    Td,
    dt,
    S,
    dSdt,
    B,
    G,
    D,
    ret_len 
)
