using NetCDF

include("config.jl")

rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")

lons    = ncread(fn, "lon")
lats    = ncread(fn, "lat")

TOT_F = ncread(fn, "total_downward_heat_flux")
#TOT_F = ncread(fn, "hfds")[:,:,:]
SST   = ncread(fn, "tos")[:,:,:]
tp = eltype(TOT_F)

missing_value = ncgetatt(fn, "tos", "missing_value")

TOT_F[TOT_F .== ncgetatt(fn, "total_downward_heat_flux", "missing_value")] = NaN
SST[SST .== ncgetatt(fn, "tos", "missing_value")] = NaN

spatial_mask = isnan.(SST[:,:,1])
spatial_temporal_mask = isnan.(SST)

T_star = SST * ρ * c_p


nmons = length(ncread(fn, "time"))

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = nmons / 12
@printf("We have %02d years of data.\n", nyrs)
