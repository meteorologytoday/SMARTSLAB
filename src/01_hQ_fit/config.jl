using NetCDF

@printf("Running %s\n", basename(@__FILE__))

function prtArr(A)
    for i = 1:size(A)[1]
        for j = 1:size(A)[2]

            @printf("%d ", A[i,j]) 

        end
        @printf("\n")
    end
end

ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K

fn = "/home/tienyiah/projects/SMART/data/SMART_Omon_GFDL-ESM2G_historical_r1i1p1_186101-200112.nc"


rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")


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

mon_secs = 365.0 / 12.0 * 86400.0
nmons = length(ncread(fn, "time"))

dt = 1.0 * mon_secs

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = nmons / 12
@printf("We have %02d years of data.\n", nyrs)