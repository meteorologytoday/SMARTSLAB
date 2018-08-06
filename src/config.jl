using NetCDF

fn = "/surtsey/tienyiah/GFDL_ESM2G_QFlux/SMART_Omon_GFDL-ESM2G_historical_r1i1p1_186101-189112.nc"


rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")
vertices = ncread(fn, "vertices")


TOT_F = ncread(fn, "total_downward_heat_flux")
#TOT_F = ncread(fn, "hfds")[:,:,:]
SST   = ncread(fn, "tos")[:,:,:]

missing_value = ncgetatt(fn, "tos", "missing_value")

TOT_F[TOT_F .== ncgetatt(fn, "total_downward_heat_flux", "missing_value")] = NaN
SST[SST .== ncgetatt(fn, "tos", "missing_value")] = NaN


ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K


mon_secs = 365.0 / 12.0 * 86400.0
nmons = length(ncread(fn, "time"))

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

nyrs = nmons / 12

# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_dt = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) / (2.0 * mon_secs)



#=
function writeData(varnames, varattrs)

    for i in 1:length(varnames)
        varname = varnames[i]
        varattr = varattrs[i]

        


    end
end
=#
