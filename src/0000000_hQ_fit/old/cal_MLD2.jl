
push!(LOAD_PATH, "./src")
import NetCDFHelper
using NetCDF


fn = "/surtsey/tienyiah/GFDL_ESM2G_QFlux/SMART_Omon_GFDL-ESM2G_historical_r1i1p1_186101-189112.nc"

filename = "mixed-layer-depth.nc"
NetCDFHelper.specialCopyNCFile(fn, filename, ["lat", "lon", "lat_vertices", "lon_vertices"])


rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")
vertices = ncread(fn, "vertices")


TOT_F = ncread(fn, "hfds")[:,:,:]
SST   = ncread(fn, "tos")[:,:,:]


ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K


mon_secs = 365.0 / 12.0 * 86400.0
nmons = length(ncread(fn, "time"))

if nmons % 12 != 0
    error("There are $nmons months, not a multiple of 12")
end

# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_dt = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) / (2.0 * mon_secs)


println("### Constant MLD")
# constant mixed-layer depth (MLD)
h_cnst = sum(TOT_F .* dT_dt) / sum(dT_dt .^ 2.0) / (ρ * c_p)
@printf("The optimal mixed-layer thickness: %.2f m\n", h_cnst)


println("### Spatial MLD")
# spatial MLD
h_spat = (sum(TOT_F .* dT_dt, 3) ./ sum(dT_dt .^ 2.0, 3) / (ρ * c_p))[:,:,1]


println("### Temporal Spatial MLD")
# temporal-spatial MLD
h_temp_spat = Array{eltype(SST)}(length(rlons), length(rlats), 12)

#for idx in CartesianRange(size(h_temp_spat))
#for idx in CartesianRange((360, 210, 12))
for m = 1:12
    println("Doing month [$m]")
    for i = 1:length(rlons), j = 1:length(rlats)
        #i, j, m = idx
        h_temp_spat[i, j, m] = sum(TOT_F[i,j,m:12:end] .* dT_dt[i,j,m:12:end]) ./ sum(dT_dt[i,j,m:12:end] .^ 2.0) / (ρ * c_p)
    end
end


time = collect(Float64, 1:12)


varname = "h_temp_spat"

h_cnst_attr = Dict(
    "long_name"=>"Constant Mixed-layer Thickness",
    "units"=>"m"
)

h_spat_attr = Dict(
    "long_name"=>"Mixed-layer Thickness",
    "units"=>"m",
    "coordinates"=>"time lat lon"
)



h_temp_spat_attr = Dict(
    "long_name"=>"Monthly Mixed-layer Thickness",
    "units"=>"m",
    "coordinates"=>"time lat lon"
)


nccreate(
    filename,
    "h_spat",
    "rlon",
    "rlat",
    atts=h_spat_attr,
    mode=NC_CLASSIC_MODEL
)

nccreate(
    filename,
    "h_temp_spat",
    "rlon",
    "rlat",
    "time", time,
    atts=h_temp_spat_attr,
    mode=NC_CLASSIC_MODEL
)

ncwrite(h_temp_spat, filename, "h_temp_spat")
ncwrite(h_spat, filename, "h_spat")

ncclose(filename)
