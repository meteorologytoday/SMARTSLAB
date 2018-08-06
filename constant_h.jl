using NetCDF

fn = "/surtsey/tienyiah/GFDL_ESM2G_QFlux/SMART_Omon_GFDL-ESM2G_historical_r1i1p1_186101-189112.nc"

rlons    = ncread(fn, "rlon")
rlats    = ncread(fn, "rlat")
vertices = ncread(fn, "vertices")


#TOT_F = ncread(fn, "total_downward_heat_flux")
TOT_F = ncread(fn, "hfds")[:,:,1:36]
SST   = ncread(fn, "tos")[:,:,1:36]


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

#=
# constant mixed-layer depth (MLD)
h_cnst = sum(TOT_F .* dT_dt) / sum(dT_dt .^ 2.0) / (ρ * c_p)
@printf("The optimal mixed-layer thickness: %.2f m\n", h)

# spatial MLD
h_spat = (sum(TOT_F .* dT_dt, 3) ./ sum(dT_dt .^ 2.0, 3) / (ρ * c_p))[:,:,1]
println(size(h_spat))
=#

# temporal-spatial MLD
h_temp_spat = Array{eltype(SST)}(length(rlons), length(rlats), 12)

#for idx in CartesianRange(size(h_temp_spat))
#for idx in CartesianRange((360, 210, 12))
for i = 1:length(rlons), j = 1:length(rlats), m = 1:12
    if j == 1
        println("$i")
    end
    #i, j, m = idx
    h_temp_spat[i, j, m] = sum(TOT_F[i,j,m:12:end] .* dT_dt[i,j,m:12:end]) ./ sum(dT_dt[i,j,m:12:end] .^ 2.0) / (ρ * c_p)
end

# monthly static




time = collect(Float64, 1:12)

filename = "mixed-layer-depth.nc"
varname = "h_temp_spat"

var_attr = Dict(
    "long_name"=>"Mixed-layer Thickness",
    "units"=>"m",
    "coordinates"=>"lat lon"
)

h_temp_spat_attr = Dict(
    "long_name"=>"Monthly Mixed-layer Thickness",
    "units"=>"m",
    "coordinates"=>"time lat lon"
)


rlon_attr = Dict(
    "long_name"=>"longitude in rotated pole grid",
    "standard_name"=> "grid_longitude",
    "units"=>"degrees",
    "axis"=>"X"
)

rlat_attr = Dict(
    "long_name"=>"latitude in rotated pole grid",
    "standard_name"=> "grid_latitude",
    "units"=>"degrees",
    "axis"=>"Y"
)

time_attr = Dict(
    "long_name"=>"Month of the year",
    "units"=>"month"
)


lon_attr  = Dict(
    "long_name"=>"longitude coordinate",
    "standard_name"=>"longitude",
    "units"=>"degrees_east",
    "bounds"=>"lon_vertices"
)

lat_attr  = Dict(
    "long_name"=>"latitude coordinate",
    "standard_name"=>"latitude",
    "units"=>"degrees_north",
    "bounds"=>"lat_vertices"
)

lon_vertices_attr = Dict(
    "units" => "degrees_east"
)


lat_vertices_attr = Dict(
    "units" => "degrees_north"
)

nccreate(
    filename,
    varname,
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    "time", time,  time_attr,
    atts=var_attr,
    mode=NC_CLASSIC_MODEL
)

nccreate(
    filename,
    "lat",
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    atts=lat_attr
)

nccreate(
    filename,
    "lon",
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    atts=lon_attr
)


nccreate(
    filename,
    "lon_vertices",
    "vertices", vertices,
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    atts=lon_vertices_attr
)


nccreate(
    filename,
    "lat_vertices",
    "vertices", vertices, 
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    atts=lat_vertices_attr
)


ncwrite(h_temp_spat, filename, varname)
ncwrite(ncread(fn, "lon"), filename, "lon")
ncwrite(ncread(fn, "lat"), filename, "lat")
ncwrite(time, filename, "time")

ncwrite(rlons, filename, "rlon")
ncwrite(rlats, filename, "rlat")

println(size(ncread(fn, "lon_vertices")))

ncwrite(ncread(fn, "lon_vertices"), filename, "lon_vertices")
ncwrite(ncread(fn, "lat_vertices"), filename, "lat_vertices")
ncwrite(ncread(fn, "vertices")    , filename, "vertices")



ncclose(filename)
