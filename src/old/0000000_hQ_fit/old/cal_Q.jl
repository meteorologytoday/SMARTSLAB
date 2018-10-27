include("config.jl")

using NetCDF

# Discard the first and last year
TOT_F = TOT_F[:, :, 13:end-12] 
dT_dt = (SST[:, :, 14:end-11] - SST[:, :, 12:end-13]) / (2.0 * mon_secs)


println("### Temporal Spatial MLD")

# temporal-spatial Q flux
Qflux_temp_spat = Array{eltype(SST)}(length(rlons), length(rlats), 12)

# Now assume constant h
h = 32.0

for m = 1:12
    println("Doing month [$m]")
    for i = 1:length(rlons), j = 1:length(rlats)
        #i, j, m = idx
        Qflux_temp_spat[i, j, m] = sum(h .* dT_dt[i,j,m:12:end] - TOT_F[i,j,m:12:end]) ./ nyrs
    end
end

Qflux_temp_spat[isnan.(Qflux_temp_spat)] = missing_value

time = collect(Float64, 1:12)

filename = "Q_flux.nc"
varname = "Q_flux"

Qflux_temp_spat_attr = Dict(
    "long_name"=>"Qflux",
    "units"=>"J / m^2 / s",
    "coordinates"=>"time lat lon",
    "missing_value"=>missing_value
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
    "Q_flux",
    "rlon", rlons, rlon_attr,
    "rlat", rlats, rlat_attr,
    "time", time,  time_attr,
    atts=Qflux_temp_spat_attr,
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


ncwrite(Qflux_temp_spat, filename, "Q_flux")
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
