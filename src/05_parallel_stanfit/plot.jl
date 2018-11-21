using PyPlot, PyCall
@pyimport mpl_toolkits.basemap as basemap

using NCDatasets

# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
filename = joinpath(data_dir, "HMC_NCAR_5deg_c4_s1000.nc")

ds = Dataset(filename, "r")

proj = ccrs.PlateCarree(central_longitude=0)

ax = plt[:axes](projection=proj)

lon = ds["lon"][:]
lat = ds["lat"][:]
data = ds["Td_mean"][:]

data = convert(Array{Float64}, nomissing(data, NaN))


# setting maps
plt[:figure](figsize=(12,10))
m = basemap.Basemap(projection="merc",llcrnrlat=-80,urcrnrlat=80,
            llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution="c",lon_0=180.0)


lab_lons = collect(0:30:360)
lab_lats = collect(-90:30:90)


m[:drawcoastlines](linewidth=1.5)
#m[:fillcontinents](color="#888888", lake_color="aqua")
m[:drawmeridians](lab_lons, labels=repeat([true, ], outer=(length(lab_lons),)))
m[:drawparallels](lab_lats, labels=repeat([true, ], outer=(length(lab_lats),)))


x, y = map(lon, lat)

m[:contourf](x, y, data)

plt[:show]()




