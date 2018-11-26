using PyPlot, PyCall
@pyimport mpl_toolkits.basemap as basemap
@pyimport mpl_toolkits.basemap.cm as cm

using NCDatasets

cvt = x -> convert(Array{Float64}, nomissing(x, NaN))


# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
filename = joinpath(data_dir, "HMC_NCAR_5deg_c4_s1000.nc")

ds = Dataset(filename, "r")

lon = cvt(ds["lon"][:])
lat = cvt(ds["lat"][:])
data = cvt(ds["Td_mean"][:])

clevs = collect(range(260; stop=300, length=9))

llon = repeat(lon ; outer=(1, length(lat)))
llat = repeat(lat'; outer=(length(lon), 1))


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


x, y = m(llon, llat)

cs = m[:contourf](x, y, data, clevs; cmap=plt[:get_cmap]("Spectral_r"))
cbar = m[:colorbar](cs, location="bottom", pad="5%", ticks=clevs)

plt[:show]()




