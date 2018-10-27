include("regions.jl")


lat = collect(range(-90, stop=90, length=50))
lon = collect(range(0, stop=360, length=100))

using PyPlot, PyCall
@pyimport mpl_toolkits.basemap as basemap

m = basemap.Basemap(projection="merc",llcrnrlat=-80,urcrnrlat=80,
            llcrnrlon=0,urcrnrlon=359,lat_ts=20,resolution="c",lon_0=180.0)

m[:drawcoastlines](linewidth=0.5)
#m[:fillcontinents](color="#ffffff", lake_color="aqua")
m[:drawmeridians](collect(0:30:360))
m[:drawparallels](collect(-90:30:90))

m[:drawmapboundary](fill_color="aqua")

#m[:scatter](190, 0, color="r", s=200, marker="o", latlon=true)

for k in keys(regions)
    for i=1:length(lat), j=1:length(lon)
        if in_region(k, lat[i], lon[j])
            
            m[:scatter](m(lon[j], lat[i])... ,color="r", s=5, marker="o")
        end
    end
end




plt[:show]()

#=
x, y = map(rad2deg(lons), rad2deg(lats))

        # Contour data over the map.
cs = map[:contour](x, y, wave+mean, 15, latsinewidths=1.5)
=#
