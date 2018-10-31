using PyPlot, PyCall
@pyimport mpl_toolkits.basemap as basemap

include("general_config.jl")
include("regions.jl")

lat = collect(range(-90, stop=90, length=50))
lon = collect(range(0, stop=360, length=100))


plt[:figure](figsize=(12,10))
m = basemap.Basemap(projection="merc",llcrnrlat=-80,urcrnrlat=80,
            llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution="c",lon_0=180.0)


lab_lons = collect(0:30:360)
lab_lats = collect(-90:30:90)


m[:drawcoastlines](linewidth=1.5)
m[:fillcontinents](color="#888888", lake_color="aqua")
m[:drawmeridians](lab_lons, labels=repeat([true, ], outer=(length(lab_lons),)))
m[:drawparallels](lab_lats, labels=repeat([true, ], outer=(length(lab_lats),)))

#m[:drawmapboundary](fill_color="aqua")
#m[:scatter](190, 0, color="r", s=200, marker="o", latlon=true)

for k in keys(regions)
    println("Region: ", k)
    v = regions[k]
    r = []
    push!(r, m(v[1][1], v[2][1]))
    push!(r, m(v[1][1], v[2][2]))
    push!(r, m(v[1][2], v[2][2]))
    push!(r, m(v[1][2], v[2][1]))

    m[:plot]([r[1][1], r[2][1], r[3][1], r[4][1], r[1][1]], [r[1][2], r[2][2], r[3][2], r[4][2], r[1][2]], linewidth=1.0)
    plt[:text]((r[1][1] + r[3][1])/2.0, (r[1][2] + r[2][2])/2.0, k, va="center", ha="center")

end

plt[:savefig](joinpath(img_path, "regions.png"), dpi=101)
plt[:show]()
