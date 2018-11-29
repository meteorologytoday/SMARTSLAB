using PyCall
@pyimport matplotlib as mpl
mpl.use("Agg")

using PyPlot

@pyimport mpl_toolkits.basemap as basemap
@pyimport mpl_toolkits.basemap.cm as cm

using Formatting

include("config_plot.jl")

clevs_h  = collect(range(0; stop=200, length=11))
clevs_Q  = collect(range(-100; stop=100, length=11))
clevs_Td = collect(range(-10; stop=15, length=11))

llon = repeat(lon ; outer=(1, length(lat)))
llat = repeat(lat'; outer=(length(lon), 1))

### plot Td ###
println("Plotting Td")

# setting maps
filename = normpath(joinpath(
    img_dir,
    format("{}-maps-Td.png", splitext(basename(nc_filename))[1])
))

fig, ax = plt[:subplots](1, 1, figsize=(8, 6))

lab_lons = collect(  0:60:360)
lab_lats = collect(-90:30: 90)

m = basemap.Basemap(ax=ax,  projection="merc",llcrnrlat=-80,urcrnrlat=80,
         llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution="c",lon_0=180.0)
m[:drawcoastlines](linewidth=1.5)
#m[:fillcontinents](color="#888888", lake_color="aqua")
m[:drawmeridians](lab_lons, labels=repeat([true, ], outer=(length(lab_lons),)))
m[:drawparallels](lab_lats, labels=repeat([true, ], outer=(length(lab_lats),)))

x, y = m(llon, llat)

#cmap = plt[:get_cmap]("gist_ncar_r")
#cmap = plt[:get_cmap]("GnBu")
#cmap = plt[:get_cmap]("gnuplot2_r")
cmap = plt[:get_cmap]("bwr")
cs = m[:contourf](x, y, data_Td[:, :], clevs_Td; cmap=cmap, extend="both")
cbar = m[:colorbar](cs, location="bottom", pad="8%", ticks=clevs_Td)
cbar[:set_label]("Deep Ocean Temperature \$T_\\mathrm{d}\$ [\$ {}^\\circ\\mathrm{C}\$]")

fig[:suptitle](format("File: {}", basename(nc_filename)))
plt[:show]()
fig[:savefig](filename, dpi=100)
println("Output file:", filename)
plt[:close](fig)

### plot h and Q ###
for i = 1:size(data_h)[3]
    println("Plotting month ", i)
    # setting maps
    filename = normpath(joinpath(
        img_dir,
        format("{}-maps-{:02d}.png", splitext(basename(nc_filename))[1], i)
    ))

    fig, axes = plt[:subplots](1, 2, figsize=(12, 6))
    mm = []

    lab_lons = collect(  0:60:360)
    lab_lats = collect(-90:30: 90)

    for ax in axes
        m = basemap.Basemap(ax=ax,  projection="merc",llcrnrlat=-80,urcrnrlat=80,
             llcrnrlon=0,urcrnrlon=360,lat_ts=20,resolution="c",lon_0=180.0)
        push!(mm, m)

        m[:drawcoastlines](linewidth=1.5)
        #m[:fillcontinents](color="#888888", lake_color="aqua")
        m[:drawmeridians](lab_lons, labels=repeat([true, ], outer=(length(lab_lons),)))
        m[:drawparallels](lab_lats, labels=repeat([true, ], outer=(length(lab_lats),)))
        
    end

    x, y = mm[1](llon, llat)

    cmap = plt[:get_cmap]("gist_ncar_r")
    cmap = plt[:get_cmap]("GnBu")
    cmap = plt[:get_cmap]("gnuplot2_r")
    cs = mm[1][:contourf](x, y, data_h[:, :, i], clevs_h; cmap=cmap, extend="both")
    cbar = mm[1][:colorbar](cs, location="bottom", pad="8%", ticks=clevs_h)
    cbar[:set_label]("MLD [\$\\mathrm{m}\$]")

    cmap = plt[:get_cmap]("Spectral_r")
    cs = mm[2][:contourf](x, y, data_Q[:, :, i], clevs_Q; cmap=cmap, extend="both")
    cbar = mm[2][:colorbar](cs, location="bottom", pad="8%", ticks=clevs_Q)
    cbar[:set_label]("Q flux [\$\\mathrm{W} \\, \\mathrm{m}^{-2} \$]")

    fig[:suptitle](format("File: {} \n Month: {:02d}", basename(nc_filename), i))


    #plt[:show]()
    fig[:savefig](filename, dpi=100)
    println("Output file:", filename)
    plt[:close](fig)
end
