

regions = Dict(
    "NPAC-1" => ([150.0, 225.0], [ 10.0,  25.0]),
    "NPAC-2" => ([150.0, 225.0], [ 25.0,  50.0]),
    "EPAC"   => ([215.0, 270.0], [-10.0,  10.0]),
    "WPAC"   => ([135.0, 215.0], [-10.0,  10.0]),
    "SPAC-1" => ([170.0, 255.0], [-25.0, -10.0]),
    "SPAC-2" => ([170.0, 255.0], [-50.0, -25.0]),
    "NATL-1" => ([290.0, 340.0], [ 10.0,  25.0]),
    "NATL-2" => ([290.0, 340.0], [ 25.0,  50.0]),
    "MATL"   => ([325.0, 350.0], [-10.0,  10.0]),
    "SATL-1" => ([325.0, 360.0], [-25.0, -10.0]),
    "SATL-2" => ([325.0, 360.0], [-50.0, -25.0]),
    "NIND"   => ([ 60.0,  90.0], [ 10.0,  15.0]),
    "MIND"   => ([ 60.0,  90.0], [-10.0,  10.0]),
    "SIND-1" => ([ 60.0,  90.0], [-25.0, -10.0]),
    "SIND-2" => ([ 60.0,  90.0], [-50.0, -25.0]),
)


function in_region(name, lat, lon)
    lon = mod(lon, 360.0)
    lon_rng, lat_rng = regions[name]
    return lon_rng[1] <= lon && lon < lon_rng[2] && lat_rng[1] <= lat && lat < lat_rng[2] 
end

function meshgrid(x, y)
    #=
    xx = zeros(eltype(x), length(x), length(y))
    yy = copy(xx)

    for i = 1:length(x), j = 1:length(y)
        xx[i, j] = x[i]
        yy[i, j] = y[i]
    end
    =#

    xx = repeat(x, outer=(1, length(y)))
    yy = repeat(y, outer=(1, length(x)))'

    return xx, yy
end


function region_mask(lon, lat, name)
    llon,llat = meshgrid(lon, lat)
    lon_rng, lat_rng = regions[name]

    return (lon_rng[1] .<= llon      ) .& 
           (llon       .<  lon_rng[2]) .&
           (lat_rng[1] .<= llat      ) .&
           (llat       .< lat_rng[2] )
end

for k in keys(regions)
    println(k, " => ", regions[k])
end
