

regions = Dict(
    "PDO" => ([115.0, 255.0],     [20.0, 65.0]),
    "Nino34" => ([190.0, 240.0], [-5.0,  5.0]),
    "AMO" => ([280.0, 360.0],     [ 0.0, 65.0]),
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

function nanmean(x)
    mask = isfinite.(x)
    return sum(x[mask]) / sum(mask)
end


function region_mean(name, data)
    println("Region:", name)
    mask = region_mask(lon, lat, name)
    time_len = size(x)[3]
    new_x = zeros(eltype(x), time_len)
    
    for i in 1:length(new_x)
        new_x[i] = nanmean(x[:, :, i][mask])
    end

    return new_x
end



for k in keys(regions)
    println(k, " => ", regions[k])
end
