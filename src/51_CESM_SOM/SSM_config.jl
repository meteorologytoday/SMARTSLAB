wdir = "/home/tienyiah/cesm/scratch/SOM_AQUAP/run"

nlon = 24
nlat = 19

lons = collect(Float64, range(  0, 360, length=nlon+1))[1:end-1]
lats = collect(Float64, range(-90,  90, length=nlat))

zs = collect(Float64, range(0, -500, step=-5))

K = 1e-5

