using PyPlot, PyCall
@pyimport cartopy.crs as ccrs


using NCDatasets

data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
filename = joinpath(data_dir, "HMC_NCAR_5deg_c4_s1000.nc")

ds = Dataset(filename, "r")

proj = ccrs.PlateCarree(central_longitude=0)

ax = plt[:axes](projection=proj)

lon = ds["lon"][:]
lat = ds["lat"][:]
data = ds["Td_mean"][:]

data = convert(Array{Float64}, nomissing(data, NaN)')

ax[:contourf](lon, lat, data, transform=proj)
ax[:coastlines]()

plt[:show]()




