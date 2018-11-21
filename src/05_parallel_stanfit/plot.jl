using PyPlot, PyCall
@pyimport cartopy.crs as ccrs


using NCDatasets

filename = ""

ds = Dataset(filename, "r")



