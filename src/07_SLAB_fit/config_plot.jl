using NCDatasets
using Formatting

cvt = x -> convert(Array{Float64}, nomissing(x, NaN))
stat = x -> (mean(x, dims=(1,))[1,:, :], std(x, dims=(1,))[1,:,:])



# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
img_dir  = normpath(joinpath( dirname(@__FILE__), "..", "..", "img"))
nc_filename = joinpath(data_dir, "fit_SLAB_bayesian.nc")

ds = Dataset(nc_filename, "r")

lon = cvt(ds["lon"][:])
lat = cvt(ds["lat"][:])
data_h = cvt(ds["h"][:])
data_Q = cvt(ds["Q"][:])

close(ds)


