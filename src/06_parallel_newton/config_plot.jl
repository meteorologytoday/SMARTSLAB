using NCDatasets
using Formatting

cvt = x -> convert(Array{Float64}, nomissing(x, NaN))
stat = x -> (mean(x, dims=(1,))[1,:, :], std(x, dims=(1,))[1,:,:])



# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
img_dir  = normpath(joinpath( dirname(@__FILE__), "..", "..", "img"))
nc_filename = joinpath(data_dir, "Newton_NCAR_5deg_init-zero.nc")

ds = Dataset(nc_filename, "r")

lon = cvt(ds["lon"][:])
lat = cvt(ds["lat"][:])
data_h = cvt(ds["h_bestfit"][:])
data_Q = cvt(ds["Q_bestfit"][:])
#data_Td = cvt(ds["Td_mean"][:]) .- 273.15

close(ds)


