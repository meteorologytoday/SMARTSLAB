using NCDatasets
using Formatting

cvt = x -> convert(Array{Float64}, nomissing(x, NaN))
stat = x -> (mean(x, dims=(1,))[1,:, :], std(x, dims=(1,))[1,:,:])



# Reading data
data_dir = normpath(joinpath( dirname(@__FILE__), "..", "..", "data"))
img_dir  = normpath(joinpath( dirname(@__FILE__), "..", "..", "img"))
nc_filename = joinpath(data_dir, "HMC_NCAR_5deg_init-omlmax_c4_s1000_w50.nc")
#nc_filename = joinpath(data_dir, "HMC_NCAR_2deg_c4_s1000_w200.nc")

ds = Dataset(nc_filename, "r")

lon = cvt(ds["lon"][:])
lat = cvt(ds["lat"][:])
data_h = cvt(ds["h_mean"][:])
data_Q = cvt(ds["Q_mean"][:])
data_Td = cvt(ds["Td_mean"][:]) .- 273.15

close(ds)


