using NCDatasets

src_data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data", "NCAR_LENS")

F_filename = joinpath(src_data_path, "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc")
SST_filename = joinpath(src_data_path, "b.e11.B1850C5CN.f09_g16.005.pop.h.SST.100001-109912.nc")

Dataset(F_filename, "r") do ds
    global nlat = ds.dim["nlat"]
    global nlon = ds.dim["nlon"]
end



