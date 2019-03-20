using NCDatasets

src_data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data", "NCAR_LENS")

F_filename = joinpath(src_data_path, "b.e11.B1850C5CN.f09_g16.005.pop.h.SHF.100001-109912.nc")
SST_filename = joinpath(src_data_path, "b.e11.B1850C5CN.f09_g16.005.pop.h.SST.100001-109912.nc")


src_data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data", "NCAR_CESM1-WACCM")
F_filename = joinpath(src_data_path, "SMART_5deg_hfds_Omon_CESM1-WACCM_piControl_r1i1p1_009601-029512.nc")
SST_filename = joinpath(src_data_path, "SMART_5deg_tos_Omon_CESM1-WACCM_piControl_r1i1p1_009601-029512.nc")


Dataset(F_filename, "r") do ds
    global nlat = ds.dim["lat"]
    global nlon = ds.dim["lon"]
end



