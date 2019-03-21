using NCDatasets

src_data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data", "NCAR_LENS")

F_filename = joinpath(src_data_path, "LENS_B1850C5CN_gx3v7_SHF.nc")
SST_filename = joinpath(src_data_path, "LENS_B1850C5CN_gx3v7_SST.nc")

Dataset(F_filename, "r") do ds
    global nlon = ds.dim["Nx"]
    global nlat = ds.dim["Ny"]
end



