
include("../01_config/paths.jl")
using NCDatasets
# Read PDO mode


PDO_fn = joinpath(data_path, "PDO_EOFs_5deg.nc")
ds = Dataset(PDO_fn,"r")
PDO_mode = nomissing(ds["EOFs"][:, :, 1], NaN)
close(ds)

println(size(PDO_mode))

# Read SST
# Inner product
# Output
