include("/home/tienyiah/projects/SMARTSLAB/src/22_dev_MLMML/SSM.jl")
include("julia_lib/Mailbox.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")

using Formatting
using Printf
using JSON
using .Mailbox
using .BinaryIO
using .SSM
using .NetCDFIO
using Statistics: std, mean


# 10 deg
domain_file = "/home/tienyiah/projects/cesm2_test/inputdata/share/domains/domain.ocn.fv10x15_gx3v7.180321.nc"

# 2 deg
#domain_file = "~/projects/cesm2_test/inputdata/share/domains/domain.ocn.1.9x2.5_gx1v6_090403.nc"

wdir = "/home/tienyiah/cesm/scratch/SOM_AQUAP_CICE/run"

zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5

output_record_length = 365
