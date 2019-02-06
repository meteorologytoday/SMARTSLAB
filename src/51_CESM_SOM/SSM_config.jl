include("/home/tienyiah/projects/SMARTSLAB/src/22_dev_MLMML/SSM.jl")
include("julia_lib/Mailbox.jl")
include("julia_lib/BinaryIO.jl")
include("julia_lib/NetCDFIO.jl")
using Formatting
using Printf
using .Mailbox
using .BinaryIO
using .SSM
using .NetCDFIO

# 10 deg
domain_file = "/home/tienyiah/projects/cesm2_test/inputdata/share/domains/domain.ocn.fv10x15_gx3v7.180321.nc"

# 2 deg
#domain_file = "~/projects/cesm2_test/inputdata/share/domains/domain.ocn.1.9x2.5_gx1v6_090403.nc"

wdir = "/home/tienyiah/cesm/scratch/SOM_AQUAP/run"

nlon = 24
nlat = 19
lsize = nlon * nlat

lons = collect(Float64, range(  0, 360, length=nlon+1))[1:end-1]
lats = collect(Float64, range(-90,  90, length=nlat))

zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5

