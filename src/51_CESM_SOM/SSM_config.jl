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

caseroot = "/home/tienyiah/projects/cesm1_test/CICE_f45"
wdir     = "/home/tienyiah/cesm1/scratch/CICE_f45/run"

domain_file = "/home/tienyiah/cesm_inputdata/cesm1/share/domains/domain.ocn.gx3v7.120323.nc"
zs = collect(Float64, range(0, -500, step=-5))
K = 1e-5

output_record_length = 365
