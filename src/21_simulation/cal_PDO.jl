include("../lib/AnalyizeTimeseries.jl")
include("../lib/nanops.jl")
using NCDatasets
using Statistics: mean, std
using .AnalyzeTimeseries

model_name = "NCAR_5deg"
include("../01_config/general_config.jl")
include("config.jl")


# Read PDO mode

PDO_fn = joinpath(data_path, "PDO_EOFs_5deg.nc")
ds = Dataset(PDO_fn,"r")
PDO_mode = nomissing(ds["EOFs"][:, :, 1], NaN)
close(ds)

function calPDOIndex(SST)
    index = zeros(dtype, size(SST)[3])
    for i = 1:length(index)
        index[i] = nansum(SST[:, :, i] .* PDO_mode)
    end
    return (index .- mean(index)) / std(index)
end


obs_SST = readModelVar("tos", (:, :, init_time:init_time+sim_len-1))
ds = Dataset(sim_nc_filename, "r")
sim_SST = nomissing(ds["SST"][:], NaN)
close(ds)

obs_PDO = calPDOIndex(obs_SST)
sim_PDO = calPDOIndex(sim_SST)

println(size(obs_PDO))
println(size(sim_PDO))

using PyPlot

plt[:plot](obs_PDO)
plt[:plot](sim_PDO)
