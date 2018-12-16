include("../lib/AnalyizeTimeseries.jl")
include("../lib/LinearRegression.jl")
include("../lib/nanops.jl")
using NCDatasets
using Statistics: mean, std
using .AnalyzeTimeseries

include("../01_config/paths.jl")

dtype = Float64
# Read PDO mode
PDO_fn = joinpath(data_path, "PDO_EOFs.nc")
ds = Dataset(PDO_fn,"r")
lat = nomissing(ds["latitude"][:], NaN)
lon = nomissing(ds["longitude"][:], NaN)
weight  = repeat(cos.(lat' * π / 180.0);outer=(length(lon), 1)) 
sum_weight = nansum(weight)

PDO_mode = weight .* nomissing(ds["EOFs"][:, :, 1], NaN)
close(ds)

# Read Reanalysis SST data
SST_fn = joinpath(data_path, "SST_1870-2017_0-360.nc")
ds = Dataset(SST_fn,"r")
SST = nomissing(ds["SST"][:, :, :], NaN)
close(ds)


function rmMeanStates!(SST; period=12)

    ntime = size(SST)[3]
    yrs  = Int(ntime / period)

    if(mod(ntime, period) != 0) 
        throw(ErrorException("Length of time is not multiple of period: ", period))
    end

    # detrend and remove annual cycle

    # Find global trend
    gSST = zeros(dtype, size(SST)[3])
    for i = 1:length(gSST)
        weight .* SST[:, :, i]
        gSST[i] = nansum(weight .* SST[:, :, i]) / sum_weight
    end 

    # Detrend
    x = collect(dtype, 1:length(gSST))
    β = LinearRegression(x, gSST)
    trend = β[1] .+ β[2]*x
    for i=1:length(lon), j=1:length(lat)
        SST[i, j, :] -= trend
    end
    
    # Remove seasonal cycle
    for i=1:length(lon), j=1:length(lat), t=1:period
        SST[i, j, t:period:ntime] .-= mean(SST[i, j, t:period:ntime])
    end
     
end


function calPDOIndex(SST)
    index = zeros(dtype, size(SST)[3])
    for i = 1:length(index)
        index[i] = nansum(SST[:, :, i] .* PDO_mode)
    end
    return (index .- mean(index)) / std(index)
end


rmMeanStates!(SST; period=12)
PDO = calPDOIndex(SST)

using PyPlot
ntime = size(SST)[3]
nyrs  = Int(ntime / 12)
t = 1870.0 .+ collect(Float64, 0:ntime-1) / 12.0
plt[:plot](t, PDO)

t = 1870.0 .+ collect(Float64, 0:nyrs-1)
PDO_annual_avg = mean(reshape(PDO, 12, :), dims=1)[1,:]
plt[:plot](t, PDO_annual_avg)

t = 1870.0 .+ collect(Float64, 0:nyrs-1)
PDO_annual_avg = mean(reshape(PDO, 12, :), dims=1)[1,:]
plt[:plot](t, PDO_annual_avg)

