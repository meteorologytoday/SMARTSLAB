include("../lib/AnalyizeTimeseries.jl")
include("../lib/LinearRegression.jl")
include("../lib/nanops.jl")
using NCDatasets
using Statistics: mean, std
using .AnalyzeTimeseries

model_name = "NCAR_5deg"
include("../01_config/general_config.jl")
include("config.jl")


# Read PDO mode
weight  = repeat(cos.(lat' * π / 180.0);outer=(length(lon), 1)) 
sum_weight = nansum(weight)
PDO_fn = joinpath(data_path, "PDO_EOFs_5deg.nc")
ds = Dataset(PDO_fn,"r")
PDO_mode = weight .* nomissing(ds["EOFs"][:, :, 1], NaN)
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


obs_SST = readModelVar("tos", (:, :, init_time:init_time+sim_len-1))
ds = Dataset(sim_nc_filename, "r")
sim_SST = nomissing(ds["SST"][:], NaN)
close(ds)

rmMeanStates!(obs_SST; period=12)
rmMeanStates!(sim_SST; period=12)

obs_PDO = calPDOIndex(obs_SST)
sim_PDO = calPDOIndex(sim_SST)


using JLD
jld_fn = joinpath(data_path, "simulated_and_cpld_PDO.jld")
save(jld_fn, Dict(
    "obs_PDO" => obs_PDO,
    "sim_PDO" => sim_PDO,
))

println(size(obs_PDO))
println(size(sim_PDO))

using PyPlot

nmons  = size(obs_SST)[3]
nyrs   = Int(nmons / 12)
t_mon  = collect(Float64, 0:nmons-1) / 12.0
t_yr   = collect(Float64, 0:nyrs-1)


fig, ax = plt[:subplots](2, 1,sharex=true)

obs_PDO_annual_avg = mean(reshape(obs_PDO, 12, :), dims=1)[1,:]
sim_PDO_annual_avg = mean(reshape(sim_PDO, 12, :), dims=1)[1,:]

ax[1][:plot](t_mon, obs_PDO)
ax[1][:plot](t_yr,  obs_PDO_annual_avg)

ax[2][:plot](t_mon, sim_PDO)
ax[2][:plot](t_yr,  sim_PDO_annual_avg)

plt[:show]()
