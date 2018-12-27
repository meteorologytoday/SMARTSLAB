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
    
    #=
    # Find global trend
    gSST = zeros(dtype, size(SST)[3])
    for i = 1:length(gSST)
        weight .* SST[:, :, i]
        gSST[i] = nansum(weight .* SST[:, :, i]) / sum_weight
    end
    =# 
    
    # Detrend

    x = collect(dtype, 1:ntime)
    for i=1:length(lon), j=1:length(lat)

        β = LinearRegression(x, SST[i, j, :])
        SST[i, j, :] -= β[1] .+ β[2] * x
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

obs_F   = readModelVar("hfds", (:, :, init_time:init_time+sim_time-1))
obs_SST = readModelVar("tos", (:, :, init_time:init_time+sim_time-1))

sim_nc_filename = joinpath(data_path, "21_KTsimulate_HMC_SST_Td-fixed_NCAR_5deg_init-30m_c4_s1000_w200.nc")
ds = Dataset(sim_nc_filename, "r")
sim_SST = nomissing(ds["SST"][:], NaN)
close(ds)

rmMeanStates!(obs_F  ; period=12)
rmMeanStates!(obs_SST; period=12)
rmMeanStates!(sim_SST; period=12)

obs_F_PDO = calPDOIndex(obs_F)
obs_PDO = calPDOIndex(obs_SST)
sim_PDO = calPDOIndex(sim_SST)


using JLD
jld_fn = joinpath(data_path, "simulated_and_cpld_PDO.jld")
save(jld_fn, Dict(
    "obs_F_PDO" => obs_F_PDO,
    "obs_PDO" => obs_PDO,
    "sim_PDO" => sim_PDO,
))


            

