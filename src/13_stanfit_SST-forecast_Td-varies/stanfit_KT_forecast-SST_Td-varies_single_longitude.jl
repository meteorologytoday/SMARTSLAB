program_beg_time = Base.time()

include("config.jl")
using Printf
using Formatting
using NCDatasets
import Statistics: mean, std

if length(ARGS) != 1 
    throw(ErrorException("Length of ARGS must be 1. That is the longitude index."))
end

lon_i = parse(Int, ARGS[1])

FILE_noext = splitext(basename(@__FILE__))[1]

# construct data folder
main_dir = joinpath(data_path, FILE_noext, exp_name)
tmp_dir = joinpath(main_dir, "stan_tmp", format("{:03d}", lon_i))
mkpath(main_dir)
mkpath(tmp_dir)

filename = format("{:03d}.jld", lon_i)
filename = joinpath(main_dir, filename)

println("This program is going to fit lon[", lon_i, "] = ", lon[lon_i])
if isfile(filename)
    println("File ", filename, " already exists. End program.")
    exit()
end

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan
@printf("done\n")

script_path = normpath(joinpath(dirname(@__FILE__), "..", "STAN_code", "KT_forecast-SST_Td-varies.stan"))
model_script = read(script_path, String)
stanmodel = Stanmodel(
    name=FILE_noext,
    nchains=nchains,
    num_samples=num_samples,
    num_warmup=num_warmup,
    model=model_script,
    pdir=tmp_dir,
)

#display(stanmodel)

h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q.{}", i) for i = 1:12 ]

data = Dict()
time_stat = Dict()

total_time = 0.0

β_mean = zeros(dtype, length(lat), 25)
β_std = zeros(dtype, length(lat), 25)

β_mean .= NaN
β_std  .= NaN

crop_range = (lon_i, :, 1:1200) # first 100 years

omlmax = readModelVar("omlmax", crop_range)
θ      = readModelVar("tos", crop_range) * (ρ * c_p)
F      = readModelVar("hfds", crop_range)

N = size(θ)[2]

for j = 1:length(lat)

    if isnan(θ[j, 1]) # skip land
        continue
    end

    writelog(lon_i, "[lat_i = {:03d}] Doing lat[{:d}] = {:.2f}", j, lon_i, j, lat[j] )

    beg_time = Base.time()

    omlmax_mean = mean( reshape(omlmax[j, :], 12, :); dims=(2,) )[:] # need squeeze dimension
    println("size of init_omlmax: ", size(omlmax_mean))
    println("init_omlmax: ", omlmax_mean)

    data = Dict(
        "raw_N"       => N, 
        "period"      => 12, 
        "dt"          => Δt, 
        "theta_d"     => θd, 
        "raw_theta"   => θ[j, :], 
        "raw_F"       => F[j, :],
        "theta_std"   => σ_θ,
        "Q_std"       => σ_Q,
    )

    init = Dict(
        "h" => omlmax_mean * 0.0 .+ 30.0
    )
    
    rc, sim1 = stan(
        stanmodel,
        [data];
        init = [init],
        CmdStanDir=ENV["CMDSTAN_HOME"]
    )

    if rc != 0
        println("There are errors!!")
        continue
    end
    
    println("Extracting result...")
    h_mean = zeros(12)
    h_std  = zeros(12)
    Q_mean = zeros(12)
    Q_std  = zeros(12)

    data_h = sim1[:, h_key, :].value
    data_Q = sim1[:, Q_key, :].value
    data_Td = sim1[:, ["theta_d"], :].value / ρ / c_p

    for i = 1:12
        h_mean[i] = mean(data_h[:, i, :])
        h_std[i]  = std(data_h[:, i, :])

        Q_mean[i] = mean(data_Q[:, i, :])
        Q_std[i]  = std(data_Q[:, i, :])
    end

    Td_mean = mean(data_Td)
    Td_std  = std(data_Td)

    println("h_mean", h_mean)
    println("Q_mean", Q_mean)

    β_mean[j,  1:12] = h_mean
    β_mean[j, 13:24] = Q_mean
    β_mean[j,    25] = Td_mean

    β_std[j,  1:12] = h_std
    β_std[j, 13:24] = Q_std
    β_std[j,    25] = Td_std

    time_stat = Base.time() - beg_time

    global total_time += time_stat
    writelog(lon_i, "[lat_i = {:03d}] Stan fit: {:.2f} min. Total time: {:.2f} min. ", j, time_stat / 60.0, total_time / 60.0 )

end

using JLD

println("Output filename: ", filename)
save(filename, Dict("β_mean" => β_mean, "β_std" => β_std))


program_end_time = Base.time()

@printf("Total time used: %.2f min for %d points, with nchains = %d, num_samples = %d\n",
    (program_end_time-program_beg_time)/ 60.0,
    length(lat),
    nchains,
    num_samples
)
