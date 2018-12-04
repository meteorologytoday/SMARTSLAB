program_beg_time = Base.time()

include("config.jl")
include("../01_config/general_config.jl")

using Printf
using Formatting
using NCDatasets
import Statistics: mean, std

if length(ARGS) != 1
    throw(ErrorException("Length of ARGS must be 1. That is the lon index."))
end

lon_i = parse(Int, ARGS[1])

# construct data folder
main_dir = joinpath(data_path, splitext(basename(@__FILE__))[1], exp_name)
tmp_dir = joinpath(main_dir, "stan_tmp", format("{:03d}", lon_i))
mkpath(main_dir)
mkpath(tmp_dir)

filename = format("{:03d}_{:03d}.jld", lon_i, lat_i)
filename = joinpath(main_dir, filename)

println(format(
    "This program is going to fit lon[{:d}] = {:.2f}, lat[{:d}] = {:.2f}", lon_i, lon[lon_i], lat_i, lat[lat_i]
))

if isfile(filename)
    println("File ", filename, " already exists. End program.")
    exit()
end

println("ENV[\"CMDSTAN_HOME\"] = ", ENV["CMDSTAN_HOME"])
@printf("Importing Stan library...")
using Stan
@printf("done\n")

model_script = read(joinpath(dirname(@__FILE__), "..", "lib", "STAN_code", "KT.stan"), String)
stanmodel = Stanmodel(
    name="KT",
    nchains=nchains,
    num_samples=num_samples,
    num_warmup=num_warmup,
    model=model_script,
    pdir=tmp_dir,
)

#display(stanmodel)

h_key = [ format("h.{}", i) for i = 1:12 ]
Q_key = [ format("Q_s.{}", i) for i = 1:12 ]

data = Dict()
time_stat = Dict()


β = zeros(dtype, num_samples, 25, nchains)
β .= NaN

crop_range = (lon_i, lat_i, :)


omlmax = readModelVar("omlmax", crop_range)
θ      = readModelVar("tos", crop_range) * (ρ * c_p)
F      = readModelVar("hfds", crop_range)

N = length(θ)

if isnan(θ[1]) # skip land
    throw(ErrorException("Invalid data on this grid point."))
end

omlmax_mean = mean( reshape(omlmax, 12, :); dims=(2,) )[:] # need squeeze dimension
println("size of init_omlmax: ", size(omlmax_mean))
println("init_omlmax: ", omlmax_mean)

data = Dict(
    "N"           => N, 
    "period"      => 12, 
    "dt"          => Δt, 
    "theta"       => θ, 
    "F"           => F,
    "epsilon_std" => 10.0,
    "Q_std"       => 100.0,
)

init = Dict(
    "h" => omlmax_mean
)

rc, sim1 = stan(
    stanmodel,
    [data];
    init = [init],
    CmdStanDir=ENV["CMDSTAN_HOME"]
)

if rc != 0
    println("There are errors!!")
    throw(ErrorException("There are errors when executing STAN program!!"))
end

println("Extracting result...")
h_mean = zeros(12)
h_std  = zeros(12)
Q_mean = zeros(12)
Q_std  = zeros(12)

data_h = sim1[:, h_key, :].value
data_Q = sim1[:, Q_key, :].value
data_Td = sim1[:, ["theta_d"], :].value / ρ / c_p

β[:,  1:12, :] = data_h
β[:, 13:24, :] = data_Q
β[:,    25, :] = data_Td

# output data...

using JLD

println("Output filename: ", filename)
save(filename, Dict("β" => β))

program_end_time = Base.time()

@printf("Total time used: %.2f min for %d points, with nchains = %d, num_samples = %d\n",
    (program_end_time-program_beg_time)/ 60.0,
    length(lat),
    nchains,
    num_samples
)
