
using Formatting

model_name = "NCAR_5deg"
include("../01_config/general_config.jl")



init_time      = 1  # initial condition read from coupled model output
warmup_time    = 12 * 100
sim_time       = 12 * 100
total_sim_time = warmup_time + sim_time

