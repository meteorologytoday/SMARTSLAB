using Formatting

# Check:
#
# model_name
# nchains, num_samples, num_warmup
# exp_name [contains initial condition status]
# 
# sbatch --job-name, --output, --cpus-per-task, --time, --array

include("setup_src_data.jl")

data_path = normpath(joinpath(dirname(@__FILE__)), "..", "..", "..", "data")

nchains     = 2
num_samples = 1000
num_warmup  = 50


ρ    = 1027.0  # kg / m^3
c_p  = 3985.0  # J / kg / K
θd   = 273.15 * ρ * c_p
σ_θ  = 1.0 * ρ * c_p
σ_Q  = 100.0

max_pts_per_task = 100.0


model_name = "NCAR_LENS_g37"
exp_name = format("stanfit_MLM2L_weak_{}_c{:d}_s{:d}_w{:d}", model_name, nchains, num_samples, num_warmup)

main_dir = joinpath(data_path, exp_name)
mkpath(main_dir)

default_log_file = normpath(joinpath(dirname(@__FILE__), format("progress.log", exp_name)))
function writelog(lon_i, fmt, params...; ending="\n", log_file=default_log_file)
    open(log_file, "a") do f
        write(f, format("[lon_i = {}] ", lon_i))
        write(f, format(fmt, params...))
        write(f, ending)
    end
end
