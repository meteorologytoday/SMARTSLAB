using Formatting

# Check:
#
# model_name
# nchains, num_samples, num_warmup
# exp_name [contains initial condition status]
# 
# sbatch --job-name, --output, --cpus-per-task, --time, --array



model_name = "NCAR_5deg"
include("../../01_config/general_config.jl")

nchains     = 4
num_samples = 1000
num_warmup  = 200

θd = 273.15 * ρ * c_p
σ_θ = 1.0 * ρ * c_p
σ_Q = 100.0


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
