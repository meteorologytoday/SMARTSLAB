using Formatting

# Check:
#
# model_name
# nchains, num_samples, num_warmup
# exp_name [contains initial condition status]
# 
# sbatch --job-name, --output, --cpus-per-task, --time, --array



model_name = "NCAR_5deg"

nchains     = 4
num_samples = 100
num_warmup  = 20

#nchains     = 2
#num_samples = 10
#num_warmup  = 2


exp_name = format("HMC_{}_init-omlmax_c{:d}_s{:d}_w{:d}", model_name, nchains, num_samples, num_warmup)


default_log_file = joinpath("/export/home/tienyiah/projects/SMARTSLAB", format("{}.log", exp_name))
function writelog(lon_i, fmt, params...; ending="\n", log_file=default_log_file)
    open(log_file, "a") do f
        write(f, format("[lon_i = {}] ", lon_i))
        write(f, format(fmt, params...))
        write(f, ending)
    end
end
