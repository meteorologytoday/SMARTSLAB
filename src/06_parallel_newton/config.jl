model_name = "NCAR_5deg"

include("../01_config/general_config.jl")
using Printf
using Formatting

# Check:
#
# model_name
# nchains, num_samples, num_warmup
# exp_name [contains initial condition status]
# 
# sbatch --job-name, --output, --cpus-per-task, --time, --array

function nanmean(x; dims)
    n = isnan.(x)

    y = copy(x)
    y[n] .= 0.0

    return sum(y; dims=dims) ./ sum(1 .- n; dims=dims)
end


# Newton setting
newton_fail_max = 20
newton_η = 1e-2
fail_count_max = 5

# Posterior control
σ_ϵ          = 10.0
σ_Q          = 100.0
σ_h          = 1.0
h_rng        = [0, 5000.0]
verbose      = true

param_path = collect(range(0.0, stop=1.0, length=11))

exp_name = format("Newton_{}_init-SLAB", model_name)


@printf("newton_fail_max : %d\n",   newton_fail_max)
@printf("newton_η        : %.2e\n", newton_η)
@printf("fail_count_max  : %d\n",   fail_count_max)

default_log_file = joinpath(root_path, "logs", format("{}.log", exp_name))

function writelog(lon_i, fmt, params...; ending="\n", log_file=default_log_file)
    open(log_file, "a") do f
        write(f, format("[lon_i = {}] ", lon_i))
        write(f, format(fmt, params...))
        write(f, ending)
    end
end
