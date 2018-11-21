using Formatting

model_name = "NCAR_2deg"

nchains     = 4
num_samples = 1000
num_warmup  = 200

#nchains     = 2
#num_samples = 10
#num_warmup  = 2


exp_name = format("HMC_{}_c{:d}_s{:d}_w{:d}", model_name, nchains, num_samples, num_warmup)


default_log_file = joinpath("/export/home/tienyiah/projects/SMARTSLAB", format("{}.log", exp_name))
function writelog(job_id, fmt, params...; ending="\n", log_file=default_log_file)
    open(log_file, "a") do f
        write(f, format("[{}] ", job_id))
        write(f, format(fmt, params...))
        write(f, ending)
    end
end
