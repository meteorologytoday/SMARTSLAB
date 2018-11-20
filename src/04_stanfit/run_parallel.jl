include("config_a_lon.jl")
include("../01_config/general_config.jl")
using Dates
using Formatting
ch = 100
max_workers = 7


script_dir = dirname(@__FILE__)
script_file = joinpath(script_dir, "stanfit_KT_a_lon.jl")
log_file = format("{}.log", exp_name)

jobs    = Channel(ch)
results = Channel(ch)

function make_jobs()
    for i=1:length(lon)
       push!(jobs, i)
    end
    close(jobs)
end

function do_jobs()
    for i in jobs
        #Base.run(`echo "$exp_name $i"`)
        try
            Base.run(`julia $script_file $exp_name $i`)
        catch e
            throw(e)
        end
        push!(results, "$i done")
    end
end

@async make_jobs()

for i=1:max_workers
    @async do_jobs()
end


function print_log(s)
    open(log_file, "a") do f
        write(f, s)
    end
end

print_log(format("[{:.2f} hr] Begin at {}.\n", 0.0, Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")))
n = 0
beg_time = Base.time()

@elapsed while true
    result = take!(results)
    global n += 1
    consumed_time = Base.time() - beg_time
    print_log(format("[{:.2f} hr] Progress: {:.2f} % ( {:d} / {:d} )\n", consumed_time / 3600.0, 100.0 * n / length(lon), n, length(lon)))

    if n == length(lon)
        print_log("All the work is done.\n")
        break
    end
end

