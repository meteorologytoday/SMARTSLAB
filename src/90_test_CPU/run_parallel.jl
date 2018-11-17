include("../01_config/general_config.jl")

ch = 100
max_workers = 4

script_dir = dirname(@__FILE__)
script_file = joinpath(script_dir, "test_cores.jl")

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
            Base.run(`julia $script_file`)
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

n = 0
@elapsed while true
    result = take!(results)
    global n += 1
    println(format("Progress: {:.2f} % ( {:d} / {:d} )", 100.0 * n / length(lon), n, length(lon)))
    if n == length(lon)
        println("All the work is done.")
        break
    end
end

