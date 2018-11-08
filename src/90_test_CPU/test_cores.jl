

jobs_size = 1000
channel_size = 100
mtx_size = 1000

jobs = Channel(channel_size)
results = Channel(channel_size)

function do_work()
    for job_id in jobs

        beg_t = Base.time()
        m = zeros(mtx_size, mtx_size)
        m * m
        end_t = Base.time()
        put!(results, (job_id, end_t - beg_t))
    end
end

function make_jobs(n)
    for i = 1:n
        put!(jobs, i)
    end
end

@async make_jobs(jobs_size)
for i in 1:jobs_size
    @async do_work()
end

jobs_left = jobs_size

@elapsed while jobs_left > 0
    job_id, exec_time = take!(results)
    println("$job_id finished in $(round(exec_time; digits=2)) seconds")
    global jobs_left -= 1
end



