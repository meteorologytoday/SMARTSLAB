
using Printf

total_jobs = 10
ptr = 1

undone = total_jobs
function getUndone()

    println("getUndone")
    if undone == 0
        return "finished"
    end

    job_id = ptr
    global ptr += 1
    global undone -= 1
    println(job_id)
    return job_id
end

jobs    = Channel(2)
results = Channel(2)

function make_jobs()

    while true
        println("find getUndone")
        assign = getUndone()
        println("find getUndone2")
        println(assign)
        if assign == "finished"
            println("finished")
            break
        end
        println("find assign: ", assign)
        put!(jobs, Dict("job_id" => assign))
    end
    
    close(jobs)
end


function do_jobs()
    println("do_jobs")
    for job in jobs
        println(job)
#       @printf("Now we do job_id = %d\n", job["job_id"])
        put!(results, "done")
        println("job done")
#        ptr+=1
    end
end

println("Call make_jobs()")
@async make_jobs()


for i in 1:1
    println("Call do_jobs() with i = ", i)
    @async do_jobs()
end



n = 0
@elapsed while true
    result = take!(results)
    global n += 1
    println("Get Result: ", n)
    mod(n, 10) == 0 && println("Now we have done ", n, " jobs.")
    if n == total_jobs
        break
    end
end


