program_beg_time = Base.time()

include("../01_config/general_config.jl")
include("../01_config/regions.jl")
include("status_code.jl")

using Formatting
using NCDatasets

#=
θ_data      = readModelVar("tos") * ρ * c_p
F_data      = readModelVar("hfds")
omlmax_data = readModelVar("omlmax")
B_data = zeros(dtype, length(time))
=#

status_attrib = Dict(
    "_FillValue" => missing_value,
    "_NAN"       => STATUS[:NAN],
    "_UNDONE"    => STATUS[:UNDONE],
    "_DOING"     => STATUS[:DOING],
    "_DONE"      => STATUS[:DONE],
    "_ERROR"      => STATUS[:ERROR],
)

verbose = false
parallel_n = 6

# Create a file for checking status
status_filename = "status.nc"
status_file_lock = false
need_create = true
if isfile(status_filename)
    println("Status file ", status_filename, " exists.")
    ds = Dataset(status_filename, "r")
    if haskey(ds, "status")
        need_create = false
    else
        need_create = true
        println("File exists but no variable [status].")
    end
end

if need_create
    println("Status file need to be (re)created: ", status_filename)
    mask = readModelVar("tos", (:, :, 1))
    mask[isfinite.(mask)] .= STATUS[:UNDONE]
    mask[isnan.(mask)]    .= STATUS[:NAN]

    ds = Dataset(status_filename, "c")
    defDim(ds, "lat", length(lat))
    defDim(ds, "lon", length(lon))

    defVar(ds, "lat",  dtype, ("lat",))[:]  = lat
    defVar(ds, "lon",  dtype, ("lon",))[:]  = lon

    v = defVar(ds, "status", dtype, ("lon", "lat"))
    for key in keys(status_attrib)
        v.attrib[key] = status_attrib[key]
    end
    v[:] = mask

    close(ds)
end


ds = Dataset(status_filename,"r")
status = convert(Array{dtype}, nomissing(ds["status"][:], NaN))
close(ds)

status[(status .== STATUS[:DOING]) .| (status .== STATUS[:ERROR])] .= STATUS[:UNDONE]

done_count, undone_count, total_count = NaN, NaN, NaN
function reportStatus()
    global done_count = sum( status .== STATUS[:DONE])
    global undone_count = sum( status .== STATUS[:UNDONE])
    global total_count = sum( status .!= STATUS[:NAN])

    println(
        format("Progress: {:.1f} % ({:d} / {:d})",
            done_count / total_count * 100.0,
            done_count,
            total_count
        )
    )
end

function markStatus(i, j, stat)
    status[i, j] = STATUS[stat]
    
    if sum(status[i, :] .== STATUS[:UNDONE]) == 0
        try
            println("lon[", i, "] = ", lon[i], " is all complete and now update status.")
            updateStatus()
        catch e
            throw(e)
        end
    end
end

function updateStatus()
    verbose && println("backup")
    backupStatus()
    verbose && println("write")
    writeStatus()
    verbose && println("report")
    reportStatus()
end


function backupStatus()
    try
        cp(status_filename, "tmp.nc", force=true)
    catch e
        println("Copying file failed. Keep going...")
        println(e)
    end
end

function writeStatus()

    try
        verbose && println("create")
        ds = Dataset(status_filename, "c")
        defDim(ds, "lat", length(lat))
        defDim(ds, "lon", length(lon))

        defVar(ds, "lat",  dtype, ("lat",))[:]  = lat
        defVar(ds, "lon",  dtype, ("lon",))[:]  = lon
        
        v = defVar(ds, "status", dtype, ("lon", "lat"))
        verbose && println("defined variable")
        for key in keys(status_attrib)
            v.attrib[key] = status_attrib[key]
        end
        verbose && println("defined attrib")
        v[:] = status
        verbose && println("output variable")


    catch e
        throw(e)
    finally
        close(ds)        
    end

end


reportStatus()
if undone_count == 0
    println("Nothing to do. All jobs are complete. Now exiting the program...")
    exit()
end


jobs    = Channel(2)
results = Channel(2)

job_assigned = 0
function make_jobs()
    println("Start making jobs")
    for i = 1:length(lon), j = 1:length(lat)
        try
        if isnan(status[i,j]) || status[i,j] != STATUS[:UNDONE]
            continue
        end

        #println("i, j = ", i, ", ", j)
        status[i, j] = STATUS[:DOING]

        # load data
        #=
        θ      = θ_data[i, j, :]
        S      = F_data[i, j, :]
        B      = B_data
        omlmax = omlmax_data[i, j, :]
        =#
        θ      = 1
        S      = 2
        B      = 3
        omlmax = 4
        
        put!(jobs,
            Dict(
                "lonlat_i" => (i, j),
                "θ"        => θ,
                "S"        => S,
                "B"        => B,
                "omlmax"   => omlmax
            )
        )

        catch e
            throw(e)
        end
    end
    close(jobs)
    println("finished making jobs.")
end


function do_jobs(worker_id)
    println("[", worker_id, "] Ready to do jobs.")
    for job in jobs
        i, j = job["lonlat_i"]

        markStatus(i, j, :DONE)
        put!(results, "done")
    end
    println("[", worker_id, "] Work done.")
end

println("Call make_jobs()")
@async make_jobs()



for i in 1:parallel_n
    println("Call do_jobs() with i = ", i)
    @async do_jobs(i)
end

n = undone_count
@elapsed for result in results
    global n-= 1
    if n == 0
        break
    end
end

