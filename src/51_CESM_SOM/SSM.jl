include("julia_lib/Mailbox.jl")
include("julia_lib/BinaryIO.jl")

using Formatting
using Printf
using .Mailbox
using .BinaryIO

wdir_file = "./working_directory"
if isfile(wdir_file)
    wdir =  strip(read(wdir_file, String))
    if isdir(wdir)
        cd(wdir)
    else
        throw(ErrorException("Working directory [ " * wdir * " ] does not exist."))
    end
end

function parseMsg(msg::AbstractString)
    pairs = split(msg, ";")
    d = Dict{AbstractString, Any}()
    for i = 1:length(pairs)

        if strip(pairs[i]) == ""
            continue
        end

        key, val = split(pairs[i], ":")
        d[key] = val
    end
    return d
end

function writeBinaryField(filename::AbstractString, arr::Array{Float64})
    n_bytes = write(filename, arr)
    if n_bytes != length(arr) * 8
        throw(ErrorException(format("Only {:d} bytes were written, while there are {:d} bytes needed to be written", n_bytes, length(arr)*8)))
    end
end




nlon = 24
nlat = 19
lsize = nlon * nlat
stage = :INIT
MI = MailboxInfo()


buffer2d  = BitArray(undef, lsize * 64)
sst       = zeros(Float64, lsize)
hflx      = zeros(Float64, lsize)
swflx     = zeros(Float64, lsize)

while true

    global stage

    println("Try to recv new msg")
    msg = parseMsg(recv(MI))
    println("Msg recv: ", msg)

    if stage == :INIT && msg["MSG"] == "INIT"
        readBinaryField!(msg["SST"], sst, buffer2d)

        send(MI, "READY")

        stage = :RUN

    elseif stage == :RUN && msg["MSG"] == "RUN"
        readBinaryField!(msg["HFLX"], hflx, buffer2d)
        readBinaryField!(msg["SWFLX"], swflx, buffer2d)

        println("Do complicated, magical calculations...")

        writeBinaryField(msg["SST_NEW"], sst)

        send(MI, msg["SST_NEW"])

    elseif stage == :RUN && msg["MSG"] == "END"

        println("Simulation ends peacefully.")
        send(MI, "END")
        break
    else
        
        send(MI, "CRASH")
        throw(ErrorException("Unknown status: stage " * stage * ", MSG: " * msg["MSG"]))
    end

end



