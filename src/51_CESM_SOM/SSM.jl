include("MailboxMod.jl")

using Formatting
using Printf
using .MailboxMod




function parseMsg(msg::String)
    pairs = split(msg, ";")
    d = Dict{String, Any}()
    for i = 1:length(pairs)

        if strip(pairs[i]) == ""
            continue
        end

        key, val = split(pairs[i], ":")
        d[key] = val
    end
    return d
end

function readBinaryField!(filename::String, arr::Array{Float32}, buffer::BitArray)
    if length(arr) * 32 != length(buffer)
        throw ErrorException("Length of buffer should be exactly 32 times of data array's.")
    end

    read!(filename, buffer)
    arr[:] = convert(Array{Float32}, reinterpret(buffer))
end

function writeBinaryField(filename::String, arr::Array{Float32})
    n_bytes = write(filename, arr)
    if n_bytes != length(arr) * 4
        throw ErrorException(format("Only {:d} bytes were written, while there are {:d} bytes needed to be written", n_bytes, length(arr)*4))
    end
end




nlon = 24
nlat = 19
lsize = nlon * nlat
stage = :INIT
MI = MailboxInfo()


buffer2d  = BitArray(undef, lsize * 32)
sst       = Array(Float32, lsize)
hflx      = Array(Float32, lsize)
swflx     = Array(Float32, lsize)

while true
    println("Try to recv new msg")
    msg = parse(recv(MI))
    println("Msg recv: ", msg)

    if stage == :INIT && msg["MSG"] == "INIT"
        readBinaryField!(msg["SST_FILE"], sst, buffer2d)

        send(MI, "READY")

        stage = :RUN

    else if stage == :RUN && msg["MSG"] == "RUN"
        readBinaryField!(msg["HEAT_FLUX"], hflx, buffer2d)
        readBinaryField!(msg["SWAV_FLUX"], swflx, buffer2d)

        println("Do complicated, magical calculations...")


        sst_filename = format("sst.{:s}.bin", msg["TIME"])

        send(MI, msg["EXPECTED_SST_FILE"])

    else if stage == :RUN && msg["MSG"] == "END"

        println("Simulation ends!")
        break
    else
        
        send(MI, "CRASH")
        throw ErrorException("Unknown status: stage " * stage * ", MSG: " * msg["MSG"])
    end

end



