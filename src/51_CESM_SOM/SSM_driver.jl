include("SSM_config.jl")

using Statistics: std, mean

function parseMsg(msg::AbstractString)
    pairs = split(msg, ";")
    d = Dict{AbstractString, Any}()
    for i = 1:length(pairs)

        if strip(pairs[i]) == ""
            continue
        end

        key, val = split(pairs[i], ":")
        key = String(key)
        val = String(val)
        d[key] = val
    end
    return d
end




if isdir(wdir)
    cd(wdir)
else
    throw(ErrorException("Working directory [ " * wdir * " ] does not exist."))
end

println("Initializing SSM")


stage = :INIT
MI = MailboxInfo()

include("init_ocean.jl")

ncio = NetCDFIO.MapInfo(domain_file)
buffer2d  = zeros(UInt8, lsize * 8)
sst       = zeros(Float64, lsize)
hflx      = zeros(Float64, lsize)
swflx     = zeros(Float64, lsize)
u = sst * 0.0
v = sst * 0.0

function reshape_2to1!(fr::Array{Float64, 2}, to::Array{Float64, 1})
    
    if size(fr) != (nlon, nlat)
        throw(ErrorException("Dimension does not match"))
    end

    if length(to) != lsize
        throw(ErrorException("Dimension does not match"))
    end

    for j=1:nlat
        to[1+(j-1)*nlon:j*nlon] = fr[:, j]
    end
end

function reshape_1to2!(fr::Array{Float64, 1}, to::Array{Float64, 2})

    if size(to) != (nlon, nlat)
        throw(ErrorException("Dimension does not match"))
    end

    if length(fr) != lsize
        throw(ErrorException("Dimension does not match"))
    end

    for j=1:nlat
        to[:, j] = fr[1+(j-1)*nlon:j*nlon]
    end

end



println("Ready to work")
while true

    global stage

    println("Try to recv new msg")
    msg = parseMsg(recv(MI))
    println("Msg recv: ", msg)

    if stage == :INIT && msg["MSG"] == "INIT"

        SSM.getSST!(occ=occ, sst=sst)
        writeBinary!(msg["SST"], sst, buffer2d; endianess=:little_endian)
        send(MI, msg["SST"])

        stage = :RUN

    elseif stage == :RUN && msg["MSG"] == "RUN"
        readBinary!(msg["HFLX"],   hflx, buffer2d; endianess=:little_endian)
        readBinary!(msg["SWFLX"], swflx, buffer2d; endianess=:little_endian)

        hflx  .*= -1.0
        swflx .*= -1.0
        
        println("Do complicated, magical calculations...")
        #println("Avg(swflx) = ", mean(swflx), "; std: ", std(swflx))
        #println("Avg(hflx) = ", mean(hflx), "; std: ", std(hflx))
        SSM.stepOceanColumnCollection!(
            occ   = occ,
            u     = u,
            v     = v,
            hflx  = hflx,
            swflx = swflx,
            Δt    = parse(Float64, msg["DT"])
        )

 
        SSM.getSST!(occ=occ, sst=sst)
        writeBinary!(msg["SST_NEW"], sst, buffer2d; endianess=:little_endian)

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



