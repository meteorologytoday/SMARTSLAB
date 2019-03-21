include("SSM_config.jl")


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

println("===== INITIALIZING SSM =====")


stage = :INIT
mail = MailboxInfo()
map = NetCDFIO.MapInfo{Float64}(domain_file)

include("init_ocean.jl")

ncio = NetCDFIO.MapInfo{Float64}(domain_file)
buffer2d  = zeros(UInt8, map.lsize * 8)
sst       = zeros(Float64, map.lsize)
mld       = copy(sst)
hflx      = copy(sst)
swflx     = copy(sst)
taux      = copy(sst)
tauy      = copy(sst)
qflux2atm = copy(sst)

sumflx      = copy(sst)

output_vars = Dict(
    "mld"       => reshape(mld,       map.nx, map.ny),
    "sst"       => reshape(sst,       map.nx, map.ny),
    "sumflx"    => reshape(sumflx,    map.nx, map.ny),
    "qflux2atm" => reshape(qflux2atm, map.nx, map.ny),
)


# Mask data
SSM.maskData!(occ, sst)
SSM.maskData!(occ, mld)
SSM.maskData!(occ, qflux2atm)

time_i = 1 
wrap_time = i -> ((time_i-1) % output_record_length) + 1
output_filename = ""
println("===== SSM IS READY =====")

beg_time = Base.time()
while true

    global stage, time_i, output_filename

    end_time = Base.time()

    println(format("Execution time: {:d}", floor(end_time - beg_time)))
    if wrap_time(time_i) == 1
        output_filename = format("SSM_output_{:03d}.nc", convert(Integer, 1+floor((time_i-1) / output_record_length)))
        
        NetCDFIO.createNCFile(map, output_filename)
    end

    println(format("# Time counter : {:d}", time_i))
    println(format("# Stage        : {}", String(stage)))

    msg = parseMsg(recv(mail))
    println("==== MESSAGE RECEIVED ====")
    print(json(msg, 4))
    println("==========================")

    if stage == :INIT && msg["MSG"] == "INIT"

        SSM.getInfo!(occ=occ, sst=sst, mld=mld)
        writeBinary!(msg["SST"], sst, buffer2d; endianess=:little_endian)
        send(mail, msg["SST"])

        NetCDFIO.write2NCFile(
            map, output_filename, output_vars;
            time=wrap_time(time_i),
            missing_value=map.missing_value
        )

        time_i += 1

        stage = :RUN
        
    elseif stage == :RUN && msg["MSG"] == "RUN"
        readBinary!(msg["HFLX"],   hflx, buffer2d; endianess=:little_endian, delete=false)
        readBinary!(msg["SWFLX"], swflx, buffer2d; endianess=:little_endian, delete=false)
        readBinary!(msg["TAUX"],   taux, buffer2d; endianess=:little_endian, delete=false)
        readBinary!(msg["TAUY"],   tauy, buffer2d; endianess=:little_endian, delete=false)

        hflx  .*= -1.0
        swflx .*= -1.0
        sumflx .= hflx + swflx
        SSM.maskData!(occ, sumflx)
        
        println("Do complicated, magical calculations...")
        
        SSM.stepOceanColumnCollection!(
            occ   = occ,
            u     = taux,
            v     = tauy,
            hflx  = hflx,
            swflx = swflx,
            Δt    = parse(Float64, msg["DT"])
        )

 
        SSM.getInfo!(occ=occ, sst=sst, mld=mld)

        writeBinary!(msg["SST_NEW"], sst, buffer2d; endianess=:little_endian)
        send(mail, msg["SST_NEW"])

        NetCDFIO.write2NCFile(map, output_filename, output_vars; time=wrap_time(time_i), missing_value=map.missing_value)
        
        time_i += 1

    elseif stage == :RUN && msg["MSG"] == "END"

        println("Simulation ends peacefully.")
        send(mail, "END")
        break
    else
        send(mail, "CRASH")
        throw(ErrorException("Unknown status: stage " * stage * ", MSG: " * String(msg["MSG"])))
    end

end


