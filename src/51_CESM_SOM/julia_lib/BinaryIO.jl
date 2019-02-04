
module BinaryIO

function readBinary!(
    filename::AbstractString,
    arr::Array{Float64,1},
    buffer::Array{UInt8,1};
    endianess::Symbol=:little_endian
)
    local func, need_swap, actual_nbread

    if length(arr) * 8 != length(buffer)
        throw(ErrorException("Length of buffer should be exactly 8 times of data array's."))
    end


    nbread = length(buffer)
    
    open(filename, "r") do io
        actual_nbread = readbytes!(io, buffer)
    end
    if nbread != actual_nbread 
        throw(ErrorException(format("Number of bytes read: {:d} does not equal to {:d}", actual_nbread, nbread)))
    end

    arr[:] = reinterpret(Float64, buffer)

    if endianess == :little_endian
        func = ltoh
    elseif endianess == :big_endian
        func = ntoh
    else
        throw(ErrorException("Unknown symbol: " * String(endianess)))
    end

    arr[:] = func.(arr)


end

function writeBinary!(
    filename::AbstractString,
    arr::Array{Float64,1},
    buffer::Array{UInt8,1};
    endianess::Symbol=:little_endian
)
    local func

    if length(arr) * 8 != length(buffer)
        throw(ErrorException("Length of buffer should be exactly 8 times of data array's."))
    end

    nbwrite = length(buffer)

    if endianess == :little_endian
        func = htol
    elseif endianess == :big_endian
        func = hton
    else
        throw(ErrorException("Unknown symbol: " * String(endianess)))
    end

    buffer[:] = reinterpret(UInt8, func.(arr))
    actual_nbwrite = write(filename, buffer)

    if nbwrite != actual_nbwrite 
        throw(ErrorException(format("Number of bytes write: {:d} does not equal to {:d}", actual_nbwrite, nbwrite)))
    end

end


end
