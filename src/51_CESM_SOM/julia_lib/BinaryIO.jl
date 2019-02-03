
module BinaryIO

function readBinaryField!(
    filename::AbstractString,
    arr::Array{Float64,1},
    buffer::Array{UInt8,1};
    endianess::Symbol=:little_endian
)
    local actual_nbread


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

    if endianess == :little_endian
        for i = 1:length(arr)
            local val
            ii = (i-1) * 8
            for j = 1:8

                val = (
                    (UInt64(buffer[ii + 1]) << 0)  +
                    (UInt64(buffer[ii + 2]) << 8)  +
                    (UInt64(buffer[ii + 3]) << 16) +
                    (UInt64(buffer[ii + 4]) << 24) +
                    (UInt64(buffer[ii + 5]) << 32) +
                    (UInt64(buffer[ii + 6]) << 40) +
                    (UInt64(buffer[ii + 7]) << 48) +
                    (UInt64(buffer[ii + 8]) << 56)
                )
            end
            arr[i] = reinterpret(Float64, val)
        end

    elseif endianess == :big_endian
        for i = 1:length(arr)
            local val
            ii = (i-1) * 8
            for j = 1:8

                val = (
                    (UInt64(buffer[ii + 1]) << 56) +
                    (UInt64(buffer[ii + 2]) << 48) +
                    (UInt64(buffer[ii + 3]) << 40) +
                    (UInt64(buffer[ii + 4]) << 32) +
                    (UInt64(buffer[ii + 5]) << 24) +
                    (UInt64(buffer[ii + 6]) << 16) +
                    (UInt64(buffer[ii + 7]) << 8) +
                    (UInt64(buffer[ii + 8]) << 0)
                )
            end
            arr[i] = reinterpret(Float64, val)
        end
    else

        throw(ErrorException("Unknown endianess: " * String(endianess)))

    end
end

end
