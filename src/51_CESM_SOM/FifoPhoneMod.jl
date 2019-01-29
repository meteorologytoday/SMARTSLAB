
module FifoPhoneMod

export FifoPhone, hello, recv, send

mutable struct FifoPhone

    recv_fn :: String
    send_fn :: String

    function FifoPhone(;
        recv :: String,
        send :: String,
    )
        return new(recv, send)
    end
end

function hello(fp::FifoPhone)
    while !isfifo(fp.recv_fn) || !isfifo(fp.send_fn)
        println("Cannot find two fifos, sleep and do it again.")
        sleep(1)
    end

    send(fp, "<<TEST>>")
    recv_msg = recv(fp) 
    if recv_msg != "<<TEST>>"
        throw(ErrorException("Weird message: " * recv_msg))
    end
end


function recv(fp::FifoPhone)
    local result

    open(fp.recv_fn, "r") do io
        result = strip(read(io, String))
    end

    return result
end

function send(fp::FifoPhone, msg::String)

    open(fp.send_fn, "w") do io
        write(io, msg)
    end
end

end
