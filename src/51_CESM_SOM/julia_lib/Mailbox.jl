
module Mailbox

export MailboxInfo, hello, recv, send

mutable struct MailboxInfo

    recv_fn :: AbstractString
    send_fn :: AbstractString

    lock_fn :: AbstractString

    function MailboxInfo(;
        recv :: AbstractString,
        send :: AbstractString,
        lock :: AbstractString,
    )
        return new(recv, send, lock)
    end

    function MailboxInfo()
        return new("cesm2mymodel.info", "mymodel2cesm.info", "lock")
    end
end


function appendPath(MI::MailboxInfo, path::AbstractString)
    MI.recv_fn = joinpath(path, MI.recv_fn)
    MI.send_fn = joinpath(path, MI.send_fn)
    MI.lock_fn = joinpath(path, MI.lock_fn)
end

function lock(fn::Function, MI::MailboxInfo)
    obtainLock(MI)
    fn()
    releaseLock(MI)
end


function obtainLock(MI::MailboxInfo)
    while isfile(MI.lock_fn)
        sleep(1)
    end

    open(MI.lock_fn, "w") do io
    end

end

function releaseLock(MI::MailboxInfo)
    rm(MI.lock_fn, force=true)
end

function recv(MI::MailboxInfo)
    local result

    while !isfile(MI.recv_fn)
        sleep(1)
    end

    lock(MI) do
        open(MI.recv_fn, "r") do io
            result = strip(read(io, String))
        end
        rm(MI.recv_fn, force=true)
        releaseLock(MI)
    end

    return result
end

function send(MI::MailboxInfo, msg::AbstractString)
    lock(MI) do
        open(MI.send_fn, "w") do io
            write(io, msg)
        end
    end
end



function hello(MI::MailboxInfo)
    send(MI, "<<TEST>>")
    recv_msg = recv(MI) 
    if recv_msg != "<<TEST>>"
        throw(ErrorException("Weird message: " * recv_msg))
    end
end


end
