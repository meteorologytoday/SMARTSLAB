include("FifoPhoneMod.jl")

using Formatting
using Printf
using .FifoPhoneMod

recv_fifo = "cesm2mymodel.fifo"
send_fifo = "mymodel2cesm.fifo"


fp = FifoPhone(recv=recv_fifo, send=send_fifo)


hello(fp)

# Mimic

while true

    println("Try to recv new msg")
    msg = recv(fp)
    println("Msg recv: ", msg)

    if msg == "<<END>>"
        println("Simulation ends!")
        break
    end



    println("Now I am doing some magical SMARTSLAB computation.")
    sst_fn = format("SST_{:03d}.nc", convert(Integer, floor(rand() * 1000)))

    sleep(2)

    println("Gonna send SST file back : ", sst_fn)

    send(fp, sst_fn)

end



