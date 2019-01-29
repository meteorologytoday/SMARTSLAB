include("FifoPhoneMod.jl")

using Formatting
using Printf
using .FifoPhoneMod

recv_fifo = "cesm2mymodel.fifo"
send_fifo = "mymodel2cesm.fifo"


fp = FifoPhone(recv=recv_fifo, send=send_fifo)
