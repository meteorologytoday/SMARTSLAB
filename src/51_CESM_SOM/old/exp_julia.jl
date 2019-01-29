using Formatting
using Printf

fifo = "phone.fifo"

while !isfifo(fifo)

    println("Cannot find fifo, sleep and do it again.")
    sleep(2)

end


open(fifo, "w") do io
    println("file opened..")
    while true
        s = format("Helloiiii[{:f}]", rand())
        println(s)
        write(io, s)
        sleep(1)
        break
    end
end
