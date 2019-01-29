using Formatting
using Printf

fifo = "phone.fifo"

while !isfifo(fifo)

    println("Cannot find fifo, sleep and do it again.")
    sleep(2)

end


open(fifo, "r") do io
    println("file opened..")
    while true
        println("try to read...")
        println("[" * read(io, String) * "]")
        sleep(1)
    end
end
