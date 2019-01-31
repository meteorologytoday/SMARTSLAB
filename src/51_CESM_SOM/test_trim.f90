program exp_FIFO
implicit none
character(1024)  :: msg1, msg2


    msg1 = "abc"
    msg2 = "abc"
    msg2(4:5) = CHAR(0)

    print *, "[", msg1, "]"
    print *, "[", msg2, "]"
    print *, "[", LEN(trim(msg2)), "]"

end program
