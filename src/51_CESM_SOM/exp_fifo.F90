program exp_FIFO
implicit none

integer :: fd1 = 10, io
logical :: file_exists
character(len=1024) :: input, fifo="phone.fifo"

    file_exists = .false.

    do while (file_exists .eqv. .false.)
        print *, "Try to check if file is there..."
        inquire(FILE=fifo, EXIST=file_exists)
        call Sleep(2)
    end do

    io = 1

    do while (io /= 0)
        print *, "Let's open it... "
        io = 0
        !open(fd1, file=fifo, form="formatted" , action="write", iostat=io)
        open(fd1, file=fifo, form="formatted", access="stream" , action="read", iostat=io)
        print *, io
        call Sleep(1)
    end do

    do
        print *, "Try to read" 
        io = 0
        input = ""
        read (fd1, *, iostat=io) input

        if (io == 0) then
            print *, "Got input: [", input, "]."
        else
            print *, "Weird io state: ", io

        end if
        
        call Sleep(2)
    end do


    close(fd1)






end program
