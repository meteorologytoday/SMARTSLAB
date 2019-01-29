program exp_FIFO

use FifoPhoneMod

implicit none
integer :: i
type(FifoPhone) :: fp
character(1024)  :: msg
    fp%recv_fd = 10
    fp%send_fd = 11
    fp%recv_fn = "mymodel2cesm.fifo"
    fp%send_fn = "cesm2mymodel.fifo"

    call hello(fp)


    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)
        call send(fp, msg)
        
        call recv(fp, msg)
        print *, "SST file recv: ", trim(msg)

    end do


    msg = "<<END>>"
    call send(fp, msg)

    print *, "Program ends."

end program
