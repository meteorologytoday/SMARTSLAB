program exp_FIFO

use FifoPhoneMod

implicit none

type(FifoPhone) :: fp
character(256)  :: msg
    fp%recv_fd = 10
    fp%send_fd = 11
    fp%recv_fn = "mymodel2cesm.fifo"
    fp%send_fn = "cesm2mymodel.fifo"

    call hello(fp)
    !call recv(fp, msg)

    print *, "Program ends."

end program
