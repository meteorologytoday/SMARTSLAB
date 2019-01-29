module FifoPhoneMod
implicit none

contains

type FifoPhone
   character(len = 256) :: recv_fn
   character(len = 256) :: send_fn
end type FifoPhone


subroutine send(fd, fp, msg)

    type(FifoPhone)  :: fp
    character(len=*) :: msg

    integer :: io


    open(fd, file=fp%send_fn, form="formatted", access="stream", action="write", iostat=io)

    io = 0
    write (fd1, *, iostat=io) msg

    if (io == 0) then
        print *, "Successfully sent."
    else
        print *, "Weird io state: ", io

    end if

    close(fd)

end subroutine












end module FifoPhoneMod
