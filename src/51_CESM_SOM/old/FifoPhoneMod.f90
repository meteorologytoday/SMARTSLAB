module FifoPhoneMod
implicit none

type FifoPhone
    Integer :: recv_fd
    Integer :: send_fd
    character(len = 256) :: recv_fn
    character(len = 256) :: send_fn
end type FifoPhone


contains

subroutine recv(fp, msg)
    implicit none
    type(FifoPhone)  :: fp
    character(len=*) :: msg

    integer :: io


    open(fp%recv_fd, file=fp%recv_fn, form="formatted", access="stream", action="read", iostat=io)

    io = 0
    read (fp%recv_fd, '(A)', iostat=io) msg

    if (io == 0 .or. io == -1) then
        print *, "Successfully read."
    else
        print *, "Weird io state: ", io

    end if

    close(fp%recv_fd)

end subroutine



subroutine send(fp, msg)
    implicit none
    type(FifoPhone)  :: fp
    character(len=*) :: msg

    integer :: io


    open(fp%send_fd, file=fp%send_fn, form="formatted", access="stream", action="write", iostat=io)

    io = 0
    write (fp%send_fd, *, iostat=io) msg

    if (io == 0) then
        print *, "Successfully sent."
    else
        print *, "Weird io state: ", io

    end if

    close(fp%send_fd)

end subroutine


subroutine hello(fp)
    implicit none
    type(FifoPhone) :: fp
    character(256) :: msg

    call recv(fp, msg)
    msg = trim(msg)

    if ((msg .eq. "<<TEST>>") .or. (len(msg) .eq. len("<<TEST>>"))) then
        print *, "Recv hello!"
    else
        print *, len(msg), " : ", len("<<TEST>>")
        print *, "Weird msg: [", msg, "]"
    end if

    call send(fp, "<<TEST>>")

end subroutine



end module FifoPhoneMod
