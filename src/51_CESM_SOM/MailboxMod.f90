module MailboxMod
implicit none

type MailboxInfo
    Integer :: recv_cnt
    Integer :: send_cnt
    Integer :: recv_fd
    Integer :: send_fd
    Integer :: lock_fd

    character(len = 256) :: recv_fn
    character(len = 256) :: send_fn

    character(len = 256) :: lock_fn

    character(len = 256) :: log_file
end type MailboxInfo


contains

subroutine setDefault(MI)
    implicit none
    type(MailboxInfo) :: MI

    MI%recv_fn  = "mymodel2cesm.info"
    MI%send_fn  = "cesm2mymodel.info"
 
    MI%lock_fn  = "lock"
   
    MI%log_file = "log"

    MI%recv_cnt = 0
    MI%send_cnt = 0
 
    MI%recv_fd = 10
    MI%send_fd = 11
    MI%lock_fd = 12
end subroutine setDefault

subroutine appendPath(MI, path)
    implicit none
    type(MailboxInfo) :: MI
    character(len=256) :: path

    MI%recv_fn  = path // "/" // MI%recv_fn 
    MI%send_fn  = path // "/" // MI%send_fn 
    MI%lock_fn  = path // "/" // MI%lock_fn
    MI%log_file = path // "/" // MI%log_file

end subroutine appendPath

subroutine obtainLock(MI)
    type(MailboxInfo) :: MI
    logical :: file_exists
    integer :: io

    do
!        print *, "try get lock"
        inquire(file=MI%lock_fn, exist=file_exists)
!        print *, "file_exists: ", file_exists
        if (file_exists .eqv. .true.) then
            call sleep(1)
            cycle
        end if
        
        io = 0
        open(unit=MI%lock_fd, file=MI%lock_fn, form="formatted", access="stream", action="write", iostat=io)
        close(MI%lock_fd)

        if (io == 0) then
            exit
        else
            call sleep(1)
            cycle
        end if
    end do 
!    print *, "Lock got"

end subroutine

subroutine releaseLock(MI)
    type(MailboxInfo) :: MI
    call delFile(MI%lock_fn, MI%lock_fd)
end subroutine

subroutine delFile(fn, fd)
    implicit none
    integer :: fd
    character(len=*) :: fn
    
    open(unit=fd, file=fn, status="old")
    close(unit=fd, status="delete")
end subroutine

subroutine recv(MI, msg)
    implicit none
    type(MailboxInfo)  :: MI
    character(len=*) :: msg

    integer :: io
    logical :: file_exists

    do
        inquire(file=MI%recv_fn, exist=file_exists)
        if (file_exists .eqv. .true.) then
            exit
        else
            cycle
        end if
    end do

    call obtainLock(MI)
    
    io = 0
    open(unit=MI%recv_fd, file=MI%recv_fn, form="formatted", access="stream", action="read", iostat=io)
    
    read (MI%recv_fd, '(A)', iostat=io) msg
    close(MI%recv_fd)

    call delFile(MI%recv_fn, MI%recv_fd)

    call releaseLock(MI)
    
end subroutine


subroutine send(MI, msg)
    implicit none
    type(MailboxInfo)  :: MI
    character(len=*) :: msg

    integer :: io

    call obtainLock(MI)
    print *, "Lock get" 
    io = 0
    open(unit=MI%send_fd, file=MI%send_fn, form="formatted", access="stream", action="write", iostat=io)
    if (io /= 0) then
        print *, "Create send file iostat: ", io
    end if

    io = 0
    write (MI%send_fd, *, iostat=io) msg
    if (io /= 0) then
        print *, "Output send file iostat: ", io
    end if
    
    close(MI%send_fd)
    call releaseLock(MI)

end subroutine


subroutine hello(MI)
    implicit none
    type(MailboxInfo) :: MI
    character(256) :: msg

    call recv(MI, msg)
    msg = trim(msg)

    if ((msg .eq. "<<TEST>>") .or. (len(msg) .eq. len("<<TEST>>"))) then
        print *, "Recv hello!"
    else
        print *, len(msg), " : ", len("<<TEST>>")
        print *, "Weird msg: [", msg, "]"
    end if

    call send(MI, "<<TEST>>")

end subroutine



end module MailboxMod
