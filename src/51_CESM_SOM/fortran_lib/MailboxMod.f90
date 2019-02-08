module MailboxMod
implicit none

type mbm_MailboxInfo
    Integer :: recv_cnt
    Integer :: send_cnt
    Integer :: recv_fd
    Integer :: send_fd
    Integer :: lock_fd

    character(len = 256) :: recv_fn
    character(len = 256) :: send_fn

    character(len = 256) :: lock_fn

    character(len = 256) :: log_file
end type


contains

integer function mbm_get_file_unit()
    
    integer :: lu, iostat
    logical :: opened
      
    do lu = 99, 1,-1
       inquire (unit=lu, opened=opened, iostat=iostat)
       if (iostat.ne.0) cycle
       if (.not.opened) exit
    end do
    
    mbm_get_file_unit = lu
    return
end function 

subroutine mbm_setDefault(MI)
    implicit none
    type(mbm_MailboxInfo) :: MI

    MI%recv_fn  = "mymodel2cesm.info"
    MI%send_fn  = "cesm2mymodel.info"
 
    MI%lock_fn  = "lock"
   
    MI%log_file = "log"

    MI%recv_cnt = 0
    MI%send_cnt = 0
 
    MI%recv_fd = mbm_get_file_unit()
    MI%send_fd = mbm_get_file_unit()
    MI%lock_fd = mbm_get_file_unit()
end subroutine 

subroutine mbm_appendPath(MI, path)
    implicit none
    type(mbm_MailboxInfo) :: MI
    character(len=256) :: path

    MI%recv_fn  = path // "/" // MI%recv_fn 
    MI%send_fn  = path // "/" // MI%send_fn 
    MI%lock_fn  = path // "/" // MI%lock_fn
    MI%log_file = path // "/" // MI%log_file

end subroutine 

subroutine mbm_obtainLock(MI)
    type(mbm_MailboxInfo) :: MI


    logical :: file_exists
    integer :: io

    do
        ! try to get lock
        inquire(file=MI%lock_fn, exist=file_exists)
        
        if (file_exists .eqv. .true.) then
            call sleep(1)
            cycle
        end if
        
        io = 0
        open(unit=MI%lock_fd, file=MI%lock_fn, form="formatted", access="stream", action="write", iostat=io)
        close(MI%lock_fd)

        ! If open file fails then try again
        if (io == 0) then
            exit
        else
            call sleep(1)
            cycle
        end if
    end do 

end subroutine

subroutine mbm_releaseLock(MI)
    type(mbm_MailboxInfo) :: MI
    call mbm_delFile(MI%lock_fn, MI%lock_fd)
end subroutine

subroutine mbm_delFile(fn, fd)
    implicit none
    integer :: fd
    character(len=*) :: fn
    logical :: file_exists

    inquire(file=fn, exist=file_exists)

    if (file_exists .eqv. .true.) then
        open(unit=fd, file=fn, status="old")
        close(unit=fd, status="delete")
    end if

end subroutine

subroutine mbm_clean(MI)
    implicit none
    type(mbm_MailboxInfo) :: MI

    call mbm_delFile(MI%recv_fn, MI%recv_fd)
    call mbm_delFile(MI%send_fn, MI%send_fd)
end subroutine


subroutine mbm_recv(MI, msg)
    implicit none
    type(mbm_MailboxInfo)  :: MI
    character(len=*) :: msg

    integer :: io
    logical :: file_exists

    do
        inquire(file=MI%recv_fn, exist=file_exists)
        if (file_exists .eqv. .true.) then
            exit
        else
            call sleep(1)
            cycle
        end if
    end do

    call mbm_obtainLock(MI)
    
    io = 0
    open(unit=MI%recv_fd, file=MI%recv_fn, form="formatted", access="stream", action="read", iostat=io)
    
    read (MI%recv_fd, '(A)', iostat=io) msg
    close(MI%recv_fd)
    
    msg = trim(msg)

    call mbm_delFile(MI%recv_fn, MI%recv_fd)

    call mbm_releaseLock(MI)
    
end subroutine


subroutine mbm_send(MI, msg)
    implicit none
    type(mbm_MailboxInfo)  :: MI
    character(len=*) :: msg

    integer :: io
    
    call mbm_obtainLock(MI)
    !print *, "Lock get" 
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
    call mbm_releaseLock(MI)

end subroutine


subroutine mbm_hello(MI)
    implicit none
    type(mbm_MailboxInfo) :: MI
    character(256) :: msg

    call mbm_recv(MI, msg)

    if ((msg .eq. "<<TEST>>") .and. (len(msg) .eq. len("<<TEST>>"))) then
        print *, "Recv hello!"
    else
        print *, len(msg), " : ", len("<<TEST>>")
        print *, "Weird msg: [", msg, "]"
    end if

    call mbm_send(MI, "<<TEST>>")

end subroutine

logical function mbm_messageCompare(msg1, msg2)
    implicit none
    character(*) :: msg1, msg2

    if (msg1 .eq. msg2) then
        mbm_messageCompare = .true.
    else
        mbm_messageCompare = .false.
    end if

end function


end module MailboxMod
