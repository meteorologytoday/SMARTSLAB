
include "MailboxMod.f90"


program exp_FIFO

use MailboxMod

implicit none
integer :: i
type(MailboxInfo) :: MI
character(1024)  :: msg

    call setDefault(MI)

    call hello(MI)
    print *, "Hello finish."

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: ", msg
        call send(MI, msg)
        print *, "Receiving message ... "
        call recv(MI, msg)
        
        print *, "Message received: ", trim(msg)

    end do


    msg = "<<END>>"
    call send(MI, msg)

    print *, "Program ends."

end program
