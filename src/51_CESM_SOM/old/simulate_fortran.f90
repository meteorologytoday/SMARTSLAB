
include "lib/MailboxMod.f90"


program exp_FIFO

use MailboxMod

implicit none
integer :: i
type(mbm_MailboxInfo) :: MI
character(1024)  :: msg

    call mbm_setDefault(MI)

    call mbm_hello(MI)
    print *, "Hello finish."

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: ", msg
        call mbm_send(MI, msg)
        print *, "Receiving message ... "
        call mbm_recv(MI, msg)
        
        print *, "Message received: ", trim(msg)

    end do


    msg = "<<END>>"
    call mbm_send(MI, msg)

    print *, "Program ends."

end program
