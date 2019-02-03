include "MailboxMod.f90"

module CommunicateMod

use MailboxMod

implicit none


contains


type BlackBox
    type(MailboxInfo) :: MI
    
end type BlackBox





subroutine init_external_model()
    




end subroutine


subroutine run_external_model()
end subroutine


subroutine end_external_model()
end subroutine

end module CommunicateMod
