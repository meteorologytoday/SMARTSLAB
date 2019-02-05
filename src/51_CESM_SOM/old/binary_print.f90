include "lib/field_tools.f90"


program binary_print

use field_tools

implicit none

real(8) :: sst(24*19)
integer :: i

do i = 1 , size(sst)

    sst(i) = i*100

enddo

call write_1Dfield(5, "test_sst.bin", sst, size(sst))


end program
