module field_tools

contains

subroutine write_1Dfield(fd, filename, f, nx)
character(len=*) :: filename
real(8), intent(in) :: f(nx)
integer, intent(in):: fd, nx
integer :: i,eflag


open (fd, file=filename, access="DIRECT", status='REPLACE', &
&       form='UNFORMATTED', recl=8*nx, iostat=eflag, convert='LITTLE_ENDIAN')

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if


write(fd,rec=1) (f(i),i=1,nx,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if

end subroutine

subroutine read_1Dfield(fd, filename, f, nx)
character(len=*) :: filename
real(8), intent(inout) :: f(nx)
integer, intent(in)    :: fd, nx
integer :: i, eflag


open (fd, file=filename, access="DIRECT", status='OLD', &
&       form='UNFORMATTED', recl=8*nx, iostat=eflag, convert='LITTLE_ENDIAN')

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
    print *, "Error number: ", eflag
end if

read(fd, rec=1) (f(i),i=1,nx,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
end if

end subroutine


subroutine read_2Dfield(fd, filename, f, nx, ny)
character(len=*) :: filename
real(4), intent(inout) :: f(nx,ny)
integer, intent(in)    :: fd, nx, ny
integer :: i,j,eflag


open (fd, file=filename, access="DIRECT", status='OLD', &
&       form='UNFORMATTED', recl=4*nx*ny, iostat=eflag)

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
    print *, "Error number: ", eflag
end if

read(fd, rec=1) ((f(i,j),i=1,nx,1),j=1,ny,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
end if

end subroutine


subroutine write_2Dfield(fd, filename, f, nx, ny)
character(len=*) :: filename
real(4), intent(in) :: f(nx,ny)
integer, intent(in):: fd, nx, ny
integer :: i,j,eflag


open (fd, file=filename, access="DIRECT", status='REPLACE', &
&       form='UNFORMATTED', recl=4*nx*ny, iostat=eflag)

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if

write(fd,rec=1) ((f(i,j),i=1,nx,1),j=1,ny,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if

end subroutine

subroutine read_3Dfield(fd, filename, f, nx, ny, nz)
character(len=*) :: filename
real(4), intent(inout) :: f(nx,ny,nz)
integer, intent(in)    :: fd, nx, ny, nz
integer :: i,j,k,eflag


open (fd, file=filename, access="DIRECT", status='OLD', &
&       form='UNFORMATTED', recl=4*nx*ny*nz, iostat=eflag)

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
end if

read(fd, rec=1) (((f(i,j,k),i=1,nx,1),j=1,ny,1),k=1,nz,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Reading field error. File name: ", trim(filename)
end if

end subroutine


subroutine write_3Dfield(fd, filename, f, nx, ny, nz)
character(len=*) :: filename
real(4), intent(in) :: f(nx,ny,nz)
integer, intent(in):: fd, nx, ny, nz
integer :: i,j,k,eflag


open (fd, file=filename, access="DIRECT", status='REPLACE', &
&       form='UNFORMATTED', recl=4*nx*ny*nz, iostat=eflag)

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if

write(fd,rec=1) (((f(i,j,k),i=1,nx,1),j=1,ny,1),k=1,nz,1)
close(fd)

if(eflag .ne. 0) then
    print *, "Writing field error. File name: ", trim(filename)
end if

end subroutine


end module field_tools
