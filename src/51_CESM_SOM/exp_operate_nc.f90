program exp_operate_nc
use netcdf_io
implicit none

character(256) :: filename = "test.nc"

integer :: i, j
!real(4) :: lat(19) = (/(i, i=-90,  90, 10)/)
!real(4) :: lon(24) = (/(i, i=  0, 345, 15)/)

real(4) :: lat(181) = (/(i, i=-90,  90, 1)/)
real(4) :: lon(360) = (/(i, i=  0, 359, 1)/)


integer :: nlat, nlon
integer :: flag

real(4), allocatable :: SST(:, :)



nlat = size(lat)
nlon = size(lon)

allocate(SST(nlon, nlat))

do j = 1, nlat 
    do i = 1, nlon
        SST(i, j) = i + j * nlon * 5
    end do
end do


flag = create_empty_netcdf(filename, nlat, lat, nlon, lon)
flag = write_netcdf(filename, nlat, nlon, SST, "SST")


contains








end program
