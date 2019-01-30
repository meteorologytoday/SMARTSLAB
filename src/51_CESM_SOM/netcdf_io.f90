module netcdf_io
use netcdf
implicit none


contains

subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
end subroutine check  

integer function create_empty_netcdf(filename, nlat, lat, nlon, lon)
    implicit none
    character(len=*) :: filename
    real(4)         :: lon(nlon), lat(nlat)
    integer :: nlat, nlon

    integer :: ncid, lon_dimid, lat_dimid, lon_varid, lat_varid

    call check( nf90_create(filename, NF90_CLOBBER, ncid) )
    call check( nf90_def_dim(ncid, "lon", nlon, lon_dimid) )
    call check( nf90_def_dim(ncid, "lat", nlat, lat_dimid) )

    call check( nf90_def_var(ncid, "lon", NF90_FLOAT, lon_dimid, lon_varid) )
    call check( nf90_def_var(ncid, "lat", NF90_FLOAT, lat_dimid, lat_varid) )

    call check( nf90_enddef(ncid) )

    call check( nf90_put_var(ncid, lon_varid, lon) )
    call check( nf90_put_var(ncid, lat_varid, lat) )
    
    call check( nf90_close(ncid) )

    print *, "*** SUCCESS create empty netCDF file: ", filename
    
    create_empty_netcdf = 0

end function create_empty_netcdf


! Append to existing file directly
integer function write_netcdf(filename, nlat, nlon, var, varname)
    character(len=*) :: filename, varname
    integer :: nlat, nlon
    real(4) :: var(nlon, nlat)

    
    integer :: ncid, varid, lat_varid, lon_varid, dimids(2)

    call check(nf90_open(filename, NF90_WRITE, ncid))
    call check(nf90_inq_varid(ncid, "lat", lat_varid))
    call check(nf90_inq_varid(ncid, "lon", lon_varid))

    dimids = (/lon_varid, lat_varid/)

    call check(nf90_redef(ncid))

    call check( nf90_def_var(ncid, varname, NF90_FLOAT, dimids, varid) )

    call check( nf90_enddef(ncid) )

    call check( nf90_put_var(ncid, varid, var) )
    
    call check( nf90_close(ncid) )

    write_netcdf = 0

end function

end module netcdf_io
