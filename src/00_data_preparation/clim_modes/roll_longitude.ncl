begin
    in_filename = "/surtsey/ypeings/data/hadisst/monthly/SST_1870-2017.nc"
    out_filename = "SST_1870-2017_0-360.nc"
    print("Input  file: " + in_filename)
    print("Output file: " + out_filename)

    system("rm -f " + out_filename)
    fi = addfile(in_filename, "r")
    fo = addfile(out_filename, "c")

    lat = fi->latitude
    lon = fi->longitude
    time = fi->time

    expected_dim = (/dimsizes(time), dimsizes(lat), dimsizes(lon)/)

    ; find zero longitude
    lon0_i = -999
    n_lon = dimsizes(lon)

    do i = 0, dimsizes(lon)-1
        if lon(i) .ge. 0
            lon0_i = i
            break
        end if
    end do

    new_lon = lon
    new_lon(0            : n_lon-lon0_i-1) = lon(lon0_i : n_lon-1  )
    new_lon(n_lon-lon0_i : n_lon-1       ) = lon(0      : lon0_i-1 ) + 360.0
   
    varnames = getfilevarnames(fi)
    n_varnames = dimsizes(varnames)

    do i = 0, n_varnames - 1
        varname := varnames(i)
        print("Doing var: " + varname)
        var := fi->$varname$
        dim := dimsizes(var)
        if (dimsizes(dim) .ne. 3) .or. any(dim .ne. expected_dim)
            print("Var [" + varname + "] does not match have dimension expected. Skip this one.")
            continue
        end if
        
        new_var := var
        new_var(:, :, 0            : n_lon-lon0_i-1) = var(:, :, lon0_i : n_lon-1  )
        new_var(:, :, n_lon-lon0_i : n_lon-1       ) = var(:, :, 0      : lon0_i-1 )
        new_var&$new_var!2$ = new_lon 
        fo->$varname$ = new_var
    end do

    delete(fo)
    delete(fi)

    print("done.")
end
