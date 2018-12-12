begin
    ; reference: https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices

    beg_yr = 1870
    end_yr = 2014
    ref_beg_yr = 1870


    yrs = end_yr - beg_yr + 1

    ntime = 12
	beg_time_i = ntime * (beg_yr - ref_beg_yr)
    end_time_i = beg_time_i + ntime * yrs - 1

    modes = 2

    j_20N  = 110
    j_65N  = 154 ; actually 64.5N
    i_115E = 295 ; actually 115.5E
    i_105W = 74  ; actually 105.5W 

    print("Reading file.")	
    fi = addfile("/surtsey/ypeings/data/hadisst/monthly/SST_1870-2017.nc", "r")

    print("Cropping...")
    ; obs = SST(latitude|j_20N:j_65N, longitude|i_105W:i_115E, time|beg_time_i:end_time_i)

    obs = fi->SST(beg_time_i:end_time_i, :, :)
    lat = fi->latitude
    lon = fi->longitude
    time = fi->time(beg_time_i:end_time_i)

    dims = dimsizes(obs)

    llat = conform_dims(dims(1:2), lat, 0)
    llon = conform_dims(dims(1:2), lon, 1)

    mmask = (llat .gt. 20.0) .and. (llat .lt. 65.0) .and. ((llon .lt. -105.0) .or. (llon .gt. 115.0))
    
    print("Masking data... ")
    obs = mask(obs, mmask, True)
    print("Detrending... (trend and mean removed)")
    obs = dtrend_msg_n(time, obs, True, False, 0)

    print("Remove seasonal cycle")
    do i = 0, 11
        mean = dim_avg_n(obs(i::12, :, :), 0)
        do j = 0, yrs-1
            obs(j*12 + i, :, :) = obs(j*12 + i, :, :) - mean
        end do
    end do


    deg2rad = atan(1.0) / 45.0
    wgt = sqrt(cos(lat * deg2rad))

    wgt_conform = conform(obs, wgt, 1)
    print("Add area weight")
    obsw = obs * wgt_conform
    obsw!0 = "time"
    obsw!1 = "latitude"
    obsw!2 = "longitude"

    print("Doing EOF")
    ; eofunc_n_Wrap is not available in NCL 6.3 yet
    ; so we must do reindexing

    optEOF = False
    optEOF@jopt = 0
    EOFs = eofunc_Wrap(obsw(latitude|:, longitude|:, time|:), modes, optEOF)

    wgt_conform := conform(EOFs, wgt, 1)
    EOFs = EOFs / wgt_conform

    print("Eigenvalues")
    print(EOFs@eval)

    EOFs!0 = "modes"
    EOFs!1 = "latitude"
    EOFs!2 = "longitude"
    
    EOFs&modes     = ispan(0, modes-1, 1)
    EOFs&latitude  = lat
    EOFs&longitude = lon

    ;output data
    fout_name = "PDO_EOFs.nc"
    system("rm -f " + fout_name)
    fout = addfile(fout_name, "c")
    fout->EOFs = EOFs
    fout->SST  = obs

    delete(fout)
    delete(fi)

    print("Output file: " + fout_name)
end
