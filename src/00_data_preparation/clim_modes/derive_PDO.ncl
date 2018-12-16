load "remap.ncl"
begin
    ; reference: https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices

    beg_yr = 1870
    end_yr = 2014
    ref_beg_yr = 1870


    yrs = end_yr - beg_yr + 1

    period = 12
	beg_time_i = period * (beg_yr - ref_beg_yr)
    end_time_i = beg_time_i + period * yrs - 1

    modes = 2

    print("Reading file.")	
    fi = addfile("./SST_1870-2017_0-360.nc", "r")

    print("Cropping...")

    obs = fi->SST(beg_time_i:end_time_i, :, :)
    lat = fi->latitude
    lon = fi->longitude
    time = fi->time(beg_time_i:end_time_i)

    ;xxx  := remap_5deg(obs)
    dims = dimsizes(obs)

    llat = conform_dims(dims(1:2), lat, 0)
    llon = conform_dims(dims(1:2), lon, 1)

    mmask = (llat .gt. 20.0) .and. (llat .lt. 65.0) .and. (llon .lt. 255.0) .and. (llon .gt. 115.0)
    
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




    ; output 5deg data
    EOFs_5deg := remap_5deg(EOFs)
    obs_5deg  := remap_5deg(obs)

    ; Get PDO index
    ntime = end_time_i - beg_time_i + 1
    PDO_idx = new(ntime, float)
    PDO_idx!0 = "time"
    PDO_idx&time = ispan(0, ntime - 1, 1)

    PDO_idx_5deg = PDO_idx
    do i = 0, ntime - 1
        PDO_idx(i) = sum(EOFs(0, :, :) * obs(beg_time_i + i, : ,:))
        PDO_idx_5deg(i) = sum(EOFs_5deg(0, :, :) * obs_5deg(beg_time_i + i, : ,:))
    end do



    ; output data
    fout_name := "PDO_EOFs.nc"
    system("rm -f " + fout_name)
    fout = addfile(fout_name, "c")
    fout->EOFs = EOFs
    fout->SST  = obs
    fout->PDO_INDEX  = PDO_idx
    delete(fout)
    print("Output file: " + fout_name)


    ; output 5deg data
    
    fout_name := "PDO_EOFs_5deg.nc"
    system("rm -f " + fout_name)
    fout = addfile(fout_name, "c")
    fout->EOFs = EOFs_5deg
    fout->SST  = obs_5deg
    fout->PDO_INDEX  = PDO_idx_5deg
    ;fout->xxx  = xxx
    delete(fout)
    print("Output file: " + fout_name)
    
    delete(fi)
end
