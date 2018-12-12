begin

    beg_yr = 1979
    end_yr = 2000

    yrs = end_yr - beg_yr + 1

    ntime = 12
	beg_time_i = ntime * (beg_yr - 1948)
    end_time_i = beg_time_i + ntime * yrs
    modes = 2

    i_20N = 0


    print("Reading file.")	
    fi_GEOP = addfile("/surtsey/ypeings/data/ncep2/monthly/GEOP_1948-2017.nc", "r")

    ;print(fi_GEOP->latitude)
   
 
    print("Load variable (20N - 90N)")

    obs = fi_GEOP->GEOP(level|0, latitude|i_20N:, longitude|:, time|beg_time_i:end_time_i-1)
    obs = dtrend_n(obs, True, 2)
    print(dimsizes(obs))

    print("Remove mean")
    obs = dim_rmvmean_n_Wrap(obs, (/2/))
    
    print("Remove seasonal cycle")

    do i = 0, 11
        mean = dim_avg_n(obs(:, :, i::12), 2)
        do j = 0, yrs - 1
            obs(:, :, j*12 + i) = obs(:, :, j*12 + i) - mean
        end do
    end do

    deg2rad = atan(1.0) / 45.0
    wgt = sqrt(cos(fi_GEOP->latitude(i_20N:) * deg2rad))

    ;print(dimsizes(obs))
    ;print(dimsizes(wgt))

    wgt_conform = conform(obs, wgt, 0)

    print("Add area weight")
    obsw = obs * wgt_conform



    print("Doing EOF")

    optEOF = False
    optEOF@jopt = 0

    ;EOFs = eofcov_Wrap(obsw, modes)
    EOFs = eofunc_Wrap(obsw, modes, optEOF)

    wgt_conform := conform(EOFs, wgt, 1)
    EOFs = EOFs / wgt_conform

    print("Eigenvalues")
    print(EOFs@eval)

    EOFs!0 = "modes"
    EOFs!1 = "latitude"
    EOFs!2 = "longitude"
    
    EOFs&modes     = ispan(0, modes-1, 1)
    EOFs&latitude  = fi_GEOP->latitude(i_20N:)
    EOFs&longitude = fi_GEOP->longitude

   ; output data

    fout_name = "EOFs.nc"
    system("rm -f " + fout_name)
    fout = addfile(fout_name, "c")
    fout->EOFs = EOFs
    fout->obs  = obs

    delete(fout)
    delete(fi_GEOP)

    print("Output file: " + fout_name)


    






end
