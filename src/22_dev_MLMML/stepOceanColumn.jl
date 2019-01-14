"""

    stepOceanColumn!(;
        oc  :: OceanColumn,
        ua  :: Float64,
        B0  :: Float64,
        J0  :: Float64,
        Δt  :: Float64 
    )

# Description
This function update the OceanColumn forward in time.

"""
function stepOceanColumn!(;
    oc  :: OceanColumn,
    ua  :: Float64, # Currently assumed to be u10
    B0  :: Float64,
    J0  :: Float64,
    Δt  :: Float64 
)
    # Pseudo code
    # Current using only Euler forward scheme:
    # 1. Determine h at t+Δt
    # 2. Determine how many layers are going to be
    #    taken away by ML.
    # 3. Cal b at t+Δt for both ML and DO
    # 4. Detect if it is buoyantly stable.
    #    Correct it (i.e. convection) if it is not.
    # 5. If convection happens, redetermine h.

    # p.s.: Need to examine carefully about the
    #       conservation of buoyancy in water column

    println("### h: ", oc.h)
    Δb = oc.b_ML - oc.bs[oc.FLDO]
    fric_u = getFricU(ua=ua)
    flag, val = calWeOrMLD(; h=oc.h, B=B0+J0, fric_u=fric_u, Δb=Δb) 
    #println("Before:" , oc.bs[10], "; oc.FLDO = ", oc.FLDO, "; Δb = ", Δb)

    # 1
    if flag == :MLD
        we = 0.0
        new_h  = val
    elseif flag == :we
        we = val 
        new_h = oc.h + Δt * we
    end
    new_h = boundMLD(new_h)
    println("flag: ", flag, "; val:", val)    
    # 2
    new_FLDO = getFLDO(zs=oc.zs, h=new_h)
    #println("new_FLDO: ", new_FLDO)
    # 3

    # ML
    #      i: Calculate integrated buoyancy that should
    #         be conserved purely through entrainment
    #     ii: Add to total buoyancy

    hb_new = getIntegratedBuoyancy(
        zs = oc.zs,
        bs = oc.bs,
        b_ML = oc.b_ML,
        h  = oc.h,
        target_z = -new_h
    )
  
    hb_chg_by_F = -(B0 + J0) * Δt

    #println(new_h, "; ", hb_new, ", ")
    new_b_ML = (hb_new + hb_chg_by_F) / new_h

   
    #println("new_h: ", new_h, "; new_b_ML: ", new_b_ML, "; hb_chg_by_F: ", hb_chg_by_F, "; -(B0+J0): ", -B0-J0)

    # Update profile
    bs_new = Base.copy(oc.bs)
    if new_FLDO > 1
        bs_new[1:new_FLDO-1] .= new_b_ML
    end

    

    new_h, new_b_ML, new_FLDO = doConvectiveAdjustment!(
        zs   = oc.zs,
        bs   = bs_new,
        h    = new_h,
        b_ML = new_b_ML,
        FLDO = new_FLDO,
    ) 
 
     
    # Diffusion of all layers
    # b_flux[i] means the flux from layer i+1 to i (upward > 0)
    # the extra b_flux[end] is artificial for easier programming
    FLDO_h = -oc.zs[new_FLDO+1] - new_h
    MLD_b_flux = - oc.K_ML * (new_b_ML - bs_new[new_FLDO]) / ( FLDO_h / 2.0)

    b_flux = zeros(Float64, length(oc.bs))
    for i = new_FLDO:length(b_flux)-1
       b_flux[i] = - oc.K_DO * (bs_new[i] - bs_new[i+1]) / oc.Δzs[i] 
    end

    new_b_ML += MLD_b_flux / new_h * Δt
    if new_FLDO > 1
        bs_new[1:new_FLDO-1] .= new_b_ML
    end
      
    bs_new[new_FLDO] += (-MLD_b_flux + b_flux[new_FLDO]) / ((-new_h) - oc.zs[new_FLDO+1]) * Δt

    for i = new_FLDO+1:length(bs_new)-1
        bs_new[i] += (- b_flux[i-1] + b_flux[i]) / oc.hs[i] * Δt
    end
 
    println("FLDO_h: ", FLDO_h)
    println("MLD_b_flux: ", MLD_b_flux)
    println("b_flux[new_FLDO]: ", b_flux[new_FLDO])

   


    oc.bs[:] = bs_new
    oc.h = new_h
    oc.FLDO = new_FLDO 
    oc.b_ML = new_b_ML

    return Dict(
        :flag => flag,
        :val  => val,
        :Δb   => Δb,
    )
end


