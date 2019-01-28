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

    #println("### h: ", oc.h)
    #println("FLDO:", oc.FLDO)


    Δb = (oc.FLDO != -1) ? oc.b_ML - oc.bs[oc.FLDO] : 0.0

    # After convective adjustment, there still might
    # be some numerical error making Δb slightly negative
    # (the one I got is like -1e-15). So I set a tolarence
    # ( 0.001 K ≈ 3e-6 m/s^2 ).
    if Δb < 0.0 && abs(Δb) <= 3e-6
        Δb = 0.0
    end

    fric_u = getFricU(ua=ua)
    flag, val = calWeOrMLD(; h_ML=oc.h_ML, B=B0+J0, fric_u=fric_u, Δb=Δb) 
    #println("Before:" , oc.bs[10], "; oc.FLDO = ", oc.FLDO, "; Δb = ", Δb)

    # 1
    if flag == :MLD
        we = 0.0
        new_h_ML  = val
    elseif flag == :we
        we = val 
        new_h_ML = oc.h_ML + Δt * we
    end
    new_h_ML = boundMLD(new_h_ML; h_ML_max=min(h_ML_max, oc.zs[1]-oc.zs[end]))

    # 2
    # 3

    # ML
    #      i: Calculate integrated buoyancy that should
    #         be conserved purely through entrainment
    #     ii: Add to total buoyancy

    hb_new = getIntegratedBuoyancy(
        zs = oc.zs,
        bs = oc.bs,
        b_ML = oc.b_ML,
        h_ML = oc.h_ML,
        target_z = -new_h_ML
    )
  
    hb_chg_by_F = -(B0 + J0) * Δt

    #println(new_h, "; ", hb_new, ", ")
    new_b_ML = (hb_new + hb_chg_by_F) / new_h_ML
    
    OC_setMixedLayer!(
        oc,
        b_ML=new_b_ML,
        h_ML=new_h_ML
    )
    
    OC_doDiffusion_EulerBackward!(oc, Δt=Δt)
    if_convective_adjustment = OC_doConvectiveAdjustment!(oc)
    #println(oc.b_ML - oc.bs[oc.FLDO])
    return Dict(
        :flag => flag,
        :val  => val,
        :Δb   => Δb,
        :convective_adjustment => if_convective_adjustment
    )
end


