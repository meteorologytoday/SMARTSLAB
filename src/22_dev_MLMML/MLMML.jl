module MLMML
using Printf

α   = 3e-4     # K^-1    http://www.kayelaby.npl.co.uk/general_physics/2_7/2_7_9.html
β   = 1e-3     # Simple estimation
c_p = 3985.0   # J / kg / K
ρ   = 1027.0   # kg / m^3
g   = 9.8      # m / s^2
h_min = 1.0
h_max = 1000.0
"""
    printConstants()

# Description
This function prints all constants used in this module.
"""
function printConstants()
    @printf("α   = %8.2f. Logrithmic expansion of density ρ as function of temperature T.\n", α)
    @printf("β   = %8.2f. Logrithmic expansion of density ρ as function of salinity S.\n", β)
    @printf("c_p = %8.2f J/kg/K. Specific heat of seawater.\n", c_p)
    @printf("ρ   = %8.2f kg/m^3. Mass density of seawater.\n", ρ)

end



mutable struct OceanColumn
    N      :: Integer           # Number of layers
    zs     :: Array{Float64, 1} # Position of (N+1) grid points
    bs     :: Array{Float64, 1} # Buoyancy of N layers
    Ks     :: Array{Float64, 1} # Diffusion coes between layers
    KML    :: Float64           # Diffusion coe of ML and FLDO
    b_ML   :: Float64
    h      :: Float64           # Mixed-layer depth
    FLDO   :: Integer           # First layer of deep ocean

    function OceanColumn(zs::Array{Float64, 1})
        N  = length(zs) - 1
        bs = zeros(Float64, N)
        Ks = zeros(Float64, N-1)
        KML = 0.0 
        b_ML = 0.0
        h  = h_min
        FLDO = 1
        return new(N, zs, bs, Ks, KML, b_ML, h, FLDO)
    end
end

"""

    calWeOrMLD(; h, B, fric_u, Δb, m=0.45, n=0.2)

# Description
This function returns either the entrainment speed ``w_e`` or diagnosed MLD (Mixed-Layer Depth) if calculated ``w_e < 0``.

``h`` is the current MLD, ``B`` the total buoyancy flux at the surface, ``Δb = b_m - b(-h)`` is the difference between ML and deep ocean which is assumed to be positive (buoyantly stable), ``m`` and ``n`` the parametrization constant of dissipation. The detailed derivation and origin of the default number follow following Gasper 1988. 

Here, we also assume that all sunlight is absorbed at the surface so that ``B`` is not a function of ``h``. This makes diagnose of ``h`` during shoaling much easier.


# Return values
This function returns a list with two elements. The first is a symbol. ``:we`` indicates the second value is the entrainment speed whereas ``:MLD`` indicates the second value is the diagnosed MLD.

"""
function calWeOrMLD(;
    h      :: Float64,
    B      :: Float64, 
    fric_u :: Float64,  
    Δb     :: Float64,
    m::Float64 = 0.45,
    n::Float64 = 0.20
)


    if Δb < 0
        throw(ErrorException("Δb cannot be negative."))
    end

    Term1 = 2.0 * m * fric_u^3.0
    Term2 = 0.5 * (B * (1.0 + n) - abs(B) * (1.0 - n))
    RHS = Term1 + h * Term2

    println("Term1: ", Term1, "; Term2:", Term2)

    if RHS > 0
        k2 = getTKE(fric_u)
        we = RHS / (h * Δb + k2)
        println(":we, h: ", h, "; Δb: ", Δb, "; B: ", B)
        return :we, we
    else
        # h becomes diagnostic. Notice that we assume
        # here that all sunlight is absorbed at the
        # very surface
       
        if Term2 == 0
            h_diag = h
        else
            h_diag = - Term1 / Term2
        end

        h_diag = h_diag

        return :MLD, h_diag
    end
end

function boundMLD(h)
    return max(min(h, h_max), h_min)
end


"""
    getTKE(fric_u)

# Description
This function returns the TKE (turbulent kinetic energy) `k = 0.5 * (v'^2)` of ML. This parameterization is given by Kim 1976: "A Generalized Bulk Model of the Oceanic Mixed Layer" in its equation (11)

"""
function getTKE(
    fric_u :: Float64
)
    cm = max(3e-2, 3.0 * fric_u)
    return 0.5 * cm^2.0
end

function getIntegratedBuoyancy(;
    zs       :: Array{Float64,1},
    bs       :: Array{Float64,1},
    b_ML     :: Float64,
    h        :: Float64,
    target_z :: Float64,
)

    if target_z < zs[end]
        throw(ErrorException("target_z cannot be deeper than the minimum of zs."))
    end

    if -target_z < h
        return b_ML * ( - target_z )
    end

    sum_b = 0.0
    sum_b += h * b_ML

    FLDO = getFLDO(zs=zs, h=h)
    
    if target_z > zs[FLDO+1]
        sum_b += bs[FLDO] * ( (-h) - target_z)
        return sum_b
    end
    
    sum_b += bs[FLDO] * ( (-h) - zs[FLDO+1]) 

    # Rest layers
    if FLDO < length(bs)
        for i = FLDO+1 : length(bs)
            if target_z < zs[i+1]
                sum_b += bs[i] * (zs[i] - zs[i+1])
            else
                sum_b += bs[i] * (zs[i] - target_z)
                return sum_b
            end
        end
    else
        return sum_b
    end

end

function getFLDO(;
    zs :: Array{Float64,1},
    h  :: Float64
)
    for i = 1:length(zs)-1
        if h <= (zs[1] - zs[i+1])
            return i
        end
    end

    return -1
end

function getWindStress(;
    u10::Float64
)

    return u10 * 1e-3 * ( (u10 < 25.0) 
                    ? 2.7 + 0.142 * u10 + 0.0764 * u10^2.0
                    : u10 * (2.16 + 0.5406 * (1.0 - exp(- (u10 - 25.0) / 7.5)))
    )

end


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

    Δb = oc.b_ML - oc.bs[oc.FLDO]
    fric_u = √(getWindStress(u10=ua) / ρ)
    flag, val = calWeOrMLD(; h=oc.h, B=B0+J0, fric_u=fric_u, Δb=Δb) 
    println("Before:" , oc.bs[10], "; oc.FLDO = ", oc.FLDO, "; Δb = ", Δb)

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

    println(new_h, "; ", hb_new, ", ")
    new_b_ML = (hb_new + hb_chg_by_F) / new_h

   
    println("new_h: ", new_h, "; new_b_ML: ", new_b_ML, "; hb_chg_by_F: ", hb_chg_by_F, "; -(B0+J0): ", -B0-J0)



    # Update profile
    bs_new = copy(oc.bs)
    if new_FLDO > 1
        bs_new[1:new_FLDO-1] .= new_b_ML
    end
 
    # Convective adjustment
    if_unstable = new_b_ML .< bs_new[new_FLDO]
    if any(if_unstable)
        conv_b_new = new_b_ML
        # Determine the least layers to be mixed
        bottom_layer_to_convect = 1
        for i = length(if_unstable):-1:1
            if if_unstable[i] == true
                bottom_layer_to_convect = i
                break
            end
        end

        # Determine how many extra layers should be mixed
        for i = bottom_layer_to_convect:length(bs_new) - 1
            conv_b_new = getIntegratedBuoyancy(
                zs       =  oc.zs,
                bs       =  bs_new,
                b_ML     =  new_b_ML,
                h        =  new_h,
                target_z =  oc.zs[i+1]
            ) / (-oc.zs[i+1])
           
            if conv_b_new < bs_new[i+1]
                bottom_layer_to_convect = i+1
            else
                break
            end
        end

        println("Convective Adjustment! new_FLDO: ", new_FLDO, "; update b from ", new_b_ML, " to ", conv_b_new)

        new_b_ML = conv_b_new
        new_FLDO = bottom_layer_to_convect + 1
        bs_new[1:new_FLDO-1] .= new_b_ML
    end
    


    #= 
    # Diffusion of all layers
    # b_flux[i] means the flux from layer i+1 to i (upward > 0)
    # the extra b_flux[end] is artificial for easier programming
    MLD_b_flux = - oc.KML * (new_b_ML - bs_new[new_FLDO])
    b_flux = zeros(Float64, length(oc.bs))
    for i = new_FLDO:length(b_flux)-1
       b_flux[i] = - oc.Ks[i] * (bs_new[i] - bs_new[i+1]) 
    end

    new_b_ML += MLD_b_flux / new_h * Δt
    if new_FLDO > 1
        bs_new[1:new_FLDO-1] .= new_b_ML
    end
      
    bs_new[new_FLDO] += (-MLD_b_flux + b_flux[new_FLDO]) / ((-new_h) - oc.zs[new_FLDO+1]) * Δt

    for i = new_FLDO+1:length(bs_new)-1
        bs_new[i] += (- b_flux[i-1] + b_flux[i]) / (oc.zs[i] - oc.zs[i+1]) * Δt
    end
    =#


    oc.bs[:] = bs_new
    oc.h = new_h
    oc.FLDO = new_FLDO 
    oc.b_ML = new_b_ML

    println("After:" , oc.bs[10])
end

end
