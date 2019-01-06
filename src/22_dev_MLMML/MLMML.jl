module MLMML
using Printf

α   = 1.0
β   = 1.0
c_p = 3985.0   # J / kg / K
ρ   = 1027.0   # kg / m^3

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



struct OceanColumn
    N      :: Integer           # Number of layers
    zs     :: Array{Float64, 1} # Position of (N+1) grid points
    bs     :: Array{Float64, 1} # Buoyancy of N layers
    h      :: Float64           # Mixed-layer depth
    FLDO   :: Float64           # First layer of deep ocean

    function OceanColumn(zs::Array{Float64, 1})
        N  = length(zs) - 1
        bs = zeros(Float64, N)
        h  = 0.0
        return new(N, zs, bs, h)
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
    if RHS > 0
        we = RHS / (h * Δb)
        return :we, we
    else
        # h becomes diagnostic. Notice that we assume
        # here that all sunlight is absorbed at the
        # very surface
        
        h_diag = - Term1 / Term2
        return :MLD, h_diag
    end
end

function getIntegratedBuoyancy(;
    zs       :: Array{Float64,1},
    bs       :: Array{Float64,1},
    h        :: Float64,
    target_z :: Float64,
)

    if target_z < zs[end]
        throw(ErrorException("target_z cannot be deeper than the minimum of zs."))
    end

    if -target_z < h
       return bs[1] * ( - target_z )
    end

    sum_b = 0.0
    sum_b += h * bs[1]

    FLDO = getFLDO(zs=zs, h=h)
    
    if target_z > zs[FLDO+1]
        sum_b += bs[FLDO] * ( (-h) - target_z)
        return sum_b
    end
    
    sum_b += bs[FLDO] * ( (-h) - zs[FLOD+1]) 

    # Rest layers
    for i = FLDO+1 : length(bs)
        if target_z < zs[i+1]
            sum_b += bs[i] * (zs[i] - zs[i+1])
        else
            sum_b += bs[i] * (zs[i] - target_z)
            return sum_b
        end
    end

end

function getFLDO(;
    zs :: Array{Float64,1},
    h  :: Float64
)
    for i = 1:length(zs)-1
        if h < (zs[i+1] - zs[1])
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

    # Find Δb
    Δb = b[1] - b[oc.DO_beg]
    fric_u = √(getWindStress(u10=ua) / ρ)
    flag, val = calWeOrMLD(; h=oc.h, B=B0+J0, fric_u=fric_u, Δb=Δb) 

    # 1
    if flag == :MLD
        we = 0.0
        new_h  = val
    else if flag == :we
        we = val 
        new_h += Δt * we
    end
    
    # 2
    new_FLDO = getFLDO(zs=oc.zs, h=new_h)

    # 3

    # ML
    #      i: Calculate integrated buoyancy that should
    #         be conserved purely through entrainment
    #     ii: Add to total buoyancy

    
    dbdt = - 1.0 / oc.h * (B0 + J0 + we * Δb)
    total_buoyancy_change = - (B0 + J0) * Δt
    
    # DO
    #      i: determine gradient of b flux
    #     ii: determine gradient of 
    for i 


    # 4

end

end
