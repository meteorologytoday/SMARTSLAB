module MLMML

#α   =
#β   =  
#c_p =
#ρ   =



struct OceanColumn
    N  :: Integer     # number of layers
    zs :: Array{Float64, 1} # position of (N+1) grid points
    S  :: Array{Float64, 1} # Salinity of N layers
    T  :: Array{Float64, 1} # Temperature of N layers
    h  :: Float64           # Mixed-layer depth
    γ  :: Float64

    function OceanColumn(zs::Array{Float64, 1})
        N = length(zs) - 1
        S = zeros(Float64, N)
        T = zeros(Float64, N)
        h = 0.0
        return new(N, zs, S, T, h)
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

    # 因為Gasper的正負號跟NiilerKraus不一致，需要重新推導
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

end
