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
    h_ML   :: Float64,
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
    RHS = Term1 + h_ML * Term2

    #println("Term1: ", Term1, "; Term2:", Term2)

    if RHS > 0
        k = getTKE(fric_u=fric_u)
        we = RHS / (h_ML * Δb + k)
        #println(":we, h: ", h, "; Δb: ", Δb, "; B: ", B, "; k:", k)
        return :we, we
    else
        # h becomes diagnostic. Notice that we assume
        # here that all sunlight is absorbed at the
        # very surface
       
        if Term2 == 0
            h_ML_diag = h_ML
        else
            h_ML_diag = - Term1 / Term2
        end

        h_ML_diag = h_ML_diag

        return :MLD, h_ML_diag
    end
end


