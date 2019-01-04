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


function calWe_Gasper1988(;
    h,
    B,       # Total turbulent buoyancy flux 
    fric_u,  
    Δb,
    m = 0.45,
    n = 0.2
)
    if Δb < 0
        throw(ErrorException("Δb cannot be negative."))
    end

    # 因為Gasper的正負號跟NiilerKraus不一致，需要重新推導
    RHS = 2.0 * m * fric_u^3.0 - 0.5 * h * ((1.0 - n) * )


end
