module MLMML

#α   =
#β   =  
#c_p =
#ρ   =



struct OceanColumn{G <: AbstractFloat}
    N  :: Integer     # number of layers
    zs :: Array{G, 1} # position of (N+1) grid points
    S  :: Array{G, 1} # Salinity of N layers
    T  :: Array{G, 1} # Temperature of N layers
    h  :: Float64           # Mixed-layer depth

    function OceanColumn(zs::Array{G, 1}) where G <: AbstractFloat
        N = length(zs) - 1
        S = zeros(G, N)
        T = zeros(G, N)
        h = 0.0
        return new{G}(N, zs, S, T, h)
    end
end

end
