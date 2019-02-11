
mutable struct OceanColumnCollection
    N_ocs    :: Integer           # Number of columns
    N        :: Integer           # Number of layers
    zs       :: Array{Float64, 1} # Position of (N+1) grid points
    K        :: Float64           # Diffusion coes
    mask     :: Array{Float64}
    mask_idx :: Any

    # Derived quantities
    hs     :: Array{Float64, 1} # Thickness of layers
    Δzs    :: Array{Float64, 1} # Δz between layers

    
    ocs    :: Array{MLMML.OceanColumn, 1}

    function OceanColumnCollection(;
        N_ocs  :: Integer,
        N      :: Integer,
        zs     :: Array{Float64, 1},
        bs     :: Array{Float64, 1},
        K      :: Float64,
        b_ML   :: Float64,
        h_ML   :: Float64,
        FLDO   :: Integer,
        mask   :: Union{Array{Float64}, Nothing} = nothing,
    )

        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0

        if mask == nothing
            mask = zeros(Float64, N_ocs)
            mask .+= 1.0
        end

        mask_idx = (mask .== 0.0)
        ocs = Array{MLMML.OceanColumn}(undef, N_ocs)

        for i=1:N_ocs

            if mask[i] == 0.0
                continue
            end

            ocs[i] = MLMML.OceanColumn(
                N = N,
                zs = zs,
                bs = bs,
                K  = K,
                b_ML = b_ML,
                h_ML = h_ML,
                FLDO = FLDO,
                hs = hs,
                Δzs = Δzs
            )
        end 

        return new(N_ocs, N, zs, K, mask, mask_idx, hs, Δzs, ocs)
    end

end


