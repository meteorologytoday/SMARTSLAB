mutable struct OceanColumn
    N      :: Integer           # Number of layers
    zs     :: Array{Float64, 1} # Position of (N+1) grid points
    bs     :: Array{Float64, 1} # Buoyancy of N layers
    K      :: Float64           # Diffusion coes
    b_ML   :: Float64
    h_ML   :: Float64           # Mixed-layer depth
    FLDO   :: Integer           # First layer of deep ocean

    # Derived quantities
    hs     :: Array{Float64, 1} # Thickness of layers
    Δzs    :: Array{Float64, 1} # Δz between layers

    function OceanColumn(;
        N      :: Integer,
        zs     :: Array{Float64, 1},
        bs     :: Array{Float64, 1},
        K      :: Float64,
        b_ML   :: Float64,
        h_ML   :: Float64,
        FLDO   :: Integer,
        hs     :: Union{Array{Float64, 1}, Nothing} = nothing,
        Δzs    :: Union{Array{Float64, 1}, Nothing} = nothing,
    )

        if hs == nothing
            hs  = zs[1:end-1] - zs[2:end]
        end

        if Δzs == nothing
            Δzs = (hs[1:end-1] + hs[2:end]) / 2.0
        end

        return new(N, zs, Base.copy(bs), K, b_ML, h_ML, FLDO, hs, Δzs)
    end
end

function copy(oc::OceanColumn)
    return OceanColumn(
        N=oc.N,
        zs=Base.copy(oc.zs),
        bs=Base.copy(oc.bs),
        K=oc.K,
        b_ML=oc.b_ML,
        h_ML=oc.h_ML,
        FLDO=oc.FLDO
    )
end

#    if h < h_min
#        throw(ErrorException(Formatting("h cannot be less than h_min: {:.2f}", h_min)))
#    end
 
function OC_setBuoyancy!(
    oc  ::OceanColumn;
    bs  ::Array{Float64,1},
    b_ML::Float64,
    h_ML::Float64,
)

   
    oc.bs[:] = bs
    OC_setMixedLayer!(oc; b_ML=b_ML, h_ML=h_ML)

end


function setMixedLayer!(;
    bs  :: Array{Float64, 1},
    zs  :: Array{Float64, 1},
    b_ML::Float64,
    h_ML::Float64
)
    FLDO  = getFLDO(zs=zs, h_ML=h_ML)

    if FLDO > 1
        bs[1:FLDO-1] .= b_ML
    elseif FLDO == -1
        bs[:] .= b_ML
    end
   
    return FLDO 
end

function OC_setMixedLayer!(
    oc  ::OceanColumn;
    b_ML::Float64,
    h_ML::Float64,
)

    oc.h_ML  = h_ML
    oc.b_ML  = b_ML
    oc.FLDO  = setMixedLayer!(bs=oc.bs, zs=oc.zs, b_ML=b_ML, h_ML=h_ML)

end


function makeBlankOceanColumn(;zs::Array{Float64, 1})
    N     = length(zs) - 1
    bs    = zeros(Float64, N)
    K     = 0.0
    b_ML  = 0.0
    h_ML  = h_ML_min
    FLDO  = 1

    oc = OceanColumn(N=N, zs=zs, bs=bs, K=K, b_ML=b_ML, h_ML=h_ML, FLDO=FLDO)
    OC_updateFLDO!(oc)

    return oc
end

function makeSimpleOceanColumn(;
    zs      :: Array{Float64, 1},
    b_slope :: Float64 = 30.0 / 5000.0 * g * α,
    b_ML    :: Float64 = 1.0,
    h_ML    :: Float64 = h_ML_min,
    Δb      :: Float64 = 0.0,
    K       :: Float64 = 1e-5
)

oc = makeBlankOceanColumn(zs=zs)

bs = zeros(Float64, length(zs)-1)
for i = 1:length(bs)
    z = (zs[i] + zs[i+1]) / 2.0
    if z > -h_ML
        bs[i] = b_ML
    else
        bs[i] = b_ML - Δb - b_slope * (-z - h_ML)
    end
end

OC_setBuoyancy!(oc, bs=bs, b_ML=b_ML, h_ML=h_ML)
oc.K = K
OC_updateFLDO!(oc)

return oc
end
