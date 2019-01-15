mutable struct OceanColumn
    N      :: Integer           # Number of layers
    zs     :: Array{Float64, 1} # Position of (N+1) grid points
    bs     :: Array{Float64, 1} # Buoyancy of N layers
    K      :: Float64           # Diffusion coes
    b_ML   :: Float64
    h      :: Float64           # Mixed-layer depth
    FLDO   :: Integer           # First layer of deep ocean

    # Derived quantities
    hs     :: Array{Float64, 1} # Thickness of layers
    Δzs    :: Array{Float64, 1} # Δz between layers

    function OceanColumn(;
        N      :: Integer,
        zs     :: Array{Float64, 1},
        bs     :: Array{Float64, 1} ,
        K      :: Float64,
        b_ML   :: Float64,
        h      :: Float64,
        FLDO   :: Integer
    )
        hs  = zs[1:end-1] - zs[2:end]
        Δzs = (hs[1:end-1] + hs[2:end]) / 2.0
        return new(N, zs, bs, K, b_ML, h, FLDO, hs, Δzs)
    end
end

function copy(oc::OceanColumn)
    return OceanColumn(
        N=oc.N,
        zs=Base.copy(oc.zs),
        bs=Base.copy(oc.bs),
        K=oc.K,
        b_ML=oc.b_ML,
        h=oc.h,
        FLDO=oc.FLDO
    )
end


function setBuoyancy!(
    oc  ::OceanColumn;
    bs  ::Array{Float64,1},
    b_ML::Float64,
    h   ::Float64,
    convective_adjustment::Bool=true
)

    if h < h_min
        throw(ErrorException(Formatting("h cannot be less than h_min: {:.2f}", h_min)))
    end
    
    oc.bs[:] = bs
    setMixedLayer!(oc; b_ML=b_ML, h=h, convective_adjustment=false)

    if convective_adjustment
        doConvectiveAdjustment!(oc)
    end
end

function setMixedLayer!(
    oc  ::OceanColumn;
    b_ML::Float64,
    h   ::Float64,
    convective_adjustment::Bool=true
)

    if h < h_min
        throw(ErrorException(Formatting("h cannot be less than h_min: {:.2f}", h_min)))
    end

    oc.h = h
    oc.FLDO  = getFLDO(zs=oc.zs, h=h)
    oc.b_ML  = b_ML
    if oc.FLDO > 1
        oc.bs[1:oc.FLDO-1] .= oc.b_ML
    end
    
    if convective_adjustment
        doConvectiveAdjustment!(oc)
    end
end


function makeBlankOceanColumn(;zs::Array{Float64, 1})
    N     = length(zs) - 1
    bs    = zeros(Float64, N)
    K     = 0.0
    b_ML  = 0.0
    h     = h_min
    FLDO  = 1

    oc = OceanColumn(N=N, zs=zs, bs=bs, K=K, b_ML=b_ML, h=h, FLDO=FLDO)
    updateFLDO!(oc)

    return oc
end

function makeSimpleOceanColumn(;
    zs      :: Array{Float64, 1},
    b_slope :: Float64 = 30.0 / 5000.0 * g * α,
    b_ML    :: Float64 = 1.0,
    h       :: Float64 = h_min,
    Δb      :: Float64 = 0.0,
    K       :: Float64 = 1e-5,
)

oc = makeBlankOceanColumn(zs=zs)

bs = zeros(Float64, length(zs)-1)
for i = 1:length(bs)
    z = (zs[i] + zs[i+1]) / 2.0
    if z > -h
        bs[i] = b_ML
    else
        bs[i] = b_ML - Δb - b_slope * (-z - h)
    end
end

setBuoyancy!(oc, bs=bs, b_ML=b_ML, h=h)
oc.K = K
updateFLDO!(oc)

return oc
end
