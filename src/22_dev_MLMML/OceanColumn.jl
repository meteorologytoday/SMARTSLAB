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
    oc.FLDO  = getFLDO(zs=oc.zs, h=h)
    oc.b_ML  = b_ML
    if oc.FLDO > 1
        oc.bs[1:oc.FLDO-1] .= oc.b_ML
    end
    
    if convective_adjustment
        doConvectiveAdjustment!(oc)
    end
end
