using LinearAlgebra

"""

function doDiffusion!(oc::OceanColumn)

# Description

This function modifies `bs` and `b_ML`

"""
function OC_doDiffusion!(oc::OceanColumn; Δt::Float64)

new_b_ML = doDiffusion!(
    hs=oc.hs,
    zs=oc.zs,
    Δzs=oc.Δzs,
    bs=oc.bs,
    b_ML=oc.b_ML,
    h_ML=oc.h_ML,
    K=oc.K,
    FLDO=oc.FLDO,
    Δt=Δt
)

oc.b_ML = new_b_ML

end

"""
function doDiffusion!(;hs, zs, Δzs, bs, b_ML, h, K, FLDO, Δt)

# Description

This function modifies `bs`, and returns new b_ML.

"""
function doDiffusion!(;
    hs   :: Array{Float64, 1}, 
    zs   :: Array{Float64, 1}, 
    Δzs  :: Array{Float64, 1}, 
    bs   :: Array{Float64, 1}, 
    b_ML :: Float64,
    h_ML :: Float64,
    K    :: Float64,
    FLDO :: Integer,
    Δt   :: Float64,
)

    # Diffusion of all layers
    # b_flux[i] means the flux from layer i+1 to i (upward > 0)
    # the extra b_flux[end] is artificial for easier programming

    FLDO_h = -zs[FLDO+1] - h_ML
    ML_FLDO_Δz = min((h_ML + FLDO_h) / 2.0, FLDO_h)
    stable_diffusion = checkDiffusionStability(Δz=ML_FLDO_Δz, K=K, Δt=Δt) 


    b_flux = zeros(Float64, length(bs))
    for i = FLDO:length(b_flux)-1
       b_flux[i] = K * (bs[i+1] - bs[i]) / Δzs[i]
    end

    if stable_diffusion
        ML_b_flux = K * (bs[FLDO] - b_ML) / ML_FLDO_Δz
    else
        ML_jump_Δz = min((h_ML + hs[FLDO+1])/2.0, hs[FLDO+1])
        ML_b_flux = K * (bs[FLDO+1] - b_ML) / ML_jump_Δz
        b_flux[FLDO] = ML_b_flux
    end

    # Update buoyancy        
    for i = FLDO+1:length(bs)-1
        bs[i] += (b_flux[i] - b_flux[i-1]) / hs[i] * Δt
    end

    if stable_diffusion
        new_b_ML += ML_b_flux / h_ML * Δt
        bs[FLDO] += (b_flux[FLDO] - ML_b_flux) / FLDO_h * Δt
    else
        # let bs[FLDO] be (b_ML + bs[FLDO+1]) / 2 and 
        # integrated buoyancy of ML and FLDO remains the same.

        A = h_ML * b_ML + FLDO_h * bs[FLDO]
        new_b_ML = (A - FLDO_h * bs[FLDO+1] / 2.0) / (h_ML + FLDO_h / 2.0)
        bs[FLDO] = (new_b_ML + bs[FLDO+1]) / 2.0
    end

    if FLDO > 1
        bs[1:FLDO-1] .= new_b_ML
    end

    return new_b_ML
end



function OC_doDiffusion_EulerBackward!(oc::OceanColumn; Δt::Float64)
    
    new_b_ML = doDiffusion_BackwardEuler!(
        zs=oc.zs,
        bs=oc.bs,
        b_ML=oc.b_ML,
        h_ML=oc.h_ML,
        K=oc.K,
        FLDO=oc.FLDO,
        Δt=Δt,
    )

    oc.b_ML = new_b_ML
end


function doDiffusion_BackwardEuler!(;
    zs   :: Array{Float64, 1}, 
    bs   :: Array{Float64, 1}, 
    b_ML :: Float64,
    h_ML :: Float64,
    K    :: Float64,
    FLDO :: Integer,
    Δt   :: Float64,
)
    if FLDO == -1
        return b_ML
    end

    # Diffusion of all layers
    # b_flux[i] means the flux from layer i+1 to i (upward > 0)
    # the extra b_flux[end] is artificial for easier programming
    n   = length(bs) - (FLDO - 1) + 1
    Δzs = zeros(Float64, n-1)
    hs  = zeros(Float64, n)
    bs_RHS = zeros(Float64, n)

    bs_RHS[1] = b_ML
    bs_RHS[2:end] = bs[FLDO:end]

    #println("bs[end]", bs[end])

    #println("length(hs) = ", length(hs))
    #println("length(zs) = ", length(zs))
    #println("length(Δzs) = ", length(Δzs))

    hs[1] = h_ML
    hs[2] = -zs[FLDO+1] - h_ML
    hs[3:end] = zs[FLDO+1:end-1] - zs[FLDO+2:end]

    #println(hs)

    Δzs[1] = max(min((hs[1] + hs[2]) / 2.0, hs[2]), Δt * K / h_ML_min)
    Δzs[2:end] = (hs[2:end-1] + hs[3:end]) / 2.0

    αs = Δt * K ./ hs

    A = spzeros(Float64, n, n)

    A[1, 1] = 1.0 + αs[1] / Δzs[1]
    A[1, 2] = - αs[1] / Δzs[1]

    for i=2:n-1
        A[i, i-1] = - αs[i] / Δzs[i-1]
        A[i, i  ] = 1.0 + αs[i] * (Δzs[i] + Δzs[i-1]) / (Δzs[i] * Δzs[i-1])
        A[i, i+1] = - αs[i] / Δzs[i]
    end

    A[n, n  ] = 1.0 + αs[n] / Δzs[n-1]
    A[n, n-1] = - αs[n] / Δzs[n-1]

    #=
    println("n: ", n)
    if n <=10
        println("Det: ", det(A))
        println("A: ", A)
    end
    =#
    bs_new = A \ bs_RHS
    new_b_ML = bs_new[1]
    if FLDO > 1
        bs[1:FLDO-1] .= new_b_ML
    end
    bs[FLDO:end] = bs_new[2:end]

    return new_b_ML
end
