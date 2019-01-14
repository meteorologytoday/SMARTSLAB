module MLMML
using Printf
using Formatting
using SparseArrays

include("constants.jl")
include("OceanColumn.jl")
include("calWeOrMLD.jl")
include("doConvectiveAdjustment.jl")
include("getIntegratedBuoyancy.jl")
include("stepOceanColumn.jl")


"""
This function checks if CFL criteria is satisfied which is required by Euler Forward Scheme. Explicitly,

 K Δt       1
------  <= ---
(Δz)^2      2

for every layer. This function returns true every layer is satisfied, returns false if any of the layers is not.

"""
function checkDiffusionStability(;
    Δz :: Array{Float64, 1},
    K  :: Float64,
    Δt :: Float64,
)

    return all( Δz .>= √(2.0 * K * Δt) )
end

function minΔz(;
    K :: Float64,
    Δt:: Float64,
)
    return √(2.0 * K * Δt)
end



function boundMLD(h)
    return max(min(h, h_max), h_min)
end


"""
    getTKE(fric_u)

# Description
This function returns the TKE (turbulent kinetic energy) `k = 0.5 * (v'^2)` of ML. This parameterization is given by Kim 1976: "A Generalized Bulk Model of the Oceanic Mixed Layer" in its equation (11)

"""
function getTKE(;
    fric_u :: Float64
)
    cm = max(3e-2, 3.0 * fric_u)
    return  0.5 * cm^2.0
end


function updateFLDO!(oc::OceanColumn)
    oc.FLDO = getFLDO(zs=oc.zs, h=oc.h)
end

function getFLDO(;
    zs :: Array{Float64,1},
    h  :: Float64
)
    for i = 1:length(zs)-1
        if h < (zs[1] - zs[i+1])  # I don't use equality in order to avoid Δb = 0 during some initialization
            return i
        end
    end

    if h == (zs[1] - zs[end])
        return length(zs)-1
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

function getFricU(;
    ua::Float64
)
    return √(getWindStress(u10=ua) / ρ)
end

end
