
include("MLMML.jl")


module SSM

using Printf
using Formatting
using ..MLMML

include("OceanColumnCollection.jl")

function stepOceanColumnCollection!(;
    occ   :: OceanColumnCollection,
    u     :: Array{Float64, 1},
    v     :: Array{Float64, 1},
    hflx  :: Array{Float64, 1},
    swflx :: Array{Float64, 1},
    Δt    :: Float64,

)

    ua   = (u.^2 + v.^2).^0.5
    hflx  *= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)
    swflx *= (MLMML.α * MLMML.g / MLMML.ρ / MLMML.c_p)
    
    for l = 1:occ.N_ocs
        
        MLMML.stepOceanColumn!(
            oc = occ.ocs[l],
            ua = ua[l],
            B0 = hflx[l],
            J0 = swflx[l],
            Δt = Δt,
        )
    end

end


function getSST!(;
    occ   :: OceanColumnCollection,
    sst   :: Array{Float64, 1}
)
    for l = 1:occ.N_ocs
      sst[l] = occ.ocs[l].b_ML / (MLMML.α * MLMML.g)
    end
end

function getInfo!(;
    occ   :: OceanColumnCollection,
    sst   :: Union{Array{Float64}, Nothing} = nothing,
    mld   :: Union{Array{Float64}, Nothing} = nothing,
)
    if mld != nothing
        for l = 1:occ.N_ocs
          mld[l] = occ.ocs[l].h_ML
        end
    end

    if sst != nothing

        for l = 1:occ.N_ocs
          sst[l] = occ.ocs[l].b_ML / (MLMML.α * MLMML.g)
        end

    end
end



end
