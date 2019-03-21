
include("MLMML.jl")


module SSM

using Printf
using Formatting
using ..MLMML

missing_value = 1e20


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

        if occ.mask[l] == 0.0
            continue
        end        

        MLMML.stepOceanColumn!(
            oc = occ.ocs[l],
            ua = ua[l],
            B0 = hflx[l],
            J0 = swflx[l],
            Δt = Δt,
        )
    end

end

function maskData!(occ::OceanColumnCollection, arr::Array{Float64})
    for i = 1:occ.N_ocs
        if occ.mask[i] == 0.0
            arr[i] = missing_value
        end
    end
end

function getInfo!(;
    occ   :: OceanColumnCollection,
    sst   :: Union{Array{Float64}, Nothing} = nothing,
    mld   :: Union{Array{Float64}, Nothing} = nothing,
)
    if mld != nothing
        for l = 1:occ.N_ocs
            if occ.mask[l] == 0.0
                continue
            end
            mld[l] = occ.ocs[l].h_ML
        end
    end

    if sst != nothing

        for l = 1:occ.N_ocs

          if occ.mask[l] == 0.0
              continue
          end
          sst[l] = occ.ocs[l].b_ML / (MLMML.α * MLMML.g)
        end

    end
end



end
