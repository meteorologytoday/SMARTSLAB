function getIntegratedBuoyancy(oc::OceanColumn; target_z::Float64 = NaN)

    if isnan(target_z)
        target_z = oc.zs[end]
    end

    return getIntegratedBuoyancy(
        zs = oc.zs,
        bs = oc.bs,
        b_ML = oc.b_ML,
        h_ML = oc.h_ML,
        target_z = target_z
    )
end


function getIntegratedBuoyancy(;
    zs       :: Array{Float64,1},
    bs       :: Array{Float64,1},
    b_ML     :: Float64,
    h_ML     :: Float64,
    target_z :: Float64,
)

    if target_z < zs[end]
        throw(ErrorException("target_z cannot be deeper than the minimum of zs."))
    end


    # Integrate mixed layer
    if -target_z < h_ML
        return b_ML * ( - target_z )
    end

    sum_b = 0.0
    sum_b += h_ML * b_ML


    # Test if entire ocean column is mixed layer
    FLDO = getFLDO(zs=zs, h_ML=h_ML)
    if FLDO == -1
        return sum_b
    end

    # Integrate FLDO
    if target_z > zs[FLDO+1]
        sum_b += bs[FLDO] * ( (-h_ML) - target_z)
        return sum_b
    end
    
    sum_b += bs[FLDO] * ( (-h_ML) - zs[FLDO+1]) 

    # Integrate rest layers
    if FLDO < length(bs)
        for i = FLDO+1 : length(bs)
            if target_z < zs[i+1]
                sum_b += bs[i] * (zs[i] - zs[i+1])
            else
                sum_b += bs[i] * (zs[i] - target_z)
                return sum_b
            end
        end
    else
        return sum_b
    end

end

