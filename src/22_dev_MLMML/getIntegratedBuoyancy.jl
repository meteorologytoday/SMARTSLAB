function getIntegratedBuoyancy(;
    zs       :: Array{Float64,1},
    bs       :: Array{Float64,1},
    b_ML     :: Float64,
    h        :: Float64,
    target_z :: Float64,
)

    if target_z < zs[end]
        throw(ErrorException("target_z cannot be deeper than the minimum of zs."))
    end

    if -target_z < h
        return b_ML * ( - target_z )
    end

    sum_b = 0.0
    sum_b += h * b_ML

    FLDO = getFLDO(zs=zs, h=h)
    
    if target_z > zs[FLDO+1]
        sum_b += bs[FLDO] * ( (-h) - target_z)
        return sum_b
    end
    
    sum_b += bs[FLDO] * ( (-h) - zs[FLDO+1]) 

    # Rest layers
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

