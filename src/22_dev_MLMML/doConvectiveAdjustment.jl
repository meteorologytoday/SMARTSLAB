function doConvectiveAdjustment!(oc::OceanColumn)
    oc.h, oc.b_ML, oc.FLDO = doConvectiveAdjustment!(
        zs   = oc.zs,
        bs   = oc.bs,
        h    = oc.h,
        b_ML = oc.b_ML,
        FLDO = oc.FLDO,
    ) 
end


"""

This function only do convective adjustment for the upper most mixed layer.
It searches for the lowest layer that has larger buoyancy than mixed-layer then mixed all layers above it.

"""
function doConvectiveAdjustment!(;
    zs   :: Array{Float64, 1},
    bs   :: Array{Float64, 1},
    h    :: Float64,
    b_ML :: Float64,
    FLDO :: Integer,
) 
    # Convective adjustment
    if_unstable = bs .> b_ML
    if any(if_unstable)
        #println("UNSTABLE!!!!!!!!!!!!!!!!!!!!")
        #println(length(if_unstable))
        conv_b_new = b_ML
        # Determine the least layers to be mixed
        bottom_layer_to_convect = 1
        for i = length(if_unstable):-1:1
            println(i, "; z = ", zs[i])
            if if_unstable[i] == true
                bottom_layer_to_convect = i
                break
            end
        end
        #println("bottom_layer_to_convect = ", bottom_layer_to_convect)
        # Determine how many extra layers should be mixed
        for i = bottom_layer_to_convect:length(bs) - 1
            conv_b_new = getIntegratedBuoyancy(
                zs       =  zs,
                bs       =  bs,
                b_ML     =  b_ML,
                h        =  h,
                target_z =  zs[i+1]
            ) / (-zs[i+1])
           
            if conv_b_new < bs[i+1]
                bottom_layer_to_convect = i+1
            else
                break
            end
        end

        FLDO = bottom_layer_to_convect + 1
        println("Convective Adjustment! new_FLDO: ", FLDO, "; update b from ", b_ML, " to ", conv_b_new)
        b_ML = conv_b_new
        bs[1:FLDO-1] .= b_ML
    end
    #bs[:] .= 0    
    return h, b_ML, FLDO
end
