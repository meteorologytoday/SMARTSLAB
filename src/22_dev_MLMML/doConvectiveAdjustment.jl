function convective_adjustment!(oc :: OceanColumn) 
    # Convective adjustment
    if_unstable = oc.b_ML .< oc.bs[oc.FLDO]
    if any(if_unstable)
        conv_b_new = oc.b_ML
        # Determine the least layers to be mixed
        bottom_layer_to_convect = 1
        for i = length(if_unstable):-1:1
            if if_unstable[i] == true
                bottom_layer_to_convect = i
                break
            end
        end

        # Determine how many extra layers should be mixed
        for i = bottom_layer_to_convect:length(oc.bs) - 1
            conv_b_new = getIntegratedBuoyancy(
                zs       =  oc.zs,
                bs       =  oc.bs,
                b_ML     =  oc.b_ML,
                h        =  oc.h,
                target_z =  oc.zs[i+1]
            ) / (-oc.zs[i+1])
           
            if conv_b_new < oc.bs[i+1]
                bottom_layer_to_convect = i+1
            else
                break
            end
        end

        println("Convective Adjustment! new_FLDO: ", new_FLDO, "; update b from ", new_b_ML, " to ", conv_b_new)

        
        oc.b_ML = conv_b_new
        oc.FLDO = bottom_layer_to_convect + 1
        oc.bs[1:oc.FLDO-1] .= oc.b_ML
    end
    
end
